#include "gosper.hpp"

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp::detail {

namespace {

// ----- rational coefficient polynomials (ascending; [0] = constant) ----------
using QPoly = std::vector<mpq_class>;

void trim(QPoly& a) {
    while (a.size() > 1 && a.back() == 0) a.pop_back();
}

[[nodiscard]] std::optional<mpq_class> expr_to_q(const Expr& e) {
    if (e->type_id() == TypeId::Integer) {
        return mpq_class(static_cast<const Integer&>(*e).value());
    }
    if (e->type_id() == TypeId::Rational) {
        return static_cast<const Rational&>(*e).value();
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<QPoly> poly_to_q(const Poly& p) {
    QPoly out;
    out.reserve(p.coeffs().size());
    for (const auto& c : p.coeffs()) {
        auto q = expr_to_q(c);
        if (!q) return std::nullopt;
        out.push_back(*q);
    }
    if (out.empty()) out.push_back(mpq_class(0));
    trim(out);
    return out;
}

[[nodiscard]] QPoly q_mul(const QPoly& a, const QPoly& b) {
    QPoly r(a.size() + b.size() - 1, mpq_class(0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) r[i + j] += a[i] * b[j];
    }
    trim(r);
    return r;
}

[[nodiscard]] QPoly q_sub(const QPoly& a, const QPoly& b) {
    QPoly r(std::max(a.size(), b.size()), mpq_class(0));
    for (std::size_t i = 0; i < a.size(); ++i) r[i] += a[i];
    for (std::size_t i = 0; i < b.size(); ++i) r[i] -= b[i];
    trim(r);
    return r;
}

// P(k + s) for an integer shift s, via Horner over the binomial expansion of
// (k + s)^i accumulated by repeated multiplication by (k + s).
[[nodiscard]] QPoly q_shift(const QPoly& a, long s) {
    QPoly base{mpq_class(s), mpq_class(1)};  // k + s
    QPoly kpow{mpq_class(1)};                // (k + s)^0
    QPoly result{mpq_class(0)};
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < kpow.size(); ++j) {
            if (result.size() <= j) result.resize(j + 1, mpq_class(0));
            result[j] += a[i] * kpow[j];
        }
        kpow = q_mul(kpow, base);
    }
    trim(result);
    return result;
}

[[nodiscard]] Expr q_to_expr(const mpq_class& q) {
    return rational(q.get_num(), q.get_den());
}

// Solve the (possibly over-determined) exact rational system M·x = rhs. Returns
// a particular solution (free variables set to 0) when consistent, else nullopt.
[[nodiscard]] std::optional<std::vector<mpq_class>> solve_linear(
    std::vector<std::vector<mpq_class>> M, std::vector<mpq_class> rhs) {
    const std::vector<std::vector<mpq_class>> M0 = M;
    const std::vector<mpq_class> rhs0 = rhs;
    const std::size_t rows = M.size();
    const std::size_t cols = rows ? M[0].size() : 0;
    std::vector<long> pivot_row(cols, -1);
    std::size_t row = 0;
    for (std::size_t col = 0; col < cols && row < rows; ++col) {
        std::size_t sel = row;
        bool found = false;
        for (std::size_t i = row; i < rows; ++i) {
            if (M[i][col] != 0) { sel = i; found = true; break; }
        }
        if (!found) continue;
        std::swap(M[row], M[sel]);
        std::swap(rhs[row], rhs[sel]);
        mpq_class piv = M[row][col];
        for (std::size_t j = col; j < cols; ++j) M[row][j] /= piv;
        rhs[row] /= piv;
        for (std::size_t i = 0; i < rows; ++i) {
            if (i != row && M[i][col] != 0) {
                mpq_class f = M[i][col];
                for (std::size_t j = col; j < cols; ++j) M[i][j] -= f * M[row][j];
                rhs[i] -= f * rhs[row];
            }
        }
        pivot_row[col] = static_cast<long>(row);
        ++row;
    }
    std::vector<mpq_class> x(cols, mpq_class(0));
    for (std::size_t c = 0; c < cols; ++c) {
        if (pivot_row[c] >= 0) x[c] = rhs[static_cast<std::size_t>(pivot_row[c])];
    }
    // Verify against the original system (catches inconsistency).
    for (std::size_t i = 0; i < rows; ++i) {
        mpq_class acc(0);
        for (std::size_t j = 0; j < cols; ++j) acc += M0[i][j] * x[j];
        if (acc != rhs0[i]) return std::nullopt;
    }
    return x;
}

// ----- hypergeometric term ratio t(k+1)/t(k) ---------------------------------
// Accumulates the ratio of a single factor into the running numerator/denominator
// (both Exprs, later turned into polynomials). Returns false if the factor is not
// hypergeometric in var.
[[nodiscard]] bool factor_ratio(const Expr& f, const Expr& var, Expr& num,
                                Expr& den) {
    if (!has(f, var)) return true;  // constant factor: ratio 1
    auto shift1 = [&](const Expr& e) { return subs(e, var, add(var, integer(1))); };

    if (f == var) {  // k → (k+1)/k
        num = mul(num, add(var, integer(1)));
        den = mul(den, var);
        return true;
    }

    if (f->type_id() == TypeId::Pow) {
        const Expr& base = f->args()[0];
        const Expr& e = f->args()[1];
        // c^(a·k+b): ratio = c^a (constant).
        if (!has(base, var) && has(e, var)) {
            Expr a = diff(e, var);
            if (has(a, var)) return false;
            num = mul(num, pow(base, a));
            return true;
        }
        // base(k)^m, m a constant integer: ratio = (base(k+1)/base(k))^m.
        if (has(base, var) && !has(e, var) && e->type_id() == TypeId::Integer) {
            long m = static_cast<const Integer&>(*e).to_long();
            Expr b1 = shift1(base);
            if (m >= 0) {
                num = mul(num, pow(b1, integer(m)));
                den = mul(den, pow(base, integer(m)));
            } else {
                num = mul(num, pow(base, integer(-m)));
                den = mul(den, pow(b1, integer(-m)));
            }
            return true;
        }
        return false;
    }

    if (f->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*f);
        const FunctionId id = fn.function_id();
        const auto& A = fn.args();
        // factorial / gamma of a·k + c (a a positive integer, c var-free):
        //   factorial(a·k+c) ratio = ∏_{j=1}^{a} (a·k + c + j)
        //   gamma(a·k+c)     ratio = ∏_{j=0}^{a-1} (a·k + c + j)
        if ((id == FunctionId::Factorial || id == FunctionId::Gamma)
            && A.size() == 1) {
            Expr a = diff(A[0], var);
            if (a->type_id() != TypeId::Integer) return false;
            long av = static_cast<const Integer&>(*a).to_long();
            if (av <= 0) return false;
            Expr c = simplify(add(A[0], mul(integer(-1), mul(a, var))));
            if (has(c, var)) return false;
            long jlo = (id == FunctionId::Gamma) ? 0 : 1;
            long jhi = (id == FunctionId::Gamma) ? (av - 1) : av;
            Expr prod = integer(1);
            for (long j = jlo; j <= jhi; ++j) {
                prod = mul(prod, add(mul(a, var), add(c, integer(j))));
            }
            num = mul(num, prod);
            return true;
        }
        if (id == FunctionId::Binomial && A.size() == 2) {
            // C(n,k): (n−k)/(k+1);  C(k,m): (k+1)/(k+1−m).
            if (!has(A[0], var) && A[1] == var) {
                num = mul(num, add(A[0], mul(integer(-1), var)));
                den = mul(den, add(var, integer(1)));
                return true;
            }
            if (A[0] == var && !has(A[1], var)) {
                num = mul(num, add(var, integer(1)));
                den = mul(den, add(add(var, integer(1)), mul(integer(-1), A[1])));
                return true;
            }
            return false;
        }
        // RisingFactorial(a, k): ratio = a + k.
        if (id == FunctionId::RisingFactorial && A.size() == 2) {
            if (!has(A[0], var) && A[1] == var) {
                num = mul(num, add(A[0], var));
                return true;
            }
            return false;
        }
        return false;
    }

    // Polynomial factor in var: ratio = p(k+1)/p(k).
    try {
        Poly pf(expand(f), var);
        for (const auto& cc : pf.coeffs()) {
            if (has(cc, var)) return false;
        }
    } catch (const std::exception&) {
        return false;
    }
    num = mul(num, shift1(f));
    den = mul(den, f);
    return true;
}

[[nodiscard]] std::optional<std::pair<Poly, Poly>> hyper_ratio(const Expr& t,
                                                               const Expr& var) {
    Expr num = integer(1);
    Expr den = integer(1);
    std::vector<Expr> factors;
    if (t->type_id() == TypeId::Mul) {
        for (const auto& f : t->args()) factors.push_back(f);
    } else {
        factors.push_back(t);
    }
    for (const auto& f : factors) {
        if (!factor_ratio(f, var, num, den)) return std::nullopt;
    }
    try {
        Poly a(expand(num), var);
        Poly b(expand(den), var);
        Poly g = gcd(a, b);
        if (g.degree() >= 1) {
            a = a / g;
            b = b / g;
        }
        return std::make_pair(a, b);
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

// True if the term t has a pole at an integer point ≥ lo (hence possibly inside
// the range [lo, hi]) coming from a reciprocal polynomial factor — e.g.
// 1/(k(k−1)) at k=1. Such a sum is undefined and must not be closed by Gosper.
[[nodiscard]] bool has_pole_in_range(const Expr& t, const Expr& var,
                                     const Expr& lo) {
    Expr denom = integer(1);
    std::vector<Expr> facs;
    if (t->type_id() == TypeId::Mul) {
        for (const auto& f : t->args()) facs.push_back(f);
    } else {
        facs.push_back(t);
    }
    bool any = false;
    for (const auto& f : facs) {
        if (f->type_id() == TypeId::Pow
            && f->args()[1]->type_id() == TypeId::Integer
            && static_cast<const Integer&>(*f->args()[1]).is_negative()
            && has(f->args()[0], var)) {
            long e = -static_cast<const Integer&>(*f->args()[1]).to_long();
            denom = mul(denom, pow(f->args()[0], integer(e)));
            any = true;
        }
    }
    if (!any) return false;
    try {
        Poly pd(expand(denom), var);
        for (const auto& r : rational_roots(pd)) {
            if (r->type_id() != TypeId::Integer) continue;
            // r ≥ lo ⇔ lo − r is not strictly positive.
            if (is_positive(simplify(add(lo, mul(integer(-1), r))))
                != std::optional<bool>{true}) {
                return true;
            }
        }
    } catch (const std::exception&) {
        return true;  // can't analyze the denominator — stay safe, bail
    }
    return false;
}

// ----- Gosper–Petkovšek normal form a/b = (p(k+1)/p(k))·(q(k)/r(k)) ----------
struct PQR {
    Poly p, q, r;
};

[[nodiscard]] Poly poly_shift(const Poly& P, long s, const Expr& var) {
    return Poly(expand(subs(P.as_expr(), var, add(var, integer(s)))), var);
}

[[nodiscard]] PQR petkovsek(const Poly& a, const Poly& b, const Expr& var) {
    PQR out{Poly(integer(1), var), a, b};
    // Dispersion candidates h = β − α ≥ 0 (α a rational root of a, β of b).
    std::vector<Expr> ra = rational_roots(a);
    std::vector<Expr> rb = rational_roots(b);
    std::set<long> hs;
    for (const auto& al : ra) {
        for (const auto& be : rb) {
            Expr d = simplify(add(be, mul(integer(-1), al)));
            if (d->type_id() == TypeId::Integer
                && static_cast<const Integer&>(*d).fits_long()) {
                long h = static_cast<const Integer&>(*d).to_long();
                if (h >= 0) hs.insert(h);
            }
        }
    }
    for (long h : hs) {
        Poly g = gcd(out.q, poly_shift(out.r, h, var));
        if (g.degree() < 1) continue;
        out.q = out.q / g;
        out.r = out.r / poly_shift(g, -h, var);
        for (long i = 1; i <= h; ++i) out.p = out.p * poly_shift(g, -i, var);
    }
    return out;
}

// Solve q(k)·f(k+1) − r(k)·f(k) = p(k)·r(k) for a polynomial f, searching the
// degree. Returns f as an Expr in var, or nullopt if no polynomial solution.
[[nodiscard]] std::optional<Expr> solve_certificate(const Poly& p, const Poly& q,
                                                    const Poly& r,
                                                    const Expr& var) {
    auto pq = poly_to_q(p);
    auto qq = poly_to_q(q);
    auto rq = poly_to_q(r);
    if (!pq || !qq || !rq) return std::nullopt;
    const QPoly RHS = q_mul(*pq, *rq);
    const long degQ = static_cast<long>(qq->size()) - 1;
    const long degR = static_cast<long>(rq->size()) - 1;
    const long degRHS = static_cast<long>(RHS.size()) - 1;
    long dmax = std::max(degRHS - std::min(degQ, degR), 0L) + 4;
    if (dmax > 20) dmax = 20;

    for (long d = 0; d <= dmax; ++d) {
        std::vector<QPoly> basis(static_cast<std::size_t>(d) + 1);
        std::size_t maxdeg = RHS.size() ? RHS.size() - 1 : 0;
        for (long i = 0; i <= d; ++i) {
            QPoly ei(static_cast<std::size_t>(i) + 1, mpq_class(0));
            ei[static_cast<std::size_t>(i)] = 1;  // k^i
            // contribution of f_i: q(k)·(k+1)^i − r(k)·k^i
            QPoly contrib = q_sub(q_mul(*qq, q_shift(ei, 1)), q_mul(*rq, ei));
            maxdeg = std::max(maxdeg, contrib.size() ? contrib.size() - 1 : 0);
            basis[static_cast<std::size_t>(i)] = std::move(contrib);
        }
        const std::size_t rows = maxdeg + 1;
        const std::size_t cols = static_cast<std::size_t>(d) + 1;
        std::vector<std::vector<mpq_class>> M(rows,
                                              std::vector<mpq_class>(cols, mpq_class(0)));
        std::vector<mpq_class> rhs(rows, mpq_class(0));
        for (std::size_t i = 0; i < cols; ++i) {
            for (std::size_t m = 0; m < basis[i].size(); ++m) M[m][i] = basis[i][m];
        }
        for (std::size_t m = 0; m < RHS.size(); ++m) rhs[m] = RHS[m];

        auto sol = solve_linear(M, rhs);
        if (!sol) continue;
        Expr f = integer(0);
        bool nonzero = false;
        for (std::size_t i = 0; i < cols; ++i) {
            if ((*sol)[i] == 0) continue;
            nonzero = true;
            f = add(f, mul(q_to_expr((*sol)[i]),
                           pow(var, integer(static_cast<long>(i)))));
        }
        if (nonzero) return f;
    }
    return std::nullopt;
}

}  // namespace

std::optional<Expr> gosper_summation(const Expr& t, const Expr& var,
                                     const Expr& lo, const Expr& hi) {
    if (!t || !has(t, var)) return std::nullopt;
    if (hi->type_id() == TypeId::Infinity || lo->type_id() == TypeId::Infinity) {
        return std::nullopt;  // Gosper handles finite ranges
    }
    // Gosper works over ℚ on a term univariate in var. A summand carrying any
    // other symbol (a parameter, e.g. binomial(n, k)) would drive Petkovšek's
    // symbolic gcd / root finding into a blow-up — bail and let Zeilberger take
    // the parametric case.
    for (const auto& s : free_symbols(t)) {
        if (!(s == var)) return std::nullopt;
    }
    if (has_pole_in_range(t, var, lo)) return std::nullopt;
    auto ratio = hyper_ratio(t, var);
    if (!ratio) return std::nullopt;
    Poly a = ratio->first;
    Poly b = ratio->second;
    if (a.degree() == 0 && b.degree() == 0) return std::nullopt;  // constant ratio

    PQR pqr = petkovsek(a, b, var);
    auto fopt = solve_certificate(pqr.p, pqr.q, pqr.r, var);
    if (!fopt) return std::nullopt;

    // Rational certificate y(k) = f(k)/p(k); antidifference z(k) = y(k)·t(k).
    Expr y = mul(*fopt, pow(pqr.p.as_expr(), integer(-1)));
    Expr aE = a.as_expr();
    Expr bE = b.as_expr();
    // Verify the Gosper identity y(k+1)·a/b − y(k) − 1 = 0 (purely rational, so
    // simplify resolves it without touching the factorials in t). This is the
    // safety gate: a wrong certificate is rejected rather than returned.
    Expr check = simplify(add(
        add(mul(subs(y, var, add(var, integer(1))), mul(aE, pow(bE, integer(-1)))),
            mul(integer(-1), y)),
        integer(-1)));
    if (!(check == S::Zero())) return std::nullopt;

    // Antidifference z(k) = y(k)·t(k), simplified as a whole so that a pole of
    // the rational y(k) cancels against a zero of t(k) (e.g. y = 1/k against the
    // k = 0 term of k·k!) before the boundaries are substituted — otherwise
    // z(0) would evaluate as (1/0)·0 = nan.
    Expr z = simplify(mul(y, t));
    Expr z_hi = subs(z, var, add(hi, integer(1)));
    Expr z_lo = subs(z, var, lo);
    Expr result = simplify(add(z_hi, mul(integer(-1), z_lo)));
    // Reject any residual singularity (a genuine pole inside the range).
    if (has(result, var) || result->type_id() == TypeId::NaN
        || result->type_id() == TypeId::ComplexInfinity) {
        return std::nullopt;
    }
    return result;
}

}  // namespace sympp::detail
