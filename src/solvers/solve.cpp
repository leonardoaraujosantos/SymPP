#include <sympp/solvers/solve.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/boolean.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

std::vector<Expr> solve(const Expr& expr, const Expr& var) {
    Poly p(expr, var);
    return p.roots();
}

std::vector<Expr> solve(const Expr& lhs, const Expr& rhs, const Expr& var) {
    return solve(lhs - rhs, var);
}

Matrix linsolve(const Matrix& A, const Matrix& b) {
    if (!A.is_square()) {
        throw std::invalid_argument("linsolve: A must be square");
    }
    if (b.rows() != A.rows()) {
        throw std::invalid_argument("linsolve: b dimension mismatch");
    }
    return A.inverse() * b;
}

// ----- solveset --------------------------------------------------------------

SetPtr solveset(const Expr& expr, const Expr& var) {
    return solveset(expr, var, reals());
}

namespace {

// Detect simple trig zero patterns and emit ImageSet representing the
// full periodic solution. Returns nullopt when expr isn't a recognized
// pattern.
//   sin(var)    → {n·π : n ∈ ℤ}
//   cos(var)    → {n·π + π/2 : n ∈ ℤ}
//   tan(var)    → {n·π : n ∈ ℤ}
//   sin(a·var)  → {n·π/a : n ∈ ℤ}
//   cos(a·var)  → {(2n+1)·π/(2a) : n ∈ ℤ}
[[nodiscard]] std::optional<SetPtr> trig_solveset(const Expr& expr,
                                                    const Expr& var) {
    if (expr->type_id() != TypeId::Function) return std::nullopt;
    const auto& fn = static_cast<const Function&>(*expr);
    if (fn.args().size() != 1) return std::nullopt;
    const Expr& inner = fn.args()[0];
    // Extract scaling: inner = a · var with a free of var.
    Expr a;
    if (inner == var) {
        a = S::One();
    } else if (inner->type_id() == TypeId::Mul) {
        std::vector<Expr> rest;
        bool found = false;
        for (const auto& f : inner->args()) {
            if (f == var) {
                if (found) return std::nullopt;
                found = true;
            } else if (has(f, var)) {
                return std::nullopt;
            } else {
                rest.push_back(f);
            }
        }
        if (!found) return std::nullopt;
        a = mul(rest);
    } else {
        return std::nullopt;
    }
    auto n = symbol("__n");
    Expr image;
    switch (fn.function_id()) {
        case FunctionId::Sin:
        case FunctionId::Tan:
            image = mul(n, S::Pi()) / a;
            break;
        case FunctionId::Cos:
            image = (mul(integer(2), n) + integer(1)) * S::Pi()
                    / (integer(2) * a);
            break;
        default:
            return std::nullopt;
    }
    return image_set(n, image, integers());
}

}  // namespace

SetPtr solveset(const Expr& expr, const Expr& var, const SetPtr& domain) {
    // Trig pattern: emit ImageSet for the full periodic solution.
    if (auto trig = trig_solveset(expr, var); trig) {
        // Filter against domain (only when domain is something other
        // than reals — the periodic solution is real-valued).
        if (domain->kind() == SetKind::Reals) return *trig;
        return set_intersection(*trig, domain);
    }
    auto roots = solve(expr, var);
    if (roots.empty()) return empty_set();
    // Filter against the domain when membership is decidable.
    std::vector<Expr> kept;
    kept.reserve(roots.size());
    for (auto& r : roots) {
        auto m = domain->contains(r);
        // Keep when membership is true OR unknown (be inclusive on
        // symbolic results).
        if (m == std::optional<bool>{false}) continue;
        kept.push_back(std::move(r));
    }
    return finite_set(std::move(kept));
}

// ----- nsolve --------------------------------------------------------------

namespace {

// Convert an Expr to a double via recursive numeric evaluation. Walks
// Add, Mul, Pow, and Function nodes, evalfs their leaves, and combines
// the results in double. Returns nullopt for symbolic-only inputs.
[[nodiscard]] std::optional<double> recursive_evalf_to_double(
    const Expr& e, int dps);

[[nodiscard]] std::optional<double> expr_to_double(const Expr& e, int dps) {
    Expr v = evalf(e, dps);
    if (v->type_id() == TypeId::Float) {
        try { return std::stod(v->str()); }
        catch (...) { return std::nullopt; }
    }
    return recursive_evalf_to_double(e, dps);
}

[[nodiscard]] std::optional<double> recursive_evalf_to_double(
    const Expr& e, int dps) {
    if (!e) return std::nullopt;
    switch (e->type_id()) {
        case TypeId::Integer:
        case TypeId::Rational:
        case TypeId::Float:
        case TypeId::NumberSymbol:
        case TypeId::RootOf: {
            Expr v = evalf(e, dps);
            if (v->type_id() != TypeId::Float) return std::nullopt;
            try { return std::stod(v->str()); }
            catch (...) { return std::nullopt; }
        }
        case TypeId::Add: {
            double sum = 0.0;
            for (const auto& a : e->args()) {
                auto v = recursive_evalf_to_double(a, dps);
                if (!v) return std::nullopt;
                sum += *v;
            }
            return sum;
        }
        case TypeId::Mul: {
            double prod = 1.0;
            for (const auto& a : e->args()) {
                auto v = recursive_evalf_to_double(a, dps);
                if (!v) return std::nullopt;
                prod *= *v;
            }
            return prod;
        }
        case TypeId::Pow: {
            auto b = recursive_evalf_to_double(e->args()[0], dps);
            auto x = recursive_evalf_to_double(e->args()[1], dps);
            if (!b || !x) return std::nullopt;
            return std::pow(*b, *x);
        }
        case TypeId::Function: {
            // Use evalf on the wrapped arg, then dispatch by name.
            // Simpler: substitute Float values into args and let evalf
            // handle the function. But that requires reconstruction;
            // use mpfr directly when possible.
            const auto& fn = static_cast<const Function&>(*e);
            if (fn.args().size() != 1) return std::nullopt;
            auto arg = recursive_evalf_to_double(fn.args()[0], dps);
            if (!arg) return std::nullopt;
            switch (fn.function_id()) {
                case FunctionId::Sin: return std::sin(*arg);
                case FunctionId::Cos: return std::cos(*arg);
                case FunctionId::Tan: return std::tan(*arg);
                case FunctionId::Exp: return std::exp(*arg);
                case FunctionId::Log: return std::log(*arg);
                case FunctionId::Abs: return std::fabs(*arg);
                default: return std::nullopt;
            }
        }
        default:
            return std::nullopt;
    }
}

}  // namespace

Expr nsolve(const Expr& expr, const Expr& var, const Expr& x0, int dps) {
    if (dps < 6) dps = 6;
    const int work_dps = dps + 10;
    const mpfr_prec_t work_prec = dps_to_prec(work_dps);
    Expr fprime = diff(expr, var);

    // Initial guess as MPFR.
    mpfr_t x_mpfr, step_mpfr, fv_mpfr, fpv_mpfr, tol_mpfr, abs_step;
    mpfr_init2(x_mpfr, work_prec);
    mpfr_init2(step_mpfr, work_prec);
    mpfr_init2(fv_mpfr, work_prec);
    mpfr_init2(fpv_mpfr, work_prec);
    mpfr_init2(tol_mpfr, work_prec);
    mpfr_init2(abs_step, work_prec);
    mpfr_set_d(tol_mpfr, 10.0, MPFR_RNDN);
    mpfr_pow_si(tol_mpfr, tol_mpfr, -dps, MPFR_RNDN);

    Expr x0_eval = evalf(x0, work_dps);
    if (x0_eval->type_id() != TypeId::Float) {
        // Try numeric coercion via str round-trip.
        try {
            mpfr_set_str(x_mpfr, x0->str().c_str(), 10, MPFR_RNDN);
        } catch (...) {
            mpfr_clear(x_mpfr); mpfr_clear(step_mpfr); mpfr_clear(fv_mpfr);
            mpfr_clear(fpv_mpfr); mpfr_clear(tol_mpfr); mpfr_clear(abs_step);
            throw std::runtime_error("nsolve: initial guess must be numeric");
        }
    } else {
        mpfr_set_str(x_mpfr, x0_eval->str().c_str(), 10, MPFR_RNDN);
    }

    const int max_iters = 200;
    bool converged = false;
    for (int i = 0; i < max_iters; ++i) {
        // Build a Float Expr from x_mpfr and substitute.
        char buf[256];
        mpfr_snprintf(buf, sizeof(buf), "%.*Re", work_dps, x_mpfr);
        Expr xe = make<Float>(std::string(buf), work_dps);
        Expr fv = evalf(subs(expr, var, xe), work_dps);
        Expr fpv = evalf(subs(fprime, var, xe), work_dps);
        if (fv->type_id() != TypeId::Float
            || fpv->type_id() != TypeId::Float) {
            mpfr_clear(x_mpfr); mpfr_clear(step_mpfr); mpfr_clear(fv_mpfr);
            mpfr_clear(fpv_mpfr); mpfr_clear(tol_mpfr); mpfr_clear(abs_step);
            throw std::runtime_error("nsolve: non-numeric step");
        }
        mpfr_set_str(fv_mpfr, fv->str().c_str(), 10, MPFR_RNDN);
        mpfr_set_str(fpv_mpfr, fpv->str().c_str(), 10, MPFR_RNDN);
        // step = fv / fpv
        if (mpfr_zero_p(fpv_mpfr)) {
            mpfr_clear(x_mpfr); mpfr_clear(step_mpfr); mpfr_clear(fv_mpfr);
            mpfr_clear(fpv_mpfr); mpfr_clear(tol_mpfr); mpfr_clear(abs_step);
            throw std::runtime_error("nsolve: derivative vanished");
        }
        mpfr_div(step_mpfr, fv_mpfr, fpv_mpfr, MPFR_RNDN);
        mpfr_sub(x_mpfr, x_mpfr, step_mpfr, MPFR_RNDN);
        mpfr_abs(abs_step, step_mpfr, MPFR_RNDN);
        if (mpfr_cmp(abs_step, tol_mpfr) < 0) { converged = true; break; }
    }
    if (!converged) {
        mpfr_clear(x_mpfr); mpfr_clear(step_mpfr); mpfr_clear(fv_mpfr);
        mpfr_clear(fpv_mpfr); mpfr_clear(tol_mpfr); mpfr_clear(abs_step);
        throw std::runtime_error("nsolve: failed to converge");
    }

    Expr result = make<Float>(static_cast<mpfr_srcptr>(x_mpfr), dps);
    mpfr_clear(x_mpfr); mpfr_clear(step_mpfr); mpfr_clear(fv_mpfr);
    mpfr_clear(fpv_mpfr); mpfr_clear(tol_mpfr); mpfr_clear(abs_step);
    return result;
}

// ----- solve_univariate_inequality -----------------------------------------

namespace {

[[nodiscard]] std::optional<double> sample_at(const Expr& f, const Expr& var,
                                                double x_val, int dps) {
    Expr xe = make<Float>(std::to_string(x_val), dps);
    Expr v = evalf(subs(f, var, xe), dps);
    return expr_to_double(v, dps);
}

}  // namespace

SetPtr solve_univariate_inequality(const Expr& lhs, const Expr& rel_op_str,
                                    const Expr& rhs, const Expr& var) {
    Expr f = lhs - rhs;
    auto roots = solve(f, var);
    // Filter to numeric roots; sort.
    std::vector<double> numeric_roots;
    for (const auto& r : roots) {
        auto d = expr_to_double(r, 30);
        if (d) numeric_roots.push_back(*d);
    }
    std::sort(numeric_roots.begin(), numeric_roots.end());
    numeric_roots.erase(std::unique(numeric_roots.begin(), numeric_roots.end()),
                         numeric_roots.end());

    // Decode the relation: rel_op_str is a Symbol with one of these names.
    std::string op = rel_op_str->str();
    bool want_pos = (op == ">" || op == ">=");
    bool want_neg = (op == "<" || op == "<=");
    bool include_eq = (op == ">=" || op == "<=");
    bool not_equal = (op == "!=");

    // Build the disjoint sample intervals.
    const double huge = 1e30;
    std::vector<double> bps;
    bps.push_back(-huge);
    for (auto r : numeric_roots) bps.push_back(r);
    bps.push_back(huge);

    SetPtr result = empty_set();
    for (std::size_t i = 0; i + 1 < bps.size(); ++i) {
        double mid = (bps[i] + bps[i + 1]) / 2.0;
        auto v = sample_at(f, var, mid, 30);
        if (!v) continue;
        bool keep = false;
        if (not_equal) keep = true;
        else if (want_pos && *v > 0) keep = true;
        else if (want_neg && *v < 0) keep = true;
        if (keep) {
            // Translate bps[i], bps[i+1] back to Exprs. -huge / +huge
            // act as ±∞ proxies; we still emit numeric Floats since we
            // don't have an Infinity singleton.
            Expr lo = (i == 0) ? make<Float>(std::to_string(-huge), 30)
                                : make<Float>(std::to_string(bps[i]), 30);
            Expr hi = (i + 1 == bps.size() - 1)
                          ? make<Float>(std::to_string(huge), 30)
                          : make<Float>(std::to_string(bps[i + 1]), 30);
            // Endpoint inclusion: at internal boundaries the boundary
            // is a root → strict inequality excludes, ≤/≥ includes.
            // For ±∞ proxies we always exclude.
            bool lo_open = (i == 0) || !include_eq;
            bool hi_open = (i + 1 == bps.size() - 1) || !include_eq;
            result = set_union(result,
                                interval(lo, hi, lo_open, hi_open));
        }
    }
    // For !=, the "kept" intervals collectively form Reals \ {roots}.
    if (not_equal && result->kind() != SetKind::Empty) {
        // Collapse to Reals \ FiniteSet(roots).
        std::vector<Expr> root_exprs;
        for (auto r : numeric_roots) {
            root_exprs.push_back(make<Float>(std::to_string(r), 30));
        }
        return set_complement(reals(), finite_set(std::move(root_exprs)));
    }
    return result;
}

// ----- rsolve --------------------------------------------------------------

Expr rsolve(const std::vector<Expr>& coeffs, const Expr& n) {
    if (coeffs.size() < 2) {
        throw std::invalid_argument("rsolve: need at least order 1 recurrence");
    }
    auto x = symbol("__rsolve_x");
    Poly p(std::vector<Expr>(coeffs.begin(), coeffs.end()), x);
    auto roots = p.roots();
    if (roots.empty()) {
        throw std::runtime_error(
            "rsolve: characteristic polynomial has no closed-form roots");
    }
    // Build Σ C_i * r_i^n. Track multiplicities for repeated roots:
    //   mult-k root r contributes (C_0 + C_1 n + ... + C_{k-1} n^{k-1}) r^n.
    std::vector<std::pair<Expr, std::size_t>> distinct;
    for (const auto& r : roots) {
        bool merged = false;
        for (auto& [val, count] : distinct) {
            if (val == r) {
                ++count;
                merged = true;
                break;
            }
        }
        if (!merged) distinct.emplace_back(r, 1);
    }
    Expr result = S::Zero();
    int next_c = 0;
    for (const auto& [root, mult] : distinct) {
        Expr poly_coef = S::Zero();
        for (std::size_t k = 0; k < mult; ++k) {
            auto C = symbol("__C" + std::to_string(next_c++));
            poly_coef = poly_coef + C * pow(n, integer(static_cast<long>(k)));
        }
        result = result + poly_coef * pow(root, n);
    }
    return result;
}

// ----- nonlinsolve --------------------------------------------------------

std::vector<std::vector<Expr>> nonlinsolve(
    const std::vector<Expr>& eqs, const std::vector<Expr>& vars) {
    if (eqs.size() != 2 || vars.size() != 2) {
        throw std::invalid_argument(
            "nonlinsolve: minimal implementation handles 2 eqs in 2 vars");
    }
    const Expr& x = vars[0];
    const Expr& y = vars[1];
    // Eliminate y via resultant in y; the resultant is a polynomial in x.
    Expr res = resultant(eqs[0], eqs[1], y);
    Expr res_simp = simplify(res);
    auto x_vals = solve(res_simp, x);
    std::vector<std::vector<Expr>> out;
    for (auto& xv : x_vals) {
        // Substitute into eqs[0] and solve for y.
        Expr eq_y = simplify(subs(eqs[0], x, xv));
        auto y_vals = solve(eq_y, y);
        for (auto& yv : y_vals) {
            // Verify the candidate (xv, yv) actually satisfies eqs[1] —
            // resultant elimination introduces extraneous solutions when
            // either equation has higher degree. Symbolic simplify can't
            // always reduce nested radicals, so fall back to a high-
            // precision numeric check.
            Expr check = simplify(subs(subs(eqs[1], x, xv), y, yv));
            if (!(check == S::Zero())) {
                auto v = expr_to_double(check, 30);
                if (!v || std::fabs(*v) > 1e-9) continue;
            }
            // De-duplicate.
            bool dup = false;
            for (const auto& prev : out) {
                if (prev[0] == xv && prev[1] == yv) { dup = true; break; }
            }
            if (!dup) out.push_back({xv, yv});
        }
    }
    return out;
}

// ----- linear_diophantine --------------------------------------------------

namespace {

// Extended Euclidean: gcd(a, b) = a*x + b*y. Operates on mpz.
[[nodiscard]] mpz_class ext_gcd_mpz(const mpz_class& a, const mpz_class& b,
                                     mpz_class& x, mpz_class& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    mpz_class x1, y1;
    mpz_class q = a / b;
    mpz_class g = ext_gcd_mpz(b, a - q * b, x1, y1);
    x = y1;
    y = x1 - q * y1;
    return g;
}

[[nodiscard]] std::optional<mpz_class> as_mpz(const Expr& e) {
    if (e->type_id() == TypeId::Integer) {
        return static_cast<const Integer&>(*e).value();
    }
    return std::nullopt;
}

}  // namespace

// ----- reduce_inequalities --------------------------------------------------

namespace {

[[nodiscard]] SetPtr reduce_single_rel(const Expr& rel, const Expr& var) {
    if (rel->type_id() != TypeId::Relational) {
        // Pass through; treat as Reals (assume vacuously true).
        return reals();
    }
    const auto& r = static_cast<const Relational&>(*rel);
    std::string op_name;
    switch (r.kind()) {
        case RelKind::Lt: op_name = "<"; break;
        case RelKind::Le: op_name = "<="; break;
        case RelKind::Gt: op_name = ">"; break;
        case RelKind::Ge: op_name = ">="; break;
        case RelKind::Ne: op_name = "!="; break;
        case RelKind::Eq:
            // Equality reduces to FiniteSet of solutions.
            return solveset(r.lhs() - r.rhs(), var);
    }
    auto op_sym = symbol(op_name);
    return solve_univariate_inequality(r.lhs(), op_sym, r.rhs(), var);
}

}  // namespace

SetPtr reduce_inequalities(const Expr& rel, const Expr& var) {
    return reduce_single_rel(rel, var);
}

SetPtr reduce_inequalities(const std::vector<Expr>& rels, const Expr& var,
                             bool conjunction) {
    if (rels.empty()) return conjunction ? reals() : empty_set();
    SetPtr acc = reduce_single_rel(rels[0], var);
    for (std::size_t i = 1; i < rels.size(); ++i) {
        SetPtr next = reduce_single_rel(rels[i], var);
        acc = conjunction ? set_intersection(acc, next)
                          : set_union(acc, next);
    }
    return acc;
}

// ----- Pythagorean Diophantine ---------------------------------------------

std::vector<std::array<Expr, 3>> pythagorean_triples(long max_z) {
    std::vector<std::array<Expr, 3>> out;
    if (max_z < 5) return out;  // smallest is (3, 4, 5)
    for (long m = 2; m * m <= max_z; ++m) {
        for (long n = 1; n < m; ++n) {
            // gcd(m, n) == 1 and opposite parity → primitive.
            mpz_class mz(m), nz(n);
            mpz_class g;
            mpz_gcd(g.get_mpz_t(), mz.get_mpz_t(), nz.get_mpz_t());
            if (g != 1) continue;
            if ((m % 2) == (n % 2)) continue;
            long x_val = m * m - n * n;
            long y_val = 2 * m * n;
            long z_val = m * m + n * n;
            if (z_val > max_z) break;
            out.push_back({integer(x_val), integer(y_val), integer(z_val)});
        }
    }
    return out;
}

// ----- Gröbner basis (Buchberger) and multivariate nonlinsolve --------------
//
// Minimal Buchberger over ℚ with lex order. Polynomials are represented
// as vectors of (coefficient, exponent_tuple) terms where the tuple is
// indexed by variable position. We build dense vectors keyed on
// monomials.

namespace {

// Multivariate monomial: vector of long exponents, one per variable.
struct Monomial {
    std::vector<long> e;
    // Lex ordering: compare exponents left-to-right. Returns -1/0/+1.
    [[nodiscard]] int cmp(const Monomial& o) const {
        for (std::size_t i = 0; i < e.size(); ++i) {
            if (e[i] != o.e[i]) return e[i] < o.e[i] ? -1 : 1;
        }
        return 0;
    }
    [[nodiscard]] bool operator==(const Monomial& o) const { return cmp(o) == 0; }
    [[nodiscard]] bool operator<(const Monomial& o) const { return cmp(o) < 0; }
};

// Multivariate polynomial as map monomial → mpq coefficient.
struct MPoly {
    std::vector<std::pair<Monomial, mpq_class>> terms;  // sorted by monomial desc

    [[nodiscard]] bool is_zero() const { return terms.empty(); }
    [[nodiscard]] const Monomial& leading_mono() const { return terms.front().first; }
    [[nodiscard]] const mpq_class& leading_coef() const { return terms.front().second; }

    void canonicalize() {
        // Sort descending by monomial.
        std::sort(terms.begin(), terms.end(),
                   [](const auto& a, const auto& b) { return b.first < a.first; });
        // Combine equal monomials.
        std::vector<std::pair<Monomial, mpq_class>> out;
        for (auto& t : terms) {
            if (!out.empty() && out.back().first == t.first) {
                out.back().second += t.second;
                if (out.back().second == 0) out.pop_back();
            } else if (!(t.second == 0)) {
                out.push_back(std::move(t));
            }
        }
        terms = std::move(out);
    }
};

// Convert an Expr polynomial in `vars` to an MPoly via expand + walk.
[[nodiscard]] MPoly expr_to_mpoly(const Expr& expr,
                                    const std::vector<Expr>& vars) {
    MPoly p;
    Expr e = expand(expr);
    auto add_term = [&](const mpq_class& coef,
                        std::vector<long> exps) {
        p.terms.push_back({Monomial{std::move(exps)}, coef});
    };
    auto handle_term = [&](const Expr& term) -> bool {
        // term must factor as (numeric coef) * (Π vars[i]^exp_i).
        std::vector<long> exps(vars.size(), 0);
        mpq_class coef = 1;
        auto handle_factor = [&](const Expr& f) -> bool {
            if (f->type_id() == TypeId::Integer) {
                coef *= mpq_class(static_cast<const Integer&>(*f).value());
                return true;
            }
            if (f->type_id() == TypeId::Rational) {
                coef *= static_cast<const Rational&>(*f).value();
                return true;
            }
            // Match var or var^k.
            for (std::size_t i = 0; i < vars.size(); ++i) {
                if (f == vars[i]) { exps[i] += 1; return true; }
                if (f->type_id() == TypeId::Pow && f->args()[0] == vars[i]
                    && f->args()[1]->type_id() == TypeId::Integer) {
                    long k = static_cast<const Integer&>(*f->args()[1]).to_long();
                    if (k < 0) return false;
                    exps[i] += k;
                    return true;
                }
            }
            return false;
        };
        if (term->type_id() == TypeId::Mul) {
            for (const auto& f : term->args()) {
                if (!handle_factor(f)) return false;
            }
        } else {
            if (!handle_factor(term)) return false;
        }
        if (coef == 0) return true;
        add_term(coef, std::move(exps));
        return true;
    };
    if (e->type_id() == TypeId::Add) {
        for (const auto& t : e->args()) {
            if (!handle_term(t)) {
                p.terms.clear();  // malformed; signal failure below
                return p;
            }
        }
    } else {
        if (!handle_term(e)) p.terms.clear();
    }
    p.canonicalize();
    return p;
}

[[nodiscard]] Expr mpoly_to_expr(const MPoly& p,
                                   const std::vector<Expr>& vars) {
    if (p.terms.empty()) return S::Zero();
    std::vector<Expr> term_exprs;
    for (const auto& [mono, coef] : p.terms) {
        std::vector<Expr> factors;
        // coefficient
        mpq_class cc = coef;
        cc.canonicalize();
        if (cc.get_den() == 1) {
            factors.push_back(make<Integer>(cc.get_num()));
        } else {
            factors.push_back(make<Rational>(cc));
        }
        for (std::size_t i = 0; i < vars.size(); ++i) {
            if (mono.e[i] == 0) continue;
            if (mono.e[i] == 1) factors.push_back(vars[i]);
            else factors.push_back(pow(vars[i], integer(mono.e[i])));
        }
        term_exprs.push_back(mul(factors));
    }
    return add(term_exprs);
}

// MPoly arithmetic helpers.
[[nodiscard]] MPoly mpoly_sub(const MPoly& a, const MPoly& b) {
    MPoly out = a;
    for (const auto& t : b.terms) {
        out.terms.push_back({t.first, -t.second});
    }
    out.canonicalize();
    return out;
}

[[nodiscard]] MPoly mpoly_scale_mono(const MPoly& a, const mpq_class& coef,
                                        const Monomial& mono) {
    MPoly out;
    out.terms.reserve(a.terms.size());
    for (const auto& t : a.terms) {
        Monomial m;
        m.e.resize(mono.e.size());
        for (std::size_t i = 0; i < mono.e.size(); ++i) {
            m.e[i] = t.first.e[i] + mono.e[i];
        }
        out.terms.push_back({std::move(m), t.second * coef});
    }
    return out;  // already sorted (multiplication preserves order)
}

[[nodiscard]] bool mono_divides(const Monomial& a, const Monomial& b) {
    for (std::size_t i = 0; i < a.e.size(); ++i) {
        if (a.e[i] > b.e[i]) return false;
    }
    return true;
}

[[nodiscard]] Monomial mono_sub(const Monomial& a, const Monomial& b) {
    Monomial m;
    m.e.resize(a.e.size());
    for (std::size_t i = 0; i < a.e.size(); ++i) m.e[i] = a.e[i] - b.e[i];
    return m;
}

[[nodiscard]] Monomial mono_lcm(const Monomial& a, const Monomial& b) {
    Monomial m;
    m.e.resize(a.e.size());
    for (std::size_t i = 0; i < a.e.size(); ++i) {
        m.e[i] = std::max(a.e[i], b.e[i]);
    }
    return m;
}

// Multivariate polynomial reduction: reduce f modulo a list G of
// polynomials. Returns the remainder.
[[nodiscard]] MPoly mpoly_reduce(MPoly f, const std::vector<MPoly>& G) {
    MPoly r;
    while (!f.is_zero()) {
        bool reduced = false;
        for (const auto& g : G) {
            if (g.is_zero()) continue;
            if (mono_divides(g.leading_mono(), f.leading_mono())) {
                Monomial diff = mono_sub(f.leading_mono(), g.leading_mono());
                mpq_class coef = f.leading_coef() / g.leading_coef();
                MPoly scaled = mpoly_scale_mono(g, coef, diff);
                f = mpoly_sub(f, scaled);
                reduced = true;
                break;
            }
        }
        if (!reduced) {
            r.terms.push_back(f.terms.front());
            f.terms.erase(f.terms.begin());
        }
    }
    return r;
}

// S-polynomial of f, g.
[[nodiscard]] MPoly s_poly(const MPoly& f, const MPoly& g) {
    Monomial l = mono_lcm(f.leading_mono(), g.leading_mono());
    Monomial diff_f = mono_sub(l, f.leading_mono());
    Monomial diff_g = mono_sub(l, g.leading_mono());
    MPoly fterm = mpoly_scale_mono(f, mpq_class(1) / f.leading_coef(), diff_f);
    MPoly gterm = mpoly_scale_mono(g, mpq_class(1) / g.leading_coef(), diff_g);
    return mpoly_sub(fterm, gterm);
}

}  // namespace

std::vector<Expr> groebner(const std::vector<Expr>& eqs,
                             const std::vector<Expr>& vars) {
    std::vector<MPoly> G;
    G.reserve(eqs.size());
    for (const auto& e : eqs) {
        auto p = expr_to_mpoly(e, vars);
        if (!p.is_zero()) G.push_back(std::move(p));
    }
    // Buchberger's algorithm.
    std::size_t i = 0;
    while (i < G.size()) {
        for (std::size_t j = 0; j < i; ++j) {
            MPoly s = s_poly(G[j], G[i]);
            MPoly r = mpoly_reduce(s, G);
            if (!r.is_zero()) G.push_back(std::move(r));
        }
        ++i;
        if (G.size() > 200) break;  // safety bound
    }
    // Reduce basis: drop polynomials whose leading mono is divisible
    // by another's, normalize each to monic.
    std::vector<MPoly> reduced;
    for (std::size_t k = 0; k < G.size(); ++k) {
        bool drop = false;
        for (std::size_t l = 0; l < G.size(); ++l) {
            if (l == k || G[l].is_zero()) continue;
            if (mono_divides(G[l].leading_mono(), G[k].leading_mono())
                && !(G[l].leading_mono() == G[k].leading_mono())) {
                drop = true;
                break;
            }
        }
        if (!drop && !G[k].is_zero()) {
            MPoly p = G[k];
            mpq_class lc = p.leading_coef();
            for (auto& t : p.terms) t.second /= lc;
            reduced.push_back(std::move(p));
        }
    }
    // Final pass: reduce each polynomial against the others.
    for (std::size_t k = 0; k < reduced.size(); ++k) {
        std::vector<MPoly> others;
        for (std::size_t l = 0; l < reduced.size(); ++l) {
            if (l != k) others.push_back(reduced[l]);
        }
        MPoly red = mpoly_reduce(reduced[k], others);
        // Re-normalize.
        if (!red.is_zero()) {
            mpq_class lc = red.leading_coef();
            for (auto& t : red.terms) t.second /= lc;
        }
        reduced[k] = std::move(red);
    }
    std::vector<Expr> out;
    out.reserve(reduced.size());
    for (auto& p : reduced) {
        if (!p.is_zero()) out.push_back(mpoly_to_expr(p, vars));
    }
    return out;
}

std::vector<std::vector<Expr>> nonlinsolve_groebner(
    const std::vector<Expr>& eqs, const std::vector<Expr>& vars) {
    auto basis = groebner(eqs, vars);
    if (basis.empty()) return {};
    // Lex order Gröbner basis with vars = [x, y, z, ...] is triangular
    // with the last variable's univariate polynomial *last* in the
    // elimination ideal. Solve in reverse order: z first, then y given
    // z, then x given y, z. Build solutions as a vector indexed by the
    // original `vars` order.
    const std::size_t n = vars.size();
    std::function<std::vector<std::vector<Expr>>(std::vector<Expr>, long)>
        rec = [&](std::vector<Expr> partial,
                  long var_idx) -> std::vector<std::vector<Expr>> {
        if (var_idx < 0) {
            return {partial};
        }
        // Substitute already-known values (for indices > var_idx) into
        // each basis element; pick the first that becomes univariate
        // in vars[var_idx].
        for (const auto& g : basis) {
            Expr g_sub = g;
            for (std::size_t k = static_cast<std::size_t>(var_idx + 1);
                 k < n; ++k) {
                g_sub = subs(g_sub, vars[k], partial[k]);
            }
            g_sub = expand(simplify(g_sub));
            // Must contain vars[var_idx] and no earlier variable
            // vars[k] for k < var_idx.
            if (!has(g_sub, vars[static_cast<std::size_t>(var_idx)])) continue;
            bool ok = true;
            for (long k = 0; k < var_idx; ++k) {
                if (has(g_sub, vars[static_cast<std::size_t>(k)])) {
                    ok = false; break;
                }
            }
            if (!ok) continue;
            auto roots = solve(g_sub,
                                vars[static_cast<std::size_t>(var_idx)]);
            std::vector<std::vector<Expr>> out;
            for (auto& r : roots) {
                auto next = partial;
                next[static_cast<std::size_t>(var_idx)] = r;
                auto sub_sols = rec(next, var_idx - 1);
                for (auto& s : sub_sols) out.push_back(std::move(s));
            }
            return out;
        }
        return {};  // no triangular polynomial found
    };
    std::vector<Expr> partial(n, S::Zero());
    return rec(partial, static_cast<long>(n) - 1);
}

// ----- linear_diophantine --------------------------------------------------

std::optional<std::pair<Expr, Expr>>
linear_diophantine(const Expr& a, const Expr& b, const Expr& c) {
    auto a_z = as_mpz(a);
    auto b_z = as_mpz(b);
    auto c_z = as_mpz(c);
    if (!a_z || !b_z || !c_z) return std::nullopt;
    mpz_class x, y;
    mpz_class g = ext_gcd_mpz(*a_z, *b_z, x, y);
    // Reduce sign: we want g > 0.
    if (sgn(g) < 0) { g = -g; x = -x; y = -y; }
    if ((*c_z) % g != 0) return std::nullopt;
    mpz_class scale = (*c_z) / g;
    x *= scale;
    y *= scale;
    return std::pair<Expr, Expr>{make<Integer>(x), make<Integer>(y)};
}

}  // namespace sympp
