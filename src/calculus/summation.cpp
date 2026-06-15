#include <sympp/calculus/summation.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/expand.hpp>
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
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/core/undefined_function.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/special.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

namespace {

// Unevaluated sum marker — an UndefinedFunction `Sum(expr, var, lo, hi)`,
// mirroring the `Integral(_, _)` marker. Returned when no closed form is
// found, so a sum is never silently dropped to its bare summand.
[[nodiscard]] Expr sum_marker(const Expr& expr, const Expr& var,
                              const Expr& lo, const Expr& hi) {
    return function_symbol("Sum")(std::vector<Expr>{expr, var, lo, hi});
}

// Telescoping sum of a rational summand c / D(var), where D is a quadratic that
// factors into two distinct linear factors whose roots differ by a nonzero
// integer. Partial fractions give c/(lead·(k−r1)(k−r2)) = c/(lead·d)·[u(k) −
// u(k+d)] with u(k) = 1/(k−r1) and d = r1 − r2 > 0, which telescopes to
//   c/(lead·d) · [ Σ_{j=0}^{d−1} u(lo+j) − Σ_{j=1}^{d} u(hi+j) ].
// Returns nullopt unless the summand has this shape. The closed form is exact
// (equivalent to SymPy's, though it may be presented differently).
[[nodiscard]] std::optional<Expr> telescope_rational(
    const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    // Split expr = c · D^{-1}: c var-free, D one or more reciprocal factors.
    Expr c = S::One();
    Expr denom;
    auto take_recip = [&](const Expr& f) -> bool {
        if (f->type_id() == TypeId::Pow && f->args()[1] == S::NegativeOne()) {
            denom = denom ? mul({denom, f->args()[0]}) : f->args()[0];
            return true;
        }
        return false;
    };
    if (expr->type_id() == TypeId::Pow) {
        if (!take_recip(expr)) return std::nullopt;
    } else if (expr->type_id() == TypeId::Mul) {
        std::vector<Expr> cf;
        for (const auto& f : expr->args()) {
            if (take_recip(f)) continue;
            if (has(f, var)) return std::nullopt;  // var in the numerator
            cf.push_back(f);
        }
        if (!cf.empty()) c = mul(std::move(cf));
    } else {
        return std::nullopt;
    }
    if (!denom || !has(denom, var)) return std::nullopt;
    // D must be a quadratic with two distinct rational roots.
    Poly p(expand(denom), var);
    if (p.degree() != 2) return std::nullopt;
    std::vector<Expr> rts = p.roots();
    if (rts.size() != 2) return std::nullopt;
    Expr r1 = rts[0], r2 = rts[1];
    auto is_rat = [](const Expr& e) {
        return e->type_id() == TypeId::Integer
               || e->type_id() == TypeId::Rational;
    };
    if (!is_rat(r1) || !is_rat(r2) || r1 == r2) return std::nullopt;
    const Expr lead = p.leading_coeff();
    // d = r1 − r2 must be a nonzero integer; orient it positive.
    Expr d_expr = simplify(r1 - r2);
    if (d_expr->type_id() != TypeId::Integer
        || !static_cast<const Integer&>(*d_expr).fits_long()) {
        return std::nullopt;
    }
    long d = static_cast<const Integer&>(*d_expr).to_long();
    if (d == 0) return std::nullopt;
    if (d < 0) {
        std::swap(r1, r2);
        d = -d;
    }
    // Guard against a pole inside the range: an integer root ≥ lo makes a
    // summand undefined. A non-integer root is never hit by an integer index.
    auto safe_root = [&](const Expr& r) {
        if (r->type_id() != TypeId::Integer) return true;
        return is_positive(simplify(lo - r)) == std::optional<bool>{true};
    };
    if (!safe_root(r1) || !safe_root(r2)) return std::nullopt;
    auto u = [&](const Expr& k) { return pow(k - r1, integer(-1)); };
    Expr first = S::Zero(), second = S::Zero();
    for (long j = 0; j < d; ++j) first = add({first, u(lo + integer(j))});
    for (long j = 1; j <= d; ++j) second = add({second, u(hi + integer(j))});
    Expr scale = pow(mul({lead, integer(d)}), integer(-1));
    return simplify(
        mul({c, scale, add({first, mul({S::NegativeOne(), second})})}));
}

// Σ_{k=lo}^{hi} P(k)·r^k for a polynomial P (degree ≥ 1) and a concrete numeric
// ratio r ≠ 1. Finds the antidifference S(k) = Q(k)·r^k, where Q is a polynomial
// of the same degree satisfying r·Q(k+1) − Q(k) = P(k); the coefficients solve
// top-down since q_j = [P_j − r·Σ_{i>j} C(i,j) q_i]/(r−1). Then the sum is
// S(hi+1) − S(lo). Returns nullopt unless the summand has this shape.
[[nodiscard]] std::optional<Expr> sum_poly_geometric(
    const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    if (expr->type_id() != TypeId::Mul) return std::nullopt;
    // Split off one geometric factor base^(c·var+d); the rest is the polynomial.
    Expr geo;
    std::vector<Expr> poly_factors;
    for (const auto& f : expr->args()) {
        if (!geo && f->type_id() == TypeId::Pow && !has(f->args()[0], var)
            && is_number(f->args()[0]) && has(f->args()[1], var)) {
            geo = f;
        } else {
            poly_factors.push_back(f);
        }
    }
    if (!geo || poly_factors.empty()) return std::nullopt;
    Expr c = diff(geo->args()[1], var);
    if (has(c, var)) return std::nullopt;  // exponent not linear in var
    Expr d = simplify(geo->args()[1] - c * var);
    if (has(d, var)) return std::nullopt;
    Expr ratio = pow(geo->args()[0], c);
    if (!is_number(ratio) || ratio == S::One()) return std::nullopt;
    Expr prefactor = pow(geo->args()[0], d);
    // Coefficients of P (P_i = coeff of var^i); must be a clean polynomial.
    std::vector<Expr> pc;
    try {
        Poly pp(expand(mul(poly_factors)), var);
        for (const auto& cc : pp.coeffs()) {
            if (has(cc, var)) return std::nullopt;
            pc.push_back(cc);
        }
    } catch (...) {
        return std::nullopt;
    }
    const long p = static_cast<long>(pc.size()) - 1;
    if (p < 1) return std::nullopt;  // degree 0 is plain geometric
    const Expr inv_rm1 = pow(simplify(ratio - integer(1)), integer(-1));
    std::vector<Expr> q(static_cast<std::size_t>(p) + 1);
    for (long j = p; j >= 0; --j) {
        Expr s = pc[static_cast<std::size_t>(j)];
        for (long i = j + 1; i <= p; ++i) {
            s = s - ratio * binomial(integer(i), integer(j))
                        * q[static_cast<std::size_t>(i)];
        }
        q[static_cast<std::size_t>(j)] = simplify(s * inv_rm1);
    }
    auto eval_q = [&](const Expr& m) -> Expr {
        std::vector<Expr> terms;
        terms.reserve(static_cast<std::size_t>(p) + 1);
        for (long i = 0; i <= p; ++i) {
            terms.push_back(q[static_cast<std::size_t>(i)]
                            * pow(m, integer(i)));
        }
        return add(std::move(terms));
    };
    Expr s_lo = eval_q(lo) * pow(ratio, lo);
    Expr s_hi;
    if (hi->type_id() == TypeId::Infinity) {
        // Σ_{k=lo}^∞ P(k)·r^k converges iff |r| < 1 (⇔ r² < 1); then the upper
        // boundary term Q(k)·r^k → 0 (geometric decay dominates the polynomial Q),
        // so the sum is −S(lo). A divergent or undecidable ratio is left alone.
        if (!(is_negative(simplify(ratio * ratio - S::One()))
              == std::optional<bool>{true})) {
            return std::nullopt;
        }
        s_hi = S::Zero();
    } else {
        s_hi = eval_q(hi + integer(1)) * pow(ratio, hi + integer(1));
    }
    return simplify(prefactor * (s_hi - s_lo));
}

// Σ_{k=lo}^∞ c·r^k / k! = c·e^r, minus the omitted head Σ_{k=0}^{lo-1} when
// lo > 0. The exponential series converges for every r, so no convergence test
// is needed. Recognizes c · P(var) · (∏ baseᵢ^(aᵢ·var + bᵢ)) · factorial(var)^(-1):
// the bases give the rate r = ∏ baseᵢ^{aᵢ} (and a constant from the bᵢ), and the
// polynomial numerator P(var) is folded in through the falling-factorial basis,
// since Σ_{k≥0} k^{(m)} rᵏ/k! = rᵐ e^r. Writing P = Σ_m c_m·k^{(m)} then gives
// Σ P(k) rᵏ/k! = (Σ_m c_m rᵐ) e^r = Q(r) e^r. So Σ k/k! = e, Σ k²/k! = 2e, etc.
//
// Falling-factorial transform: returns Q(r) for a polynomial whose power-basis
// coefficients are `coeffs` (ascending). Extracts the monic falling factorials
// k^{(m)} = k(k−1)…(k−m+1) top-down (a triangular solve, no Stirling numbers).
[[nodiscard]] Expr exp_series_poly_transform(std::vector<Expr> coeffs,
                                             const Expr& var, const Expr& rate) {
    std::vector<Expr> q_terms;
    for (long m = static_cast<long>(coeffs.size()) - 1; m >= 0; --m) {
        Expr cm = simplify(coeffs[static_cast<std::size_t>(m)]);
        if (cm == S::Zero()) continue;
        q_terms.push_back(mul(cm, pow(rate, integer(m))));
        if (m >= 1) {
            // Subtract cm · ∏_{i=0}^{m−1}(var − i) from the remaining coeffs.
            std::vector<Expr> ff_factors;
            for (long i = 0; i < m; ++i) {
                ff_factors.push_back(add(var, integer(-i)));
            }
            Poly ff(expand(mul(std::move(ff_factors))), var);
            const auto& fc = ff.coeffs();
            for (std::size_t j = 0; j < fc.size(); ++j) {
                coeffs[j] = simplify(
                    add(coeffs[j], mul(S::NegativeOne(), mul(cm, fc[j]))));
            }
        }
    }
    return q_terms.empty() ? Expr{S::Zero()} : add(std::move(q_terms));
}

[[nodiscard]] std::optional<Expr> sum_exponential_series(
    const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    if (hi->type_id() != TypeId::Infinity) return std::nullopt;
    if (lo->type_id() != TypeId::Integer) return std::nullopt;
    const auto& loz = static_cast<const Integer&>(*lo);
    if (loz.is_negative() || !loz.fits_long()) return std::nullopt;

    std::vector<Expr> factors;
    if (expr->type_id() == TypeId::Mul) {
        for (const auto& f : expr->args()) factors.push_back(f);
    } else {
        factors.push_back(expr);
    }
    bool has_factorial = false;
    Expr cst = S::One();                // var-free constant prefactor c
    Expr rate = S::One();               // r
    std::vector<Expr> poly_factors;     // var-dependent numerator P(var)
    for (const auto& f : factors) {
        // 1/var! = factorial(var)^(-1).
        if (f->type_id() == TypeId::Pow && f->args()[1] == S::NegativeOne()
            && f->args()[0]->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f->args()[0]);
            if (fn.function_id() == FunctionId::Factorial
                && fn.args().size() == 1 && fn.args()[0] == var) {
                if (has_factorial) return std::nullopt;  // (k!)^(-2): not this
                has_factorial = true;
                continue;
            }
        }
        if (!has(f, var)) {
            cst = mul(cst, f);
            continue;
        }
        // base^(a·var + b) with base var-free: base^a → rate, base^b → const.
        if (f->type_id() == TypeId::Pow && !has(f->args()[0], var)) {
            const Expr& base = f->args()[0];
            try {
                Poly p(expand(f->args()[1]), var);
                if (p.degree() > 1) return std::nullopt;
                const auto& cf = p.coeffs();
                Expr a = cf.size() > 1 ? cf[1] : Expr{S::Zero()};
                Expr b = !cf.empty() ? cf[0] : Expr{S::Zero()};
                if (has(a, var) || has(b, var)) return std::nullopt;
                rate = mul(rate, pow(base, a));
                cst = mul(cst, pow(base, b));
                continue;
            } catch (const std::exception&) {
                return std::nullopt;
            }
        }
        // Any other var-dependent factor must be polynomial in var (the
        // numerator P): k, (k+1), k², … are folded via the falling-factorial
        // transform; a non-polynomial factor (sin(k), 1/k, …) bails.
        poly_factors.push_back(f);
    }
    if (!has_factorial) return std::nullopt;
    rate = simplify(rate);
    cst = simplify(cst);

    // Build P(var) and its closed form Q(r): Σ P(k) r^k/k! = Q(r) e^r.
    Expr poly = poly_factors.empty() ? Expr{S::One()}
                                     : mul(std::move(poly_factors));
    std::vector<Expr> pcoeffs;
    try {
        Poly p(expand(poly), var);
        pcoeffs = p.coeffs();
    } catch (const std::exception&) {
        return std::nullopt;  // numerator is not a polynomial in var
    }
    Expr q_of_r = exp_series_poly_transform(pcoeffs, var, rate);
    Expr e_r = exp(rate);
    Expr tail = mul(q_of_r, e_r);  // Σ_{k≥0} P(k) r^k/k!

    const long L = loz.to_long();
    if (L == 0) return simplify(mul(cst, tail));
    // Subtract the omitted head terms Σ_{k=0}^{L-1} P(k) r^k/k!.
    std::vector<Expr> head;
    for (long kk = 0; kk < L; ++kk) {
        head.push_back(mul({subs(poly, var, integer(kk)),
                            pow(rate, integer(kk)),
                            pow(factorial(integer(kk)), integer(-1))}));
    }
    return simplify(
        mul(cst, add(tail, mul(S::NegativeOne(), add(std::move(head))))));
}

}  // namespace

Expr summation(const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    if (!expr) return expr;

    // Single-term range (hi == lo): Σ_{k=a}^{a} f(k) = f(a).
    if (hi == lo) return simplify(subs(expr, var, lo));

    // Linearity: split Add into separate sums.
    if (expr->type_id() == TypeId::Add) {
        std::vector<Expr> partial;
        partial.reserve(expr->args().size());
        for (const auto& term : expr->args()) {
            partial.push_back(summation(term, var, lo, hi));
        }
        return add(std::move(partial));
    }

    // Pull constant factors out of a Mul.
    if (expr->type_id() == TypeId::Mul) {
        std::vector<Expr> const_factors;
        std::vector<Expr> var_factors;
        for (const auto& f : expr->args()) {
            if (has(f, var)) var_factors.push_back(f);
            else const_factors.push_back(f);
        }
        if (!const_factors.empty() && !var_factors.empty()) {
            Expr inner = summation(mul(std::move(var_factors)), var, lo, hi);
            return mul(mul(std::move(const_factors)), inner);
        }
    }

    // A product or power that expands to a sum — e.g. k·(k+1) → k²+k,
    // (k+1)² → k²+2k+1, (k+1)·2ᵏ → k·2ᵏ+2ᵏ — isn't matched by the closed forms
    // below. Expand and recurse so linearity splits it into terms each of which
    // (a monomial kᵖ, or a poly·geometric) is summable.
    if (expr->type_id() == TypeId::Mul || expr->type_id() == TypeId::Pow) {
        Expr ex = expand(expr);
        if (!(ex == expr) && ex->type_id() == TypeId::Add) {
            return summation(ex, var, lo, hi);
        }
    }

    // Constant w.r.t. var: c * (hi - lo + 1).
    if (!has(expr, var)) {
        return simplify(expr * (hi - lo + integer(1)));
    }

    // Identity Σ k from lo to hi = (hi-lo+1)(lo+hi)/2.
    if (expr == var) {
        return simplify((hi - lo + integer(1)) * (lo + hi) / integer(2));
    }

    // Σ kᵖ for a positive integer p via Faulhaber's formula, using Bernoulli
    // numbers (SymPy's B₁ = +1/2 convention):
    //   Σ_{k=1}^n kᵖ = 1/(p+1) · Σ_{j=0}^{p} C(p+1, j) B_j n^(p+1−j).
    if (expr->type_id() == TypeId::Pow && expr->args()[0] == var
        && expr->args()[1]->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*expr->args()[1]);
        if (z.fits_long()) {
            long p = z.to_long();
            if (p >= 2 && p <= 200) {  // p=1 handled above; cap for cost
                auto sum_to_n = [&](const Expr& n) -> Expr {
                    std::vector<Expr> terms;
                    terms.reserve(static_cast<std::size_t>(p) + 1);
                    for (long j = 0; j <= p; ++j) {
                        terms.push_back(mul(
                            {binomial(integer(p + 1), integer(j)),
                             bernoulli(integer(j)),
                             pow(n, integer(p + 1 - j))}));
                    }
                    return mul(pow(integer(p + 1), integer(-1)),
                               add(std::move(terms)));
                };
                Expr s_hi = sum_to_n(hi);
                Expr s_lo = sum_to_n(lo - integer(1));
                return simplify(s_hi - s_lo);
            }
        }
    }

    // Convergent p-series Σ_{k=1}^∞ 1/k^s = ζ(s) for an integer s ≥ 2. zeta()
    // closes the even cases (ζ(2)=π²/6, …) and keeps odd s as a symbolic ζ(s)
    // (matching SymPy's Sum(1/k**3).doit() = zeta(3)). The divergent harmonic
    // s=1 (m=-1) is excluded and falls through unevaluated.
    if (lo == S::One() && hi->type_id() == TypeId::Infinity
        && expr->type_id() == TypeId::Pow && expr->args()[0] == var
        && expr->args()[1]->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*expr->args()[1]);
        if (z.fits_long()) {
            const long m = z.to_long();  // summand is var^m
            if (m <= -2) return zeta(integer(-m));
        }
    }

    // Alternating p-series Σ_{k=1}^∞ (−1)^k / k^s = (2^(1−s) − 1)·ζ(s) for s ≥ 2,
    // and −log 2 for s = 1 (the Dirichlet eta values η(s)). A (−1)^(a·k+b) factor
    // with odd a is (−1)^k up to the constant sign (−1)^b; a leading constant and
    // the sign multiply through.
    if (lo == S::One() && hi->type_id() == TypeId::Infinity) {
        Expr sign_exp;
        long m = 0;
        bool have_sign = false;
        bool have_pow = false;
        Expr coeff = S::One();
        bool ok = true;
        std::vector<Expr> factors;
        if (expr->type_id() == TypeId::Mul) {
            for (const auto& f : expr->args()) factors.push_back(f);
        } else {
            factors.push_back(expr);
        }
        for (const auto& f : factors) {
            if (f->type_id() == TypeId::Pow
                && f->args()[0] == integer(-1) && has(f->args()[1], var)) {
                if (have_sign) { ok = false; break; }
                sign_exp = f->args()[1];
                have_sign = true;
            } else if (f->type_id() == TypeId::Pow && f->args()[0] == var
                       && f->args()[1]->type_id() == TypeId::Integer) {
                const auto& zp = static_cast<const Integer&>(*f->args()[1]);
                if (have_pow || !zp.fits_long()) { ok = false; break; }
                m = zp.to_long();
                have_pow = true;
            } else if (!has(f, var)) {
                coeff = mul(coeff, f);
            } else {
                ok = false;
                break;
            }
        }
        if (ok && have_sign && have_pow && m <= -1) {
            // sign_exp = a·k + b, a an odd integer (so (−1)^(a·k) = (−1)^k) and b
            // an integer (so (−1)^b is a concrete ±1).
            try {
                Poly pe(expand(sign_exp), var);
                if (pe.degree() == 1
                    && pe.coeffs()[1]->type_id() == TypeId::Integer
                    && pe.coeffs()[0]->type_id() == TypeId::Integer) {
                    const long a = static_cast<const Integer&>(*pe.coeffs()[1])
                                       .to_long();
                    const long b = static_cast<const Integer&>(*pe.coeffs()[0])
                                       .to_long();
                    if (a % 2 != 0) {
                        const long s = -m;
                        Expr sign_b = (b % 2 == 0) ? S::One() : S::NegativeOne();
                        Expr base_sum =
                            (s == 1)
                                ? Expr{mul(S::NegativeOne(), log(integer(2)))}
                                : mul(simplify(pow(integer(2), integer(1 - s))
                                               - S::One()),
                                      zeta(integer(s)));
                        return simplify(mul(mul(coeff, sign_b), base_sum));
                    }
                }
            } catch (const std::exception&) {
                // not a clean affine sign exponent — fall through
            }
        }
    }

    // Dirichlet beta Σ_{k=0}^∞ (−1)^k/(2k+1)^s: β(1) = π/4 (Leibniz), β(2) =
    // Catalan's constant. Higher s have no elementary closed form (SymPy returns a
    // polylog), so only s ∈ {1, 2} are recognized. A (−1)^(a·k+b) factor with odd
    // a is (−1)^k up to the sign (−1)^b; a leading constant multiplies through.
    if (lo == S::Zero() && hi->type_id() == TypeId::Infinity) {
        Expr sign_exp;
        long s = 0;
        bool have_sign = false;
        bool have_base = false;
        Expr coeff = S::One();
        bool ok = true;
        std::vector<Expr> factors;
        if (expr->type_id() == TypeId::Mul) {
            for (const auto& f : expr->args()) factors.push_back(f);
        } else {
            factors.push_back(expr);
        }
        for (const auto& f : factors) {
            if (f->type_id() == TypeId::Pow && f->args()[0] == integer(-1)
                && has(f->args()[1], var)) {
                if (have_sign) { ok = false; break; }
                sign_exp = f->args()[1];
                have_sign = true;
            } else if (f->type_id() == TypeId::Pow
                       && f->args()[1]->type_id() == TypeId::Integer
                       && has(f->args()[0], var)) {
                // base must be the odd-denominator affine 2·var + 1.
                const auto& ze = static_cast<const Integer&>(*f->args()[1]);
                if (have_base || !ze.fits_long() || ze.to_long() >= 0) {
                    ok = false;
                    break;
                }
                try {
                    Poly pb(expand(f->args()[0]), var);
                    if (pb.degree() == 1 && pb.coeffs()[0] == S::One()
                        && pb.coeffs()[1] == integer(2)) {
                        s = -ze.to_long();
                        have_base = true;
                    } else {
                        ok = false;
                        break;
                    }
                } catch (const std::exception&) {
                    ok = false;
                    break;
                }
            } else if (!has(f, var)) {
                coeff = mul(coeff, f);
            } else {
                ok = false;
                break;
            }
        }
        if (ok && have_sign && have_base && (s == 1 || s == 2)) {
            try {
                Poly pe(expand(sign_exp), var);
                if (pe.degree() == 1
                    && pe.coeffs()[1]->type_id() == TypeId::Integer
                    && pe.coeffs()[0]->type_id() == TypeId::Integer) {
                    const long a = static_cast<const Integer&>(*pe.coeffs()[1])
                                       .to_long();
                    const long b = static_cast<const Integer&>(*pe.coeffs()[0])
                                       .to_long();
                    if (a % 2 != 0) {
                        Expr sign_b = (b % 2 == 0) ? S::One() : S::NegativeOne();
                        Expr val = (s == 1)
                                       ? Expr{S::Pi() / integer(4)}
                                       : S::Catalan();
                        return simplify(mul(mul(coeff, sign_b), val));
                    }
                }
            } catch (const std::exception&) {
                // not a clean affine sign exponent — fall through
            }
        }
    }

    // Exponential series Σ_{k=lo}^∞ r^k/k! = e^r (e.g. Σ 1/k! = e,
    // Σ x^k/k! = e^x, Σ 2^k/k! = e²). Convergent for every r.
    if (auto r = sum_exponential_series(expr, var, lo, hi)) return *r;

    // Geometric: base^(c*var + d) with base, c, d all independent of var.
    // Rewrites as base^d * (base^c)^var = A * ratio^var, a geometric series
    // with ratio = base^c and constant prefactor A = base^d. This subsumes
    // the plain base^var case (c=1, d=0) and also handles base^(-var),
    // base^(2*var), etc. — exponents the previous exact-match check missed.
    if (expr->type_id() == TypeId::Pow) {
        const Expr& base = expr->args()[0];
        const Expr& exponent = expr->args()[1];
        if (!has(base, var) && has(exponent, var)) {
            Expr c = diff(exponent, var);
            if (!has(c, var)) {  // exponent is linear in var
                Expr d = simplify(exponent - c * var);
                if (!has(d, var)) {
                    Expr ratio = pow(base, c);
                    Expr prefactor = pow(base, d);
                    // A * (ratio^lo - ratio^(hi+1)) / (1 - ratio)
                    Expr num = pow(ratio, lo) - pow(ratio, hi + integer(1));
                    Expr den = integer(1) - ratio;
                    return simplify(prefactor * num / den);
                }
            }
        }
    }

    // Arithmetic-geometric: Σ k·r^k where r = base^c is a concrete numeric
    // ratio ≠ 1 (k times a geometric term). Restricted to a numeric base so the
    // result is the plain closed form, matching SymPy (a symbolic ratio would
    // need a Piecewise on r = 1). An infinite upper bound is left to the general
    // poly·geometric path below, which has the |r| < 1 convergence handling.
    if (expr->type_id() == TypeId::Mul && expr->args().size() == 2
        && hi->type_id() != TypeId::Infinity) {
        const Expr& f0 = expr->args()[0];
        const Expr& f1 = expr->args()[1];
        // Identify the bare `var` factor and the geometric `base^(c·var+d)`.
        const Expr* geo = nullptr;
        if (f0 == var && f1->type_id() == TypeId::Pow) {
            geo = &f1;
        } else if (f1 == var && f0->type_id() == TypeId::Pow) {
            geo = &f0;
        }
        if (geo != nullptr) {
            const Expr& base = (*geo)->args()[0];
            const Expr& exponent = (*geo)->args()[1];
            if (is_number(base) && has(exponent, var)) {
                Expr c = diff(exponent, var);
                if (!has(c, var)) {  // exponent linear in var
                    Expr d = simplify(exponent - c * var);
                    Expr ratio = pow(base, c);
                    if (!has(d, var) && is_number(ratio)
                        && !(ratio == S::One())) {
                        // Σ_{k=0}^{N} k·r^k = r(1 − (N+1)r^N + N·r^{N+1})/(1−r)².
                        // The base^d prefactor of the geometric term factors out.
                        Expr prefactor = pow(base, d);
                        Expr one_minus_r_sq =
                            pow(integer(1) - ratio, integer(2));
                        auto partial = [&](const Expr& n) -> Expr {
                            return ratio
                                   * (integer(1) - (n + integer(1)) * pow(ratio, n)
                                      + n * pow(ratio, n + integer(1)))
                                   / one_minus_r_sq;
                        };
                        Expr s = partial(hi) - partial(lo - integer(1));
                        return simplify(prefactor * s);
                    }
                }
            }
        }
    }

    // Polynomial × geometric: Σ P(k)·r^k for a polynomial P and concrete ratio
    // r ≠ 1 (covers k²·2^k, (2k+1)·3^k, k³/2^k — the general arithmetic-geometric
    // case the degree-1 block above doesn't reach).
    if (auto pg = sum_poly_geometric(expr, var, lo, hi); pg) return *pg;

    // Telescoping rational summand (e.g. 1/(k(k+1)), 1/(4k²−1)).
    if (auto t = telescope_rational(expr, var, lo, hi); t) return *t;

    // No closed form found — return the unevaluated Sum marker rather than the
    // bare summand (Σ 1/k² must not collapse to 1/k²).
    return sum_marker(expr, var, lo, hi);
}

}  // namespace sympp
