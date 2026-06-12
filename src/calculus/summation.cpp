#include <sympp/calculus/summation.hpp>

#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/core/undefined_function.hpp>
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

    // Constant w.r.t. var: c * (hi - lo + 1).
    if (!has(expr, var)) {
        return simplify(expr * (hi - lo + integer(1)));
    }

    // Identity Σ k from lo to hi = (hi-lo+1)(lo+hi)/2.
    if (expr == var) {
        return simplify((hi - lo + integer(1)) * (lo + hi) / integer(2));
    }

    // Σ k^p with integer p ∈ {2, 3} via the standard formulas (telescoped
    // from 1..n form).
    if (expr->type_id() == TypeId::Pow && expr->args()[0] == var
        && expr->args()[1]->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*expr->args()[1]);
        if (z.fits_long()) {
            long p = z.to_long();
            auto sum_to_n = [&](const Expr& n) -> Expr {
                if (p == 2) {
                    // n(n+1)(2n+1)/6
                    return n * (n + integer(1))
                           * (integer(2) * n + integer(1)) / integer(6);
                }
                if (p == 3) {
                    // (n(n+1)/2)²
                    return pow(n * (n + integer(1)) / integer(2), integer(2));
                }
                return Expr{};
            };
            if (p == 2 || p == 3) {
                Expr s_hi = sum_to_n(hi);
                Expr s_lo = sum_to_n(lo - integer(1));
                return simplify(s_hi - s_lo);
            }
        }
    }

    // Convergent even p-series Σ_{k=1}^∞ 1/k^(2n) = ζ(2n) = r·π^(2n) (Basel and
    // its relatives). Only even exponents have an elementary closed form; odd
    // p>1 (ζ(3), …) and the divergent p≤1 cases fall through unevaluated.
    if (lo == S::One() && hi->type_id() == TypeId::Infinity
        && expr->type_id() == TypeId::Pow && expr->args()[0] == var
        && expr->args()[1]->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*expr->args()[1]);
        if (z.fits_long()) {
            const long m = z.to_long();  // summand is var^m
            if (m <= -2 && m % 2 == 0) {
                const long twon = -m;
                Expr coeff;  // rational coefficient of π^(2n); null if untabled
                switch (twon) {
                    case 2:  coeff = rational(1, 6); break;
                    case 4:  coeff = rational(1, 90); break;
                    case 6:  coeff = rational(1, 945); break;
                    case 8:  coeff = rational(1, 9450); break;
                    case 10: coeff = rational(1, 93555); break;
                    case 12: coeff = rational(691, 638512875); break;
                    case 14: coeff = rational(2, 18243225); break;
                    default: break;
                }
                if (coeff) return mul(coeff, pow(S::Pi(), integer(twon)));
            }
        }
    }

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
    // need a Piecewise on r = 1). Higher-degree P(k)·r^k stays unevaluated.
    if (expr->type_id() == TypeId::Mul && expr->args().size() == 2) {
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

    // No closed form found — return the unevaluated Sum marker rather than the
    // bare summand (Σ 1/k² must not collapse to 1/k²).
    return sum_marker(expr, var, lo, hi);
}

}  // namespace sympp
