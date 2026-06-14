#include <sympp/calculus/summation.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/expand.hpp>
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

    // Telescoping rational summand (e.g. 1/(k(k+1)), 1/(4k²−1)).
    if (auto t = telescope_rational(expr, var, lo, hi); t) return *t;

    // No closed form found — return the unevaluated Sum marker rather than the
    // bare summand (Σ 1/k² must not collapse to 1/k²).
    return sum_marker(expr, var, lo, hi);
}

}  // namespace sympp
