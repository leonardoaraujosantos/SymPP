#include <sympp/calculus/limit.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/infinity.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

namespace {

// Split a Mul expression of the form num * Pow(den, -1) * ... into
// (numerator, denominator) where denominator is the product of all
// negative-exponent-1 factors. Used after together() so the input is
// already a single fraction.
struct NumDen { Expr num; Expr den; };

[[nodiscard]] NumDen split_after_together(const Expr& e) {
    Expr t = together(e);
    if (t->type_id() != TypeId::Mul) return {t, S::One()};
    std::vector<Expr> num_factors;
    std::vector<Expr> den_factors;
    for (const auto& f : t->args()) {
        // Match any factor with a negative-integer exponent; collect its
        // base raised to the absolute exponent into the denominator. This
        // covers both Pow(d, -1) and Pow(d, -k) cases (e.g. Pow(x, -2)
        // after the canonical (a^m)^n → a^(m·n) fold).
        if (f->type_id() == TypeId::Pow
            && f->args()[1]->type_id() == TypeId::Integer) {
            const auto& z = static_cast<const Integer&>(*f->args()[1]);
            if (z.is_negative()) {
                Expr abs_exp = mul(S::NegativeOne(), f->args()[1]);
                den_factors.push_back(pow(f->args()[0], abs_exp));
                continue;
            }
        }
        num_factors.push_back(f);
    }
    if (den_factors.empty()) return {t, S::One()};
    return {mul(num_factors), mul(den_factors)};
}

[[nodiscard]] bool is_nan(const Expr& e) noexcept {
    return e->type_id() == TypeId::NaN;
}

// L'Hôpital on the fraction num/den. Resolves the 0/0 and ∞/∞ indeterminate
// forms by differentiating top and bottom, re-rationalising the ratio each
// step (so algebraic cancellation, e.g. x²/(x²+x) → x/(x+1), keeps the limit
// tractable). Returns the value, or std::nullopt when it cannot progress.
[[nodiscard]] std::optional<Expr> lhopital_nd(Expr num, Expr den,
                                              const Expr& var,
                                              const Expr& target) {
    // "Infinite-like": a genuine infinity, or nan arising from oo±oo when a
    // polynomial is substituted at the target (its true magnitude is ∞).
    auto infinite_like = [](const Expr& v) noexcept {
        return is_infinity(v) || v->type_id() == TypeId::NaN;
    };
    for (int iter = 0; iter < 16; ++iter) {
        Expr num_at = simplify(subs(num, var, target));
        Expr den_at = simplify(subs(den, var, target));
        const bool zero_zero = (num_at == S::Zero() && den_at == S::Zero());
        const bool inf_inf = infinite_like(num_at) && infinite_like(den_at);
        if (zero_zero || inf_inf) {
            Expr num_d = diff(num, var);
            Expr den_d = diff(den, var);
            if (den_d == S::Zero()) return std::nullopt;
            // Re-rationalise the new ratio with together() (cheap) before the
            // next step; full simplify() here is too slow over many iterations.
            auto nd = split_after_together(num_d / den_d);
            num = std::move(nd.num);
            den = std::move(nd.den);
            continue;
        }
        // Determinate denominator (finite non-zero or an infinity): the value
        // is num/den, which the infinity arithmetic folds (finite/∞ = 0, …).
        if (!is_nan(den_at) && !(den_at == S::Zero())) {
            Expr val = simplify(num_at / den_at);
            if (!is_nan(val)) return val;
        }
        return std::nullopt;
    }
    return std::nullopt;
}

// L'Hôpital on `expr` first split into a single fraction.
[[nodiscard]] std::optional<Expr> lhopital(const Expr& expr, const Expr& var,
                                           const Expr& target) {
    auto nd = split_after_together(expr);
    return lhopital_nd(nd.num, nd.den, var, target);
}

[[nodiscard]] Expr limit_impl(const Expr& expr, const Expr& var,
                              const Expr& target, int depth);

// exp() of a (possibly infinite) limit value: exp(+oo)=oo, exp(-oo)=0.
[[nodiscard]] Expr exp_of_limit(const Expr& v) {
    if (v->type_id() == TypeId::Infinity) return S::Infinity();
    if (v->type_id() == TypeId::NegativeInfinity) return S::Zero();
    return simplify(exp(v));
}

// Power indeterminate forms 1^∞, ∞^0, 0^0 via a^b = exp(b·log a).
[[nodiscard]] std::optional<Expr> try_power_form(const Expr& expr,
                                                 const Expr& var,
                                                 const Expr& target,
                                                 int depth) {
    if (expr->type_id() != TypeId::Pow) return std::nullopt;
    const Expr& base = expr->args()[0];
    const Expr& ex = expr->args()[1];
    Expr bl = limit_impl(base, var, target, depth + 1);
    Expr el = limit_impl(ex, var, target, depth + 1);
    const bool one_inf = (bl == S::One() && is_infinity(el));
    const bool inf_zero = (is_infinity(bl) && el == S::Zero());
    const bool zero_zero = (bl == S::Zero() && el == S::Zero());
    if (!(one_inf || inf_zero || zero_zero)) return std::nullopt;
    Expr inner = limit_impl(mul(ex, log(base)), var, target, depth + 1);
    if (is_nan(inner)) return std::nullopt;
    return exp_of_limit(inner);
}

// Product indeterminate form 0·∞: rewrite the factors that tend to ∞ into the
// denominator (as reciprocals) and resolve the resulting quotient.
[[nodiscard]] std::optional<Expr> try_product_form(const Expr& expr,
                                                   const Expr& var,
                                                   const Expr& target,
                                                   int depth) {
    if (expr->type_id() != TypeId::Mul) return std::nullopt;
    std::vector<Expr> num_factors;
    std::vector<Expr> den_factors;
    bool saw_zero = false;
    bool saw_inf = false;
    for (const auto& f : expr->args()) {
        Expr lf = limit_impl(f, var, target, depth + 1);
        if (is_nan(lf)) return std::nullopt;
        if (is_infinity(lf)) {
            den_factors.push_back(pow(f, S::NegativeOne()));
            saw_inf = true;
        } else {
            num_factors.push_back(f);
            if (lf == S::Zero()) saw_zero = true;
        }
    }
    if (!(saw_zero && saw_inf)) return std::nullopt;
    // f·g with f→0, g→∞ becomes f / (1/g): both numerator and (1/g) → 0, a
    // 0/0 form. den_factors already hold the 1/g reciprocals, so their product
    // is the denominator (→ 0); the remaining factors form the numerator.
    Expr num = num_factors.empty() ? S::One() : mul(num_factors);
    Expr den = den_factors.empty() ? S::One() : mul(den_factors);
    return lhopital_nd(num, den, var, target);
}

// Direct substitution at a finite pole yields zoo (the unsigned 1/0). Recover
// the signed infinity by sampling the sign of the expression just either side of
// the target: equal signs ⇒ ±oo (1/x² → +oo, −1/x² → −oo), opposite signs ⇒ the
// limit is genuinely two-sided and stays zoo (1/x). Inconclusive samples (a
// non-real value, or no definite sign) also keep zoo.
[[nodiscard]] std::optional<Expr> signed_pole(const Expr& expr, const Expr& var,
                                              const Expr& target) {
    if (!is_number(target)) return std::nullopt;
    auto side_sign = [&](const Expr& point) -> int {
        Expr v = simplify(subs(expr, var, point));
        if (v->type_id() == TypeId::Infinity) return 1;
        if (v->type_id() == TypeId::NegativeInfinity) return -1;
        if (is_infinity(v) || v->type_id() == TypeId::NaN) return 0;
        Expr f = evalf(v, 30);
        if (is_positive(f) == std::optional<bool>{true}) return 1;
        if (is_negative(f) == std::optional<bool>{true}) return -1;
        return 0;
    };
    Expr h = pow(integer(10), integer(-6));  // 1e-6
    int s_right = side_sign(add(target, h));
    int s_left = side_sign(add(target, mul(S::NegativeOne(), h)));
    if (s_right == 0 || s_left == 0 || s_right != s_left) return std::nullopt;
    return s_right > 0 ? S::Infinity() : S::NegativeInfinity();
}

Expr limit_impl(const Expr& expr, const Expr& var, const Expr& target,
                int depth) {
    Expr direct = simplify(subs(expr, var, target));

    // A finite-target pole surfaces as zoo; resolve its sign when both sides
    // agree, otherwise leave the unsigned zoo.
    if (direct->type_id() == TypeId::ComplexInfinity) {
        if (auto s = signed_pole(expr, var, target)) return *s;
    }

    if (depth < 12) {
        // Indeterminate forms surface as nan after substitution.
        if (is_nan(direct)) {
            if (auto v = try_power_form(expr, var, target, depth)) return *v;
            if (auto v = try_product_form(expr, var, target, depth)) return *v;
        }
        // 0/0 and ∞/∞ quotients (also recovers finite 0/0 where direct
        // substitution collapses to 0 or nan).
        if (auto v = lhopital(expr, var, target)) return *v;
    }
    return direct;
}

}  // namespace

Expr limit(const Expr& expr, const Expr& var, const Expr& target) {
    return limit_impl(expr, var, target, 0);
}

}  // namespace sympp
