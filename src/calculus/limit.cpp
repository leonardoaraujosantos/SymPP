#include <sympp/calculus/limit.hpp>

#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
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

}  // namespace

Expr limit(const Expr& expr, const Expr& var, const Expr& target) {
    // Direct-substitute first; if the value is finite and determinate, we're
    // done.
    Expr direct = simplify(subs(expr, var, target));

    // Detect 0/0 indeterminate forms by inspecting numerator / denominator
    // separately at the target. Apply L'Hôpital iteratively up to a fixed
    // bound (real-world textbook limits rarely need more than a couple).
    Expr current = expr;
    for (int iter = 0; iter < 10; ++iter) {
        auto nd = split_after_together(current);
        Expr num_at = simplify(subs(nd.num, var, target));
        Expr den_at = simplify(subs(nd.den, var, target));
        if (num_at == S::Zero() && den_at == S::Zero()) {
            // 0/0 — L'Hôpital: limit(f/g) = limit(f'/g') when both → 0.
            Expr num_d = diff(nd.num, var);
            Expr den_d = diff(nd.den, var);
            if (den_d == S::Zero()) break;  // can't recurse meaningfully
            current = num_d / den_d;
            continue;
        }
        // Non-zero denominator, finite numerator: we have the value.
        if (!(den_at == S::Zero())) {
            return simplify(num_at / den_at);
        }
        // Numerator nonzero, denominator zero: ±∞ (we don't model infinity,
        // so fall through and return the direct substitution which will
        // contain a division by zero).
        break;
    }
    return direct;
}

}  // namespace sympp
