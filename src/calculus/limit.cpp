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
    Expr den;
    for (const auto& f : t->args()) {
        if (f->type_id() == TypeId::Pow
            && f->args()[1] == S::NegativeOne()
            && !den) {
            den = f->args()[0];
            continue;
        }
        num_factors.push_back(f);
    }
    if (!den) return {t, S::One()};
    return {mul(num_factors), den};
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
