#include <sympp/calculus/summation.hpp>

#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

Expr summation(const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    if (!expr) return expr;

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

    // Geometric: r^var with r independent of var.
    if (expr->type_id() == TypeId::Pow && expr->args()[1] == var) {
        const Expr& r = expr->args()[0];
        if (!has(r, var)) {
            // (r^lo - r^(hi+1)) / (1 - r)
            Expr num = pow(r, lo) - pow(r, hi + integer(1));
            Expr den = integer(1) - r;
            return simplify(num / den);
        }
    }

    // Could not recognize. Return unchanged.
    return expr;
}

}  // namespace sympp
