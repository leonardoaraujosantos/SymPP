#include <sympp/calculus/series.hpp>

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
#include <sympp/functions/combinatorial.hpp>

namespace sympp {

Expr series(const Expr& expr, const Expr& var, const Expr& x0, std::size_t n) {
    if (n == 0) return S::Zero();
    std::vector<Expr> terms;
    terms.reserve(n);

    Expr current_deriv = expr;
    Expr dx = var - x0;
    for (std::size_t k = 0; k < n; ++k) {
        // Evaluate the k-th derivative at x = x0.
        Expr value_at_x0 = subs(current_deriv, var, x0);
        // Build (1/k!) * value * (x - x0)^k
        Expr fact = factorial(integer(static_cast<long>(k)));
        Expr term = mul({pow(fact, S::NegativeOne()),
                         value_at_x0,
                         pow(dx, integer(static_cast<long>(k)))});
        terms.push_back(std::move(term));
        // Prepare next derivative.
        if (k + 1 < n) {
            current_deriv = diff(current_deriv, var);
        }
    }
    return add(std::move(terms));
}

}  // namespace sympp
