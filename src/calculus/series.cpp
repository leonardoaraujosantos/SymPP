#include <sympp/calculus/series.hpp>

#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/calculus/limit.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/combinatorial.hpp>

namespace sympp {

namespace {
// A non-finite Taylor coefficient — the function is singular at the expansion
// point (e.g. log(x), 1/x at 0), or direct substitution hit an indeterminate
// 0/0 form (a removable singularity such as sin(x)/x).
[[nodiscard]] bool is_non_finite(const Expr& e) {
    switch (e->type_id()) {
        case TypeId::ComplexInfinity:
        case TypeId::Infinity:
        case TypeId::NegativeInfinity:
        case TypeId::NaN:
            return true;
        default:
            return false;
    }
}
}  // namespace

Expr series(const Expr& expr, const Expr& var, const Expr& x0, std::size_t n) {
    if (n == 0) return S::Zero();
    std::vector<Expr> terms;
    terms.reserve(n);

    Expr current_deriv = expr;
    Expr dx = var - x0;
    for (std::size_t k = 0; k < n; ++k) {
        // Evaluate the k-th derivative at x = x0.
        Expr value_at_x0 = subs(current_deriv, var, x0);
        if (is_non_finite(value_at_x0)) {
            // A removable singularity (0/0 like sin(x)/x) has a finite limit and
            // a valid Taylor coefficient; a genuine singularity (log, 1/x) does
            // not, so the series cannot be formed — return the input unexpanded,
            // matching SymPy (series(log(x), x, 0) = log(x)).
            value_at_x0 = limit(current_deriv, var, x0);
            if (is_non_finite(value_at_x0)) return expr;
        }
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
