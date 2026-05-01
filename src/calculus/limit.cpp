#include <sympp/calculus/limit.hpp>

#include <sympp/core/traversal.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

Expr limit(const Expr& expr, const Expr& var, const Expr& target) {
    // Direct substitution + canonical simplify. Indeterminate forms aren't
    // handled here — caller can call diff (l'Hôpital) manually.
    return simplify(subs(expr, var, target));
}

}  // namespace sympp
