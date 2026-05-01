#include <sympp/calculus/euler_lagrange.hpp>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

Expr euler_lagrange(const Expr& L, const Expr& y, const Expr& yprime,
                    const Expr& ydoubleprime, const Expr& x) {
    // ∂L/∂y
    Expr Ly = diff(L, y);
    // ∂L/∂y'
    Expr Lyp = diff(L, yprime);
    // D_x(∂L/∂y') = ∂Lyp/∂x + (∂Lyp/∂y)*y' + (∂Lyp/∂y')*y''
    Expr Dx = diff(Lyp, x)
              + diff(Lyp, y) * yprime
              + diff(Lyp, yprime) * ydoubleprime;
    return simplify(Ly - Dx);
}

}  // namespace sympp
