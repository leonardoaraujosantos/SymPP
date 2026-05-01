#pragma once

// euler_lagrange — variational-calculus functional derivative.
//
// For a Lagrangian L(y, y', y'', ...; x), the Euler-Lagrange equation
// for a single first-derivative dependency is
//
//     ∂L/∂y  -  D_x ( ∂L/∂y' )  =  0
//
// where D_x is the total derivative w.r.t. x (chain rule):
//
//     D_x f = ∂f/∂x + (∂f/∂y) y' + (∂f/∂y') y''
//
// API: pass L as an expression in symbols representing y, y', y'',
// and x. The function returns the LHS of the EL equation as an
// expression — set it equal to zero to get the equation of motion.
//
// Reference: sympy/calculus/euler.py::euler_equations

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr euler_lagrange(
    const Expr& L,
    const Expr& y,
    const Expr& yprime,
    const Expr& ydoubleprime,
    const Expr& x);

}  // namespace sympp
