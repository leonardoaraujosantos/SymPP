#pragma once

// MATLAB facade — PDE solving.
//
// MATLAB's PDE Toolbox is a separate product with its own API
// (`pdepe`, `solvepde`, etc.). The naming here mirrors SymPy's
// `pdsolve`, which is what most symbolic users will be reaching for
// when they ask "solve this PDE symbolically." The thin wrappers
// re-expose SymPP's PDE entry points under MATLAB-friendly names so
// callers can `using namespace sympp::matlab` and stay consistent
// with `dsolve` / `solve`.
//
// Coverage shipped:
//   pdsolve(a, b, c, x, y)       — constant-coefficient first-order
//                                   linear PDE a·u_x + b·u_y = c via
//                                   method of characteristics.
//   pdsolve_variable(a, b, c, y, yp, x)
//                                — variable-coefficient first-order
//                                   linear PDE; homogeneous case
//                                   resolved, inhomogeneous returns
//                                   an unevaluated Pdsolve marker.
//   pdsolve_heat(k, lambda, x, t) — heat equation u_t = k·u_xx via
//                                    separation of variables.
//   pdsolve_wave(c, x, t)        — wave equation u_tt = c²·u_xx via
//                                    d'Alembert.
//
// Reference: sympp/ode/dsolve.hpp (PDE block).

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>
#include <sympp/ode/dsolve.hpp>

namespace sympp::matlab {

// pdsolve(a, b, c, x, y) — solve a·u_x + b·u_y = c (constants).
[[nodiscard]] inline Expr pdsolve(const Expr& a, const Expr& b, const Expr& c,
                                     const Expr& x, const Expr& y) {
    return pdsolve_first_order_linear(a, b, c, x, y);
}

// pdsolve_variable(a, b, c, y, yp, x) — variable-coefficient
// first-order linear PDE a(x,y)·u_x + b(x,y)·u_y = c(x,y). `y` and
// `yp` are the symbolic stand-ins for the dependent variable and its
// derivative used internally to integrate the characteristic ODE.
[[nodiscard]] inline Expr pdsolve_variable(
    const Expr& a, const Expr& b, const Expr& c,
    const Expr& y, const Expr& yp, const Expr& x) {
    return pdsolve_first_order_variable(a, b, c, y, yp, x);
}

// pdsolve_heat(k, lambda, x, t) — heat equation u_t = k·u_xx;
// `lambda` is the separation constant.
[[nodiscard]] inline Expr pdsolve_heat(const Expr& k, const Expr& lambda,
                                          const Expr& x, const Expr& t) {
    return sympp::pdsolve_heat(k, lambda, x, t);
}

// pdsolve_wave(c, x, t) — wave equation u_tt = c²·u_xx via d'Alembert.
[[nodiscard]] inline Expr pdsolve_wave(const Expr& c, const Expr& x,
                                          const Expr& t) {
    return sympp::pdsolve_wave(c, x, t);
}

}  // namespace sympp::matlab
