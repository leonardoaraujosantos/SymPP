#pragma once

// MATLAB facade — ODE / DAE solving.
//
// MATLAB's `dsolve(eq, y)` introspects an unevaluated `diff(y(x), x)`
// form to detect ODE order. SymPP's AST has no unevaluated
// `Derivative` node — `diff` always returns the derivative
// immediately. So the user must pass the derivative symbols
// explicitly, matching the existing core API:
//
//     auto [x, y, yp]      = syms("x", "y", "yp");
//     auto [x, y, yp, ypp] = syms("x", "y", "yp", "ypp");  // 2nd order
//
// The four overloads below cover the common cases:
//
//   dsolve(eq, y, yp, x)       — first-order. Routes through the
//                                 first-order classifier (separable,
//                                 linear, Bernoulli, exact, Riccati,
//                                 homogeneous, autonomous).
//   dsolve(eq, y, yp, ypp, x)  — second-order. Auto-classifies as
//                                 constant-coefficient or
//                                 Cauchy-Euler when linear; falls
//                                 through to the unevaluated marker
//                                 otherwise.
//   dsolve_constant_coeff(coeffs, x)
//                              — re-export.  coeffs[k] is the
//                                 coefficient of y^(k); equation is
//                                 Σ coeffs[k] · y^(k) = 0.
//   dsolve(A, x)               — linear ODE system y' = A·y.
//
// Each has an `_ivp` companion that applies a list of initial
// conditions to the general solution via apply_ivp.
//
// Reference: sympp/ode/dsolve.hpp.

#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/api.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/core/undefined_function.hpp>
#include <sympp/fwd.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/ode/dsolve.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp::matlab {

namespace detail {

// Extract a coefficient by partial differentiation. For a linear ODE
// the derivative w.r.t. any single derivative-symbol gives that
// coefficient; for nonlinear inputs the result depends on the other
// symbols, which the caller must check.
[[nodiscard]] inline Expr partial(const Expr& eq, const Expr& d) {
    return simplify(diff(eq, d));
}

// True iff `coef` matches `factor · x^expected_pow` for some
// constant `factor`. Used by the Cauchy-Euler detector: in
// a x^2 y'' + b x y' + c y = 0, the coefficient of y^(k) is b_k · x^k.
[[nodiscard]] inline bool is_x_to_the(const Expr& coef, const Expr& x,
                                         long expected_pow) {
    if (expected_pow == 0) return !has(coef, x);
    Expr factor = simplify(coef / pow(x, integer(expected_pow)));
    return !has(factor, x);
}

}  // namespace detail

// dsolve(eq, y, yp, x) — first-order.
[[nodiscard]] inline Expr dsolve(const Expr& eq, const Expr& y,
                                    const Expr& yp, const Expr& x) {
    return dsolve_first_order(eq, y, yp, x);
}

// dsolve_ivp(eq, y, yp, x, conditions) — first-order with initial
// conditions applied.
[[nodiscard]] inline Expr dsolve_ivp(
    const Expr& eq, const Expr& y, const Expr& yp, const Expr& x,
    const std::vector<std::pair<Expr, Expr>>& conditions) {
    return apply_ivp(dsolve_first_order(eq, y, yp, x), x, conditions);
}

// dsolve(eq, y, yp, ypp, x) — second-order, auto-classified.
//
// Strategy:
//   1. ∂eq/∂ypp, ∂eq/∂yp, ∂eq/∂y give the coefficients of ypp, yp, y.
//      If any of those still depend on y/yp/ypp the equation is not
//      linear → return the unevaluated Dsolve(...) marker.
//   2. Extract the RHS  g(x) = −eq|y=yp=ypp=0  so that the equation
//      reads a₂·ypp + a₁·yp + a₀·y = g(x). When g is non-zero the
//      nonhomogeneous-variant solvers (variation of parameters) are
//      called; otherwise the homogeneous ones.
//   3. Coefficients all constant in x → dsolve_constant_coeff[_nonhomogeneous].
//   4. Coefficients match a_k · x^k → dsolve_cauchy_euler[_nonhomogeneous].
//   5. Otherwise → unevaluated marker (variable-coefficient ODEs
//      with non-Cauchy-Euler shape are out of scope here).
[[nodiscard]] inline Expr dsolve(const Expr& eq, const Expr& y,
                                    const Expr& yp, const Expr& ypp,
                                    const Expr& x) {
    Expr a2 = detail::partial(eq, ypp);
    Expr a1 = detail::partial(eq, yp);
    Expr a0 = detail::partial(eq, y);

    // Linearity check.
    if (has(a2, y) || has(a2, yp) || has(a2, ypp)
        || has(a1, y) || has(a1, yp) || has(a1, ypp)
        || has(a0, y) || has(a0, yp) || has(a0, ypp)) {
        return function_symbol("Dsolve")(std::vector<Expr>{eq, y, x});
    }

    // Extract RHS g(x): zero out y, yp, ypp in eq → −g(x), so g = −that.
    Expr eq_at_zero = subs(subs(subs(eq, y, S::Zero()),
                                   yp, S::Zero()),
                            ypp, S::Zero());
    Expr g = simplify(-eq_at_zero);

    bool const_coef = !has(a0, x) && !has(a1, x) && !has(a2, x);
    if (const_coef) {
        if (g == S::Zero()) {
            return dsolve_constant_coeff({a0, a1, a2}, x);
        }
        return dsolve_constant_coeff_nonhomogeneous({a0, a1, a2}, g, x);
    }
    bool cauchy_euler = detail::is_x_to_the(a0, x, 0)
                         && detail::is_x_to_the(a1, x, 1)
                         && detail::is_x_to_the(a2, x, 2);
    if (cauchy_euler) {
        Expr c0 = a0;
        Expr c1 = simplify(a1 / x);
        Expr c2 = simplify(a2 / pow(x, integer(2)));
        if (g == S::Zero()) {
            return dsolve_cauchy_euler({c0, c1, c2}, x);
        }
        return dsolve_cauchy_euler_nonhomogeneous({c0, c1, c2}, g, x);
    }
    return function_symbol("Dsolve")(std::vector<Expr>{eq, y, x});
}

// dsolve_ivp — second-order with initial conditions.
[[nodiscard]] inline Expr dsolve_ivp(
    const Expr& eq, const Expr& y, const Expr& yp, const Expr& ypp,
    const Expr& x,
    const std::vector<std::pair<Expr, Expr>>& conditions) {
    return apply_ivp(dsolve(eq, y, yp, ypp, x), x, conditions);
}

// dsolve_constant_coeff(coeffs, x) — re-export.
[[nodiscard]] inline Expr dsolve_constant_coeff(
    const std::vector<Expr>& coeffs, const Expr& x) {
    return sympp::dsolve_constant_coeff(coeffs, x);
}

// dsolve(A, x) — linear ODE system y' = A·y. Returns a Matrix
// general solution (column vector).
[[nodiscard]] inline Matrix dsolve(const Matrix& A, const Expr& x) {
    return dsolve_system(A, x);
}

}  // namespace sympp::matlab
