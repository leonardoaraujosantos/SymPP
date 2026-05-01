#pragma once

// dsolve — symbolic ODE solver.
//
// API convention: the unknown function y(x), its derivatives y', y'', ...
// are passed as ordinary symbols. This matches euler_lagrange(L, y, yp,
// ypp, x) and avoids introducing a separate AppliedFunction layer for now.
//
// First-order ODEs:
//     dsolve_first_order(eq, y, yp, x)
//   `eq` is an Expr in {y, yp, x} representing eq = 0. Returns the general
//   solution as an Expr in x with free constants __C0, __C1, ... or an
//   unevaluated Dsolve(...) marker if no strategy matches.
//
// Higher-order constant-coefficient linear ODEs:
//     dsolve_constant_coeff(coeffs, x)
//   Pass coefficients [a_0, a_1, ..., a_n] for
//   a_n y^(n) + ... + a_1 y' + a_0 y = 0
//
// Cauchy-Euler:
//     dsolve_cauchy_euler(coeffs, x)
//   Pass coefficients [a_0, a_1, ..., a_n] for
//   a_n x^n y^(n) + ... + a_1 x y' + a_0 y = 0
//
// Linear ODE systems:
//     dsolve_system(A, x)   for y' = A·y
//
// PDE:
//     pdsolve_first_order_linear(a, b, c, ux, uy, x, y)
//   Solves a(x,y)·u_x + b(x,y)·u_y = c(x,y).
//
// Apply IVP:
//     apply_ivp(general_solution, x, conditions)
//
// Verify:
//     checkodesol(eq, sol, y, yp, x)  — returns true if substituting
//   sol (and its derivative) into eq yields zero.
//
// Reference: sympy/solvers/ode/{ode,single,systems,nonhomogeneous}.py

#include <utility>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>
#include <sympp/matrices/matrix.hpp>

namespace sympp {

// --- First-order classifier + dispatch ---
[[nodiscard]] SYMPP_EXPORT Expr dsolve_first_order(
    const Expr& eq, const Expr& y, const Expr& yp, const Expr& x);

// --- Per-class first-order solvers (return nullopt-equivalent unevaluated
//     marker on failure) ---
[[nodiscard]] SYMPP_EXPORT Expr dsolve_separable(
    const Expr& eq, const Expr& y, const Expr& yp, const Expr& x);
[[nodiscard]] SYMPP_EXPORT Expr dsolve_linear_first_order(
    const Expr& eq, const Expr& y, const Expr& yp, const Expr& x);
[[nodiscard]] SYMPP_EXPORT Expr dsolve_bernoulli(
    const Expr& eq, const Expr& y, const Expr& yp, const Expr& x);
[[nodiscard]] SYMPP_EXPORT Expr dsolve_exact(
    const Expr& M, const Expr& N, const Expr& y, const Expr& x);

// --- Higher-order ---
// Returns y(x) general solution Σ Cᵢ exp(rᵢ x) (with x^k multipliers for
// repeated roots) for the homogeneous constant-coefficient ODE
//   coeffs[n] · y^(n) + ... + coeffs[1] · y' + coeffs[0] · y = 0.
[[nodiscard]] SYMPP_EXPORT Expr dsolve_constant_coeff(
    const std::vector<Expr>& coeffs, const Expr& x);

// Cauchy-Euler: coeffs[n] x^n y^(n) + ... + coeffs[0] y = 0.
[[nodiscard]] SYMPP_EXPORT Expr dsolve_cauchy_euler(
    const std::vector<Expr>& coeffs, const Expr& x);

// --- Systems ---
// y' = A · y. Returns a column vector (Matrix n×1) representing the
// general solution. A must be diagonalizable in the current expression
// domain.
[[nodiscard]] SYMPP_EXPORT Matrix dsolve_system(const Matrix& A, const Expr& x);

// --- IVP / verification ---
// Replace each free constant with the value that satisfies one initial
// condition. `conditions` is a list of (x_value, y_value) pairs.
[[nodiscard]] SYMPP_EXPORT Expr apply_ivp(
    const Expr& general_solution, const Expr& x,
    const std::vector<std::pair<Expr, Expr>>& conditions);

// Verify that `sol` (with derivative computed via diff) satisfies `eq`.
// Returns the residual after substitution + simplify; should be 0.
[[nodiscard]] SYMPP_EXPORT Expr checkodesol(
    const Expr& eq, const Expr& sol, const Expr& y, const Expr& yp,
    const Expr& x);

// --- PDE ---
// Solve a(x,y) u_x + b(x,y) u_y = c(x,y) via method of characteristics
// for the simplest case where a, b, c are constants (a u_x + b u_y = c).
// Returns the general solution u(x, y) = (c/a) x + F(b x - a y) where F
// is an arbitrary function (left as `function_symbol("F")`).
[[nodiscard]] SYMPP_EXPORT Expr pdsolve_first_order_linear(
    const Expr& a, const Expr& b, const Expr& c,
    const Expr& x, const Expr& y);

}  // namespace sympp
