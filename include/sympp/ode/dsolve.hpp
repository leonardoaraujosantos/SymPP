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

// Riccati: y' = P(x) + Q(x)·y + R(x)·y².
//   Linearizes via y = -u'/(R·u) into the second-order linear ODE
//     R·u'' − (R' + Q·R)·u' + R²·P·u = 0
//   Solve the linear ODE; back-substitute. P, Q, R are passed as
//   Exprs in x (with no y dependence).
[[nodiscard]] SYMPP_EXPORT Expr dsolve_riccati(
    const Expr& P, const Expr& Q, const Expr& R, const Expr& x);

// Homogeneous: y' = f(y/x). Detects when the rhs is a function of
// (y/x) only by substituting y = v·x and checking if x cancels.
// Then v + x·v' = f(v) is separable.
[[nodiscard]] SYMPP_EXPORT Expr dsolve_homogeneous(
    const Expr& eq, const Expr& y, const Expr& yp, const Expr& x);

// Lie point symmetries — minimal subset:
//   * autonomous y' = f(y) where rhs is independent of x — separable
//   * scaling-invariant ODE under x → λx, y → λᵏy for some k.
// Returns the reduced general solution when a symmetry is detected;
// otherwise an unevaluated marker.
[[nodiscard]] SYMPP_EXPORT Expr dsolve_lie_autonomous(
    const Expr& eq, const Expr& y, const Expr& yp, const Expr& x);

// Hypergeometric ODE recognition:
//   x(1-x) y'' + (c - (a+b+1)x) y' - a·b·y = 0  →  ₂F₁(a, b; c; x).
// Coefficients are extracted from the input as polynomials in x.
// Returns hyper(a, b, c, x) when matching, an unevaluated marker
// otherwise.
[[nodiscard]] SYMPP_EXPORT Expr dsolve_hypergeometric(
    const std::vector<Expr>& coeffs_y_ypp,  // [c0_of_y, c1_of_yp, c2_of_ypp]
    const Expr& x);

// hyper(a, b, c, z) — Gauss hypergeometric ₂F₁(a, b; c; z) as an
// opaque function symbol. The full power series implementation
// (hyperexpand) lands in a separate Phase 5 follow-up.
[[nodiscard]] SYMPP_EXPORT Expr hyper(const Expr& a, const Expr& b,
                                        const Expr& c, const Expr& z);

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

// Variable-coefficient method of characteristics: a(x,y)·u_x +
// b(x,y)·u_y = c(x,y). Solves dy/dx = b/a as an ODE for y(x);
// the implicit form ξ(x, y) = const becomes the characteristic.
// Returns u(x, y) = F(ξ(x, y)) on the homogeneous (c=0) case;
// otherwise emits an unevaluated Pdsolve marker (full inhomogeneous
// MOC requires solving for u along characteristics).
[[nodiscard]] SYMPP_EXPORT Expr pdsolve_first_order_variable(
    const Expr& a, const Expr& b, const Expr& c,
    const Expr& y, const Expr& yp, const Expr& x);

// Heat equation u_t = k·u_xx via separation of variables. Returns the
// product solution exp(-k·λ·t) · (A·sin(√λ·x) + B·cos(√λ·x)) for an
// arbitrary separation constant λ (passed as an Expr).
[[nodiscard]] SYMPP_EXPORT Expr pdsolve_heat(
    const Expr& k, const Expr& lambda, const Expr& x, const Expr& t);

// Wave equation u_tt = c²·u_xx via d'Alembert: F(x - c·t) + G(x + c·t).
[[nodiscard]] SYMPP_EXPORT Expr pdsolve_wave(
    const Expr& c, const Expr& x, const Expr& t);

// DAE helpers: Jacobian extraction. Given F(y, y', x) = 0 (a vector of
// equations), returns (M, K) where M = ∂F/∂y' and K = ∂F/∂y.
// A DAE is index-1 when M is non-singular.
[[nodiscard]] SYMPP_EXPORT std::pair<Matrix, Matrix> dae_jacobians(
    const std::vector<Expr>& F, const std::vector<Expr>& y,
    const std::vector<Expr>& yp);

// Structural differential index — coarse estimate. Returns 0 if M is
// already non-singular (purely algebraic), 1 if it's singular but the
// system is otherwise well-posed, and ≥2 with a brief explanation
// when full Pantelides analysis would be needed (deep-deferred).
[[nodiscard]] SYMPP_EXPORT std::size_t dae_structural_index(
    const std::vector<Expr>& F, const std::vector<Expr>& y,
    const std::vector<Expr>& yp);

}  // namespace sympp
