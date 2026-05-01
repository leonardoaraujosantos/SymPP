#pragma once

// Equation solvers.
//
// Univariate polynomial:
//   solve(expr, var)            roots of expr = 0
//   solve(lhs, rhs, var)        roots of lhs - rhs = 0
//
// Coverage in solve():
//   * deg 1, 2          — closed form (linear / quadratic formula)
//   * deg 3             — Cardano via depressed cubic
//   * deg 4             — Ferrari via resolvent cubic
//   * deg ≥ 5           — rational-roots theorem deflation; cofactor
//                         recurses (so a quintic with one rational
//                         root and a quartic remainder solves fully).
//                         No general radical formula exists by
//                         Abel-Ruffini for irreducible higher degree.
//
// Set-returning variant:
//   solveset(expr, var [, domain]) returns a Set object. For trig
//   patterns (sin(a·x), cos(a·x), tan(a·x)) it emits an ImageSet
//   over Integers representing the full periodic family. Filters
//   roots against the supplied domain when membership is decidable.
//
// Numeric:
//   nsolve(expr, var, x0, dps)  Newton's method in MPFR for a real
//                                root near x0.
//
// Linear algebra:
//   linsolve(A, b)              A·x = b for square A via Matrix::inverse.
//
// Multivariate polynomial systems:
//   nonlinsolve(eqs, vars)            two equations, two variables.
//                                       Eliminates via resultant.
//   nonlinsolve_groebner(eqs, vars)   n equations, n variables. Builds
//                                       a Gröbner basis (Buchberger,
//                                       lex order) and back-solves.
//
// Inequalities:
//   solve_univariate_inequality(lhs, op, rhs, var)
//   reduce_inequalities(rel, var) / reduce_inequalities(rels, var, conj)
//
// Recurrences:
//   rsolve(coeffs, n)           constant-coefficient linear recurrence.
//
// Diophantine:
//   linear_diophantine(a, b, c) ax + by = c over integers.
//   pythagorean_triples(max_z)  primitive triples via Euclid's formula.
//
// Reference: sympy/solvers/{solveset,solvers,polysys,inequalities,recurr}.py.

#include <array>
#include <optional>
#include <utility>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/sets/sets.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT std::vector<Expr> solve(const Expr& expr, const Expr& var);
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> solve(const Expr& lhs, const Expr& rhs,
                                                    const Expr& var);

// Solve A·x = b, returning x as a column matrix. Throws on singular A.
[[nodiscard]] SYMPP_EXPORT Matrix linsolve(const Matrix& A, const Matrix& b);

// solveset(expr, var, domain) — solveset wrapper that returns a Set
// rather than a vector. Empty result becomes EmptySet; otherwise
// returns FiniteSet of the roots solve() found. Falls through to a
// ConditionSet placeholder (currently the literal solve-output
// FiniteSet, since we don't have ConditionSet) when solve() can't.
[[nodiscard]] SYMPP_EXPORT SetPtr solveset(const Expr& expr, const Expr& var);
[[nodiscard]] SYMPP_EXPORT SetPtr solveset(const Expr& expr, const Expr& var,
                                              const SetPtr& domain);

// nsolve(expr, var, x0, dps) — Newton's method in MPFR for finding a
// real root of expr=0 near initial guess x0. Returns a Float Expr at
// the requested precision. Throws std::runtime_error on non-convergence.
[[nodiscard]] SYMPP_EXPORT Expr nsolve(const Expr& expr, const Expr& var,
                                          const Expr& x0, int dps = 15);

// solve_univariate_inequality(lhs, rel, rhs, var) — solution set of
// `lhs rel rhs` (rel ∈ {<, ≤, >, ≥, ≠}) over the reals. Computes
// roots of (lhs - rhs), sorts them, tests the sign on each interval
// between consecutive roots, and assembles a Union of Intervals.
//
// Use the relation kind enum from sympp/core/boolean.hpp.
[[nodiscard]] SYMPP_EXPORT SetPtr solve_univariate_inequality(
    const Expr& lhs, const Expr& rel_op_str, const Expr& rhs, const Expr& var);

// rsolve(coeffs, n) — closed-form solution of the linear constant-
// coefficient recurrence with coefficient list (lowest-index first):
//
//     c_k * y(n+k) + c_{k-1} * y(n+k-1) + ... + c_0 * y(n) = 0
//
// Returns Σ_i C_i * r_i^n where r_i are the roots of the characteristic
// polynomial Σ c_i x^i. Free constants are named __C0, __C1, ...
//
// Reference: sympy/solvers/recurr.py
[[nodiscard]] SYMPP_EXPORT Expr rsolve(const std::vector<Expr>& coeffs,
                                          const Expr& n);

// nonlinsolve(eqs, vars) — for two polynomial equations in two
// variables, eliminates one via resultant (using the Phase 4
// resultant), solves the resulting univariate, and back-substitutes.
// Returns a list of (val_x, val_y) pairs for the joint solutions.
[[nodiscard]] SYMPP_EXPORT std::vector<std::vector<Expr>> nonlinsolve(
    const std::vector<Expr>& eqs, const std::vector<Expr>& vars);

// linear_diophantine(a, b, c) — find integer solutions to a*x + b*y = c.
// Returns std::nullopt if gcd(a, b) does not divide c. Otherwise
// returns (x0, y0) such that a*x0 + b*y0 = c — the general solution
// being x = x0 + b/g * t, y = y0 - a/g * t for integer t.
[[nodiscard]] SYMPP_EXPORT std::optional<std::pair<Expr, Expr>>
linear_diophantine(const Expr& a, const Expr& b, const Expr& c);

// pythagorean_triples(max_z) — primitive Pythagorean triples (x, y, z)
// with z ≤ max_z, parametrized by Euclid's formula:
//   x = m² − n², y = 2·m·n, z = m² + n²
// for m > n > 0 with gcd(m, n) = 1 and opposite parity. Returns the
// list of (x, y, z) triples; up to swap of x, y per triple (one even
// leg, one odd).
[[nodiscard]] SYMPP_EXPORT std::vector<std::array<Expr, 3>>
pythagorean_triples(long max_z);

// reduce_inequalities(rel, var) — solve a single Relational over the
// reals. Thin wrapper that decodes the relation kind and dispatches to
// solve_univariate_inequality.
[[nodiscard]] SYMPP_EXPORT SetPtr reduce_inequalities(
    const Expr& rel, const Expr& var);

// reduce_inequalities(rels, var, conjunction) — combine multiple
// Relational expressions. When conjunction == true, intersects the
// per-leaf solution sets (logical AND); when false, unions them
// (logical OR). Useful for compound bounds like x > 0 AND x < 5.
[[nodiscard]] SYMPP_EXPORT SetPtr reduce_inequalities(
    const std::vector<Expr>& rels, const Expr& var, bool conjunction);

// groebner(eqs, vars) — Buchberger's algorithm with lexicographic
// monomial ordering on multivariate polynomials over ℚ. Returns the
// reduced Gröbner basis as a list of polynomials. Useful for
// triangularizing polynomial systems before back-solving.
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> groebner(
    const std::vector<Expr>& eqs, const std::vector<Expr>& vars);

// nonlinsolve_groebner(eqs, vars) — n-variable polynomial system
// solver via Gröbner triangularization. Returns a list of
// solution-vectors (one per joint root). Generalizes the basic
// 2×2 nonlinsolve.
[[nodiscard]] SYMPP_EXPORT std::vector<std::vector<Expr>>
nonlinsolve_groebner(const std::vector<Expr>& eqs,
                       const std::vector<Expr>& vars);

}  // namespace sympp
