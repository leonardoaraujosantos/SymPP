#pragma once

// solve / linsolve.
//
//   solve(expr, var) -> roots of `expr = 0` (treats expr as the LHS with
//                       implicit RHS = 0).
//   solve(lhs, rhs, var) -> roots of `lhs - rhs = 0`.
//   linsolve(A, b) -> solution of A·x = b for square A (uses Matrix::inverse).
//
// Phase 10 minimal: only polynomial equations of degree <= 2 in `var` are
// handled in closed form. Degree >= 3 / transcendental returns empty
// (deferred — full SymPy solveset is its own phase).
//
// Reference: sympy/solvers/solveset.py and sympy/solvers/solvers.py

#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>
#include <sympp/matrices/matrix.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT std::vector<Expr> solve(const Expr& expr, const Expr& var);
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> solve(const Expr& lhs, const Expr& rhs,
                                                    const Expr& var);

// Solve A·x = b, returning x as a column matrix. Throws on singular A.
[[nodiscard]] SYMPP_EXPORT Matrix linsolve(const Matrix& A, const Matrix& b);

}  // namespace sympp
