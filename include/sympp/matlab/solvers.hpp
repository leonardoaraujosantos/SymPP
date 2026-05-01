#pragma once

// MATLAB facade — solvers (multi-eq solve, linsolve, nsolve, vpasolve).
//
// Existing matlab.hpp covers the single-equation univariate
// `solve(eq, var)`. This header adds:
//
//   * solve(eqs, vars)   — multi-equation polynomial system.
//   * linsolve(A, b)     — A·x = b for a square symbolic matrix.
//   * nsolve(eq, var, x0, dps)   — Newton's method in MPFR.
//   * vpasolve(eq, var, x0, dps) — alias for nsolve at higher default
//                                     precision (32 dps), matching
//                                     MATLAB's vpasolve naming.
//
// The multi-eq solve overload picks `nonlinsolve` for the 2×2 case
// (uses resultant elimination — fewer extraneous candidates) and
// `nonlinsolve_groebner` for n ≥ 3 (Buchberger lex Gröbner basis +
// triangular back-substitution).
//
// Reference: sympp/solvers/solve.hpp.

#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/solvers/solve.hpp>

namespace sympp::matlab {

// linsolve(A, b) — solve A·x = b for a square symbolic matrix A.
// Throws on singular A.
[[nodiscard]] inline Matrix linsolve(const Matrix& A, const Matrix& b) {
    return sympp::linsolve(A, b);
}

// solve(eqs, vars) — joint roots of a polynomial system, returned as
// a vector of solution vectors (one entry per var, in the same order
// as `vars`). For 2 eqs / 2 vars uses the resultant-based path
// (filters extraneous candidates numerically); for ≥ 3 falls through
// to Gröbner triangularization.
[[nodiscard]] inline std::vector<std::vector<Expr>> solve(
    const std::vector<Expr>& eqs, const std::vector<Expr>& vars) {
    if (eqs.size() == 2 && vars.size() == 2) {
        return nonlinsolve(eqs, vars);
    }
    return nonlinsolve_groebner(eqs, vars);
}

// nsolve(eq, var, x0, dps) — Newton's method in MPFR for a real
// root near x0. Throws on derivative-vanish or non-convergence.
[[nodiscard]] inline Expr nsolve(const Expr& eq, const Expr& var,
                                    const Expr& x0, int dps = 15) {
    return sympp::nsolve(eq, var, x0, dps);
}

// vpasolve(eq, var, x0, dps) — variable-precision arithmetic solve;
// MATLAB's name for high-precision Newton. Defaults to 32 dps
// (matching MATLAB vpa default).
[[nodiscard]] inline Expr vpasolve(const Expr& eq, const Expr& var,
                                      const Expr& x0, int dps = 32) {
    return sympp::nsolve(eq, var, x0, dps);
}

}  // namespace sympp::matlab
