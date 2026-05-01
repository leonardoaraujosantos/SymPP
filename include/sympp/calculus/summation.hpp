#pragma once

// summation — closed-form symbolic sums.
//
// summation(expr, var, lo, hi) computes Σ_{var=lo}^{hi} expr when a closed
// form is recognized. Coverage:
//   * Linearity over Add and Mul-with-constant-factors
//   * Constant w.r.t. var: c → c * (hi - lo + 1)
//   * Identity (var): linear sum (hi - lo + 1)(lo + hi)/2
//   * Squares (var^2): n(n+1)(2n+1)/6 telescoped to [lo, hi]
//   * Cubes (var^3): (n(n+1)/2)² telescoped
//   * Geometric (r^var, r free of var): (r^lo - r^(hi+1))/(1 - r)
//
// Returns the input unchanged when no closed form is recognized — the
// caller can layer Gosper/hypergeometric machinery on top in a future
// phase.
//
// Reference: sympy/concrete/summations.py::Sum / .doit

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr summation(const Expr& expr, const Expr& var,
                                            const Expr& lo, const Expr& hi);

}  // namespace sympp
