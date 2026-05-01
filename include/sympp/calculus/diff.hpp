#pragma once

// Symbolic differentiation.
//
//     diff(expr, var)         — first derivative w.r.t. `var` (a Symbol Expr)
//     diff(expr, var, n)      — n-th derivative
//
// Rules:
//   * d/dx of a Symbol y where y == x  -> 1
//   * d/dx of any other atomic         -> 0
//   * Sum rule:    d/dx (a + b + c)    = da/dx + db/dx + dc/dx
//   * Product rule:d/dx (a * b * c)    = sum_i (prod_{j != i} arg_j) * d/dx arg_i
//   * Power rule:  d/dx (a^b)          = a^b * (d/dx b * log(a) + b * d/dx a / a)
//                  for the common case of constant exponent reduces to b*a^(b-1)*da/dx.
//   * Chain rule on Function:
//       d/dx f(g1, g2, ...)  =  sum_i  partial_i_f(g1, g2, ...) * d/dx gi
//     where partial_i_f is queried via the function's diff_arg() table.
//
// Reference: sympy/core/expr.py::diff and sympy/core/function.py::Function.fdiff

#include <cstddef>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr diff(const Expr& e, const Expr& var);
[[nodiscard]] SYMPP_EXPORT Expr diff(const Expr& e, const Expr& var, std::size_t order);

}  // namespace sympp
