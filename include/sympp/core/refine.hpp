#pragma once

// refine — apply assumption-justified rewrites to an expression.
//
// Phase 2c v1 ships a minimal subset:
//
//   * Pow flattening:  (x^a)^b → x^(a*b)
//       valid when x is positive (any real a, b)
//       or when both a and b are integers (any x with conventional 0^0)
//
// More rewrites land alongside the elementary-function module (Phase 3):
//
//   * |x| → x         when x >= 0
//   * |x| → -x        when x <= 0
//   * sqrt(x^2) → x   when x >= 0
//   * sqrt(x^2) → -x  when x <= 0
//   * log(exp(x)) → x   when x is real
//   * exp(log(x)) → x   when x > 0
//
// Reference: sympy/assumptions/refine.py.

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr refine(const Expr& e);

}  // namespace sympp
