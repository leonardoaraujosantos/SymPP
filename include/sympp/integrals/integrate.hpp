#pragma once

// integrate(expr, var) — antiderivative.
// integrate(expr, var, lower, upper) — definite via Newton-Leibniz.
//
// Phase 7 minimal: linearity, basic table.
//
//     ∫ a*f(x) dx       = a * ∫ f(x) dx                  (linearity)
//     ∫ x^n dx          = x^(n+1) / (n+1)   (n != -1)
//     ∫ 1/x dx          = log(x)
//     ∫ exp(x) dx       = exp(x)
//     ∫ sin(x) dx       = -cos(x)
//     ∫ cos(x) dx       = sin(x)
//     ∫ exp(c*x) dx     = exp(c*x)/c       (c constant)
//     ∫ sin(c*x) dx     = -cos(c*x)/c      (c constant)
//     ∫ cos(c*x) dx     = sin(c*x)/c       (c constant)
//
// Anything outside this table returns the input wrapped in a placeholder
// "Integral" UndefinedFunction so downstream code can recognize the
// unevaluated case. Full Risch / heuristic / meijer-g remains deferred.
//
// Reference: sympy/integrals/integrals.py and manualintegrate.py

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr integrate(const Expr& expr, const Expr& var);
[[nodiscard]] SYMPP_EXPORT Expr integrate(const Expr& expr, const Expr& var,
                                            const Expr& lower, const Expr& upper);

}  // namespace sympp
