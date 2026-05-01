#pragma once

// Laplace and Fourier transforms (table-based, linearity-aware).
//
// Phase 8 minimal: pattern-matches a small set of well-known transforms
// against the input and returns the table entry. Anything outside the
// table returns an UndefinedFunction marker.
//
// Laplace pairs (transformed variable s):
//     1        ->  1/s
//     t        ->  1/s^2
//     t^n      ->  factorial(n) / s^(n+1)              (n nonneg integer)
//     exp(a*t) ->  1/(s - a)
//     sin(a*t) ->  a / (s^2 + a^2)
//     cos(a*t) ->  s / (s^2 + a^2)
//
// Fourier pairs (with the SymPy "0,-2*pi" convention) — defer to a follow-up;
// only Laplace lands in this batch.
//
// Reference: sympy/integrals/laplace.py

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr laplace_transform(const Expr& f,
                                                    const Expr& t,
                                                    const Expr& s);

}  // namespace sympp
