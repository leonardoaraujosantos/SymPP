#pragma once

// Numeric definite integration at variable precision (MPFR).
//
// vpaintegral(f, var, lower, upper, dps) — compute ∫_lower^upper f dvar to
// dps decimal digits. Uses adaptive Simpson's rule: split the interval
// recursively when the difference between Simpson's-1/3 on [a,b] and on
// [a,m]+[m,b] exceeds the tolerance derived from dps.
//
// Endpoints `lower` and `upper` may be any numeric Expr (Integer,
// Rational, Float, NumberSymbol). The integrand is materialized to a
// Float at each sample via subs() + evalf().
//
// Returns a Float Expr at the requested dps.
//
// Reference: sympy/integrals/quadrature.py — though SymPy delegates to
// mpmath; ours is a self-contained adaptive Simpson.

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr vpaintegral(const Expr& f, const Expr& var,
                                              const Expr& lower,
                                              const Expr& upper,
                                              int dps = 15);

}  // namespace sympp
