#pragma once

// pade — Padé approximant.
//
// pade(expr, var, m, n) returns the Padé approximant of order [m/n] for
// `expr` at var = 0: a rational function P(x)/Q(x) with deg(P) = m,
// deg(Q) = n, normalized to Q(0) = 1, such that
//
//     expr - P/Q  =  O(var^(m+n+1))
//
// Computed by extracting Taylor coefficients up to order m+n and solving
// the standard n×n linear system for Q's coefficients (we have linsolve
// from Phase 9). P's coefficients fall out via the lower-triangular
// system once Q is known.
//
// Returns the input unchanged if the denominator system is singular —
// e.g. when the Padé table degenerates and a different (m, n) ordering
// would succeed instead.
//
// Reference: sympy/series/approximants.py — alternative formulation;
//            classical: Baker & Graves-Morris, "Padé Approximants".

#include <cstddef>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr pade(const Expr& expr, const Expr& var,
                                      std::size_t m, std::size_t n);

}  // namespace sympp
