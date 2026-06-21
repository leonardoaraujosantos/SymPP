#pragma once

// Berlekamp–Zassenhaus univariate factorization over ℤ: factor mod a prime
// (Berlekamp), Hensel-lift to a high prime power, then recombine. This is the
// polynomial-time replacement for the exponential Kronecker method; it is
// exposed separately so the established factor() path is untouched.
//
//   factor_zassenhaus(x^4 - 1, x)  // → [x - 1, x + 1, x^2 + 1]
//
// Reference: sympy/polys/factortools.py — dup_zz_factor / Hensel lifting.

#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// Irreducible factors over ℤ of a univariate integer polynomial `p` in `var`
// (content and constant factors dropped; the unit sign is folded into a factor).
// Throws std::invalid_argument if `p` is not a univariate integer polynomial.
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> factor_zassenhaus(const Expr& p,
                                                               const Expr& var);

}  // namespace sympp
