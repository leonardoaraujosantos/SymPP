#pragma once

// Classical orthogonal polynomials, evaluated for a nonnegative integer degree
// via their three-term recurrences and returned as an expanded polynomial in
// the given variable.
//
//   legendre(2, x)    → 3/2·x² − 1/2
//   chebyshevt(3, x)  → 4·x³ − 3·x
//
// Reference: sympy/functions/special/polynomials.py — legendre, chebyshevt,
//            chebyshevu, hermite (physicists'), laguerre.

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// Legendre polynomial Pₙ(x).
[[nodiscard]] SYMPP_EXPORT Expr legendre(const Expr& n, const Expr& x);
// Chebyshev polynomial of the first kind Tₙ(x).
[[nodiscard]] SYMPP_EXPORT Expr chebyshevt(const Expr& n, const Expr& x);
// Chebyshev polynomial of the second kind Uₙ(x).
[[nodiscard]] SYMPP_EXPORT Expr chebyshevu(const Expr& n, const Expr& x);
// Physicists' Hermite polynomial Hₙ(x).
[[nodiscard]] SYMPP_EXPORT Expr hermite(const Expr& n, const Expr& x);
// Laguerre polynomial Lₙ(x).
[[nodiscard]] SYMPP_EXPORT Expr laguerre(const Expr& n, const Expr& x);

}  // namespace sympp
