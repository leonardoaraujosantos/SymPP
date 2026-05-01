#pragma once

// S — central registry of canonical singleton expressions, modeled on
// SymPy's `sympy.S`. Every singleton returns the same shared_ptr-equal
// instance, so identity-comparison via .get() is meaningful.
//
// Reference: sympy/core/singleton.py and sympy/core/numbers.py

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

namespace S {

// Number singletons (cached Integer / Rational instances).
[[nodiscard]] SYMPP_EXPORT const Expr& Zero();
[[nodiscard]] SYMPP_EXPORT const Expr& One();
[[nodiscard]] SYMPP_EXPORT const Expr& NegativeOne();
[[nodiscard]] SYMPP_EXPORT const Expr& Half();
[[nodiscard]] SYMPP_EXPORT const Expr& NegativeHalf();

// Mathematical constants (NumberSymbol instances).
[[nodiscard]] SYMPP_EXPORT const Expr& Pi();
[[nodiscard]] SYMPP_EXPORT const Expr& E();
[[nodiscard]] SYMPP_EXPORT const Expr& EulerGamma();
[[nodiscard]] SYMPP_EXPORT const Expr& Catalan();

// Boolean atoms.
[[nodiscard]] SYMPP_EXPORT const Expr& True();
[[nodiscard]] SYMPP_EXPORT const Expr& False();

// Imaginary unit. I*I = -1; folded by Mul/Pow.
[[nodiscard]] SYMPP_EXPORT const Expr& I();

}  // namespace S

}  // namespace sympp
