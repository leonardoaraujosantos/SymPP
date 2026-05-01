#pragma once

// Number-on-Number arithmetic with coercion.
//
// Coercion ladder: Integer ⊂ Rational ⊂ Float. NumberSymbol is treated as
// non-numeric here — operations involving it return std::nullopt so callers
// keep them unevaluated as Add/Mul nodes.
//
// Reference: sympy/core/numbers.py — Integer/Rational/Float arithmetic
//            sympy/core/operations.py — AssocOp.flatten

#include <optional>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class Number;

// Returns sum / product / quotient / power if both operands are numeric and
// the result is exactly representable, else std::nullopt. Float results are
// computed at the higher of the two operands' precisions (for Float-on-Float)
// or at default precision when promoting from exact.
[[nodiscard]] SYMPP_EXPORT std::optional<Expr> number_add(const Number& a, const Number& b);
[[nodiscard]] SYMPP_EXPORT std::optional<Expr> number_mul(const Number& a, const Number& b);

// Returns base^exp when (base, exp) are both numeric and the result is exact:
//   * Integer^non-negative-Integer  → Integer
//   * Integer^negative-Integer      → Rational
//   * Rational^Integer              → Rational
//   * Float^Integer                 → Float
//   * x^0                           → 1
//   * x^1                           → x
//   * 1^x                           → 1
//   * 0^positive                    → 0
//
// Returns std::nullopt for cases that should stay unevaluated (e.g. 2^(1/2)).
[[nodiscard]] SYMPP_EXPORT std::optional<Expr> number_pow(const Number& base, const Number& exp);

// True iff the Expr is a Number subclass (Integer, Rational, Float; *not*
// NumberSymbol — Pi/E behave symbolically in arithmetic).
[[nodiscard]] SYMPP_EXPORT bool is_number(const Expr& e) noexcept;

}  // namespace sympp
