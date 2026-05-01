#pragma once

// Operator overloads for Expr — sugar over the add()/mul()/pow() factories.
//
// Examples:
//   Expr x = symbol("x");
//   auto e = pow(x, integer(2)) + integer(3) * x - integer(7);
//
// All operators short-circuit through the canonical-form factories, so you
// never have to call them directly for common usage.

#include <sympp/core/add.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] inline Expr operator+(const Expr& a, const Expr& b) {
    return add(a, b);
}

[[nodiscard]] inline Expr operator*(const Expr& a, const Expr& b) {
    return mul(a, b);
}

// Unary minus = -1 * a.
[[nodiscard]] inline Expr operator-(const Expr& a) {
    return mul(S::NegativeOne(), a);
}

[[nodiscard]] inline Expr operator-(const Expr& a, const Expr& b) {
    return add(a, mul(S::NegativeOne(), b));
}

// Division = a * b^-1. (Phase 1f.1: only auto-evaluates for numeric b that
// inverts cleanly. Symbolic divisions stay as Mul(a, Pow(b, -1)).)
[[nodiscard]] inline Expr operator/(const Expr& a, const Expr& b) {
    return mul(a, pow(b, S::NegativeOne()));
}

}  // namespace sympp
