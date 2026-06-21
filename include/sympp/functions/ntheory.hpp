#pragma once

// Number-theory primitives over arbitrary-precision integers (GMP-backed).
// These complement the combinatorial number functions (totient, mobius,
// divisor_count/sigma, primepi, fibonacci, …) already in <functions/
// combinatorial.hpp>.
//
// Reference: sympy/ntheory/{factor_,residue_ntheory}.py — factorint, divisors,
//            igcdex, jacobi_symbol.

#include <tuple>
#include <utility>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// Prime factorization of |n| as ascending (prime, exponent) pairs. n == 0 is an
// error; ±1 factors to the empty list. Uses trial division for small factors and
// Pollard's rho for larger composite cofactors.
[[nodiscard]] SYMPP_EXPORT std::vector<std::pair<Expr, Expr>> factorint(const Expr& n);

// All positive divisors of |n|, ascending (so 1 and |n| included). |n| == 0
// returns the empty list; ±1 returns {1}.
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> divisors(const Expr& n);

// Extended Euclid: returns (x, y, g) with a·x + b·y = g = gcd(a, b), g ≥ 0.
[[nodiscard]] SYMPP_EXPORT std::tuple<Expr, Expr, Expr> igcdex(const Expr& a,
                                                               const Expr& b);

// Jacobi symbol (a / n) for an odd n > 0; returns −1, 0 or 1.
[[nodiscard]] SYMPP_EXPORT Expr jacobi_symbol(const Expr& a, const Expr& n);

}  // namespace sympp
