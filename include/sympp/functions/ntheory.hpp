#pragma once

// Number-theory primitives over arbitrary-precision integers (GMP-backed).
// These complement the combinatorial number functions (totient, mobius,
// divisor_count/sigma, primepi, fibonacci, …) already in <functions/
// combinatorial.hpp>.
//
// Reference: sympy/ntheory/{factor_,residue_ntheory}.py — factorint, divisors,
//            igcdex, jacobi_symbol.

#include <optional>
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

// Simple continued-fraction expansion of a rational (or integer) as a list of
// partial quotients [a0; a1, a2, …] (floor convention, so negatives match
// SymPy). Finite for every rational.
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> continued_fraction(const Expr& x);

// Multiplicative order of a modulo n — the least k > 0 with a^k ≡ 1 (mod n).
// Requires gcd(a, n) = 1 (throws otherwise).
[[nodiscard]] SYMPP_EXPORT Expr n_order(const Expr& a, const Expr& n);

// Smallest primitive root modulo n (n > 1), or std::nullopt when none exists
// (i.e. n is not 2, 4, pᵏ or 2·pᵏ for an odd prime p).
[[nodiscard]] SYMPP_EXPORT std::optional<Expr> primitive_root(const Expr& n);

// A square root of a modulo a prime p — some r with r² ≡ a (mod p) — or
// std::nullopt when a is a non-residue. p must be prime (throws otherwise).
[[nodiscard]] SYMPP_EXPORT std::optional<Expr> sqrt_mod(const Expr& a, const Expr& p);

}  // namespace sympp
