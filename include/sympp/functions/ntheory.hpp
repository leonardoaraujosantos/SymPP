#pragma once

// Number-theory primitives over arbitrary-precision integers (GMP-backed).
// These complement the combinatorial number functions (totient, mobius,
// divisor_count/sigma, primepi, fibonacci, …) already in <functions/
// combinatorial.hpp>.
//
// Reference: sympy/ntheory/{factor_,residue_ntheory,modular}.py and
//            sympy/solvers/diophantine — factorint, divisors, igcdex,
//            jacobi_symbol, crt, discrete_log, diop_linear.

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

// Chinese Remainder Theorem. Solve x ≡ residues[i] (mod moduli[i]) by pairwise
// combination. Returns (x, M) with M = lcm(moduli) and 0 ≤ x < M, or
// std::nullopt when the congruences are inconsistent. Moduli (all > 0) need not
// be pairwise coprime; sizes must match.
[[nodiscard]] SYMPP_EXPORT std::optional<std::pair<Expr, Expr>> crt(
    const std::vector<Expr>& moduli, const std::vector<Expr>& residues);

// Discrete logarithm: the least non-negative x with bˣ ≡ a (mod n), or
// std::nullopt when no such x exists. Baby-step/giant-step, O(√n).
[[nodiscard]] SYMPP_EXPORT std::optional<Expr> discrete_log(const Expr& n, const Expr& a,
                                                            const Expr& b);

// Linear Diophantine equation a·x + b·y = c over ℤ. Parametric solution family
// x = x0 + dx·t, y = y0 + dy·t for integer t.
struct DiophantineLinear {
    Expr x0, y0, dx, dy;
};
// Returns the parametric family, or std::nullopt when gcd(a, b) ∤ c (and so
// there is no integer solution). The a = b = 0 case requires c = 0.
[[nodiscard]] SYMPP_EXPORT std::optional<DiophantineLinear> diop_linear(const Expr& a,
                                                                        const Expr& b,
                                                                        const Expr& c);

// Fundamental solution (x, y) with x, y > 0 of the Pell equation x² − D·y² = 1,
// found from the continued-fraction expansion of √D. std::nullopt when D ≤ 0 or
// D is a perfect square (only the trivial (±1, 0)).
[[nodiscard]] SYMPP_EXPORT std::optional<std::pair<Expr, Expr>> diop_pell(const Expr& D);

// Represent n as a² + b² with 0 ≤ a ≤ b, or std::nullopt when impossible (some
// prime ≡ 3 (mod 4) divides n to an odd power).
[[nodiscard]] SYMPP_EXPORT std::optional<std::pair<Expr, Expr>> sum_of_two_squares(
    const Expr& n);

// Represent n ≥ 0 as a² + b² + c² with 0 ≤ a ≤ b ≤ c (a ternary quadratic form),
// or std::nullopt when impossible — by Legendre's three-square theorem exactly
// when n = 4ᵃ·(8b+7).
[[nodiscard]] SYMPP_EXPORT std::optional<std::vector<Expr>> sum_of_three_squares(
    const Expr& n);

// Represent n ≥ 0 as a² + b² + c² + d² with a ≥ b ≥ c ≥ d ≥ 0 (Lagrange's
// four-square theorem — always possible). std::nullopt only when n < 0.
[[nodiscard]] SYMPP_EXPORT std::optional<std::vector<Expr>> sum_of_four_squares(
    const Expr& n);

}  // namespace sympp
