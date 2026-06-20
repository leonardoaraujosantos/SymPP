#pragma once

// AssumptionKey — the predicate vocabulary used by the assumption query
// system. Each key names a property an Expr may answer True / False / Unknown
// (std::optional<bool> with std::nullopt = Unknown).
//
// Phase 2 ships the practical subset. Less common predicates (algebraic,
// transcendental, hermitian, even, odd, prime, composite, ...) can be added
// without breaking the API.
//
// Reference: sympy/core/assumptions.py — _assume_defined,
//            sympy/assumptions/ask.py — Q.<name>

#include <cstdint>

namespace sympp {

enum class AssumptionKey : std::uint8_t {
    // Numeric tower membership.
    Real,
    Rational,
    Integer,

    // Sign.
    Positive,      // > 0  (implies real, nonzero, nonnegative)
    Negative,      // < 0  (implies real, nonzero, nonpositive)
    Zero,          // == 0 (implies real, integer, rational, nonpositive, nonnegative)
    Nonzero,       // != 0 (= positive | negative for real-known expressions)
    Nonnegative,   // >= 0 (= positive | zero)
    Nonpositive,   // <= 0 (= negative | zero)

    // Range.
    Finite,

    // Parity (integers only).
    Even,          // = 2k for some integer k (implies integer; 0 is even)
    Odd,           // = 2k+1 (implies integer, nonzero)

    // Complex domain.
    Complex,       // a finite complex number (real ⇒ complex, imaginary ⇒ complex)
    Imaginary,     // a nonzero real multiple of i (implies complex, ¬real, finite;
                   // 0 is NOT imaginary — it is real)

    // Number theory (integers only).
    Prime,         // a prime ≥ 2 (implies integer, positive; NOT odd — 2 is prime
                   // and even). ¬integer ⇒ ¬prime.
    Composite,     // a composite ≥ 4: integer > 1 that is not prime (implies
                   // integer, positive, ¬prime; NOT a parity). ¬integer ⇒
                   // ¬composite. Mutually exclusive with Prime; 1 is neither.
    Irrational,    // a real number that is not rational (implies real, finite,
                   // nonzero, ¬integer). irrational ⟺ real ∧ ¬rational; a
                   // rational or non-real value is not irrational.
    Algebraic,     // a root of a nonzero polynomial with rational coefficients
                   // (implies complex, finite). rational ⇒ algebraic. Does NOT
                   // imply real (i is algebraic). Mutually exclusive with
                   // Transcendental.
    Transcendental,// a complex number that is not algebraic (implies complex,
                   // finite, ¬algebraic, ¬rational). transcendental ⟺ complex ∧
                   // ¬algebraic; a real transcendental is irrational.

    // Extended real line / boundedness.
    ExtendedReal,  // a point of ℝ ∪ {−∞, +∞} (implies ¬imaginary). real ⇒
                   // extended_real, but NOT conversely (±∞ are extended-real,
                   // not real). Does not imply finite or complex.
    Infinite,      // an infinite quantity (±∞, zoo). infinite ⟺ ¬finite; implies
                   // ¬real, ¬complex, ¬zero.

    // Algebraic structure.
    Commutative,   // commutes under multiplication (a·b = b·a). DEFAULT TRUE:
                   // every number and ordinary symbol is commutative, and a
                   // compound (Add/Mul/Pow/…) is commutative iff all of its
                   // operands are. Only a symbol explicitly declared
                   // commutative=false (and anything built from it) answers
                   // false. Mirrors SymPy, where commutative defaults to True.
};

}  // namespace sympp
