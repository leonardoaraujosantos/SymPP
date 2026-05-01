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
};

}  // namespace sympp
