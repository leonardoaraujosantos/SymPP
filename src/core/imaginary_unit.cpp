#include <sympp/core/imaginary_unit.hpp>

#include <sympp/core/assumption_key.hpp>

namespace sympp {

std::size_t ImaginaryUnit::hash() const noexcept {
    // Stable per-program-run constant. Distinct from numeric hashes.
    return 0xC0FFEE'1A'1A'7E'1AULL;
}

std::optional<bool> ImaginaryUnit::ask(AssumptionKey k) const noexcept {
    // I is a finite, non-zero, purely-imaginary complex number — not real /
    // integer / rational / signed.
    switch (k) {
        case AssumptionKey::ExtendedPositive:
        case AssumptionKey::ExtendedNegative:
        case AssumptionKey::ExtendedNonnegative:
        case AssumptionKey::ExtendedNonpositive:
        case AssumptionKey::Hermitian:
        case AssumptionKey::Antihermitian:
        case AssumptionKey::Commutative:
            return std::nullopt;  // derived by the generic ask() layer
        case AssumptionKey::Finite:
        case AssumptionKey::Nonzero:
        case AssumptionKey::Imaginary:
        case AssumptionKey::Complex:
        case AssumptionKey::Algebraic:  // i is a root of x² + 1
            return true;
        case AssumptionKey::Real:
        case AssumptionKey::Integer:
        case AssumptionKey::Rational:
        case AssumptionKey::Positive:
        case AssumptionKey::Negative:
        case AssumptionKey::Nonnegative:
        case AssumptionKey::Nonpositive:
        case AssumptionKey::Zero:
        case AssumptionKey::Even:
        case AssumptionKey::Odd:
        case AssumptionKey::Prime:
        case AssumptionKey::Composite:
        case AssumptionKey::Irrational:  // I is not real, so not irrational
        case AssumptionKey::Transcendental:  // i is algebraic
        case AssumptionKey::ExtendedReal:  // i is off the real line
        case AssumptionKey::Infinite:      // i is finite
            return false;
    }
    return std::nullopt;
}

}  // namespace sympp
