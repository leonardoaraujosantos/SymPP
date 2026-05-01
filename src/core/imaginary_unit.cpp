#include <sympp/core/imaginary_unit.hpp>

#include <sympp/core/assumption_key.hpp>

namespace sympp {

std::size_t ImaginaryUnit::hash() const noexcept {
    // Stable per-program-run constant. Distinct from numeric hashes.
    return 0xC0FFEE'1A'1A'7E'1AULL;
}

std::optional<bool> ImaginaryUnit::ask(AssumptionKey k) const noexcept {
    // I is finite and non-zero, but not real / integer / rational / signed.
    switch (k) {
        case AssumptionKey::Finite:
        case AssumptionKey::Nonzero:
            return true;
        case AssumptionKey::Real:
        case AssumptionKey::Integer:
        case AssumptionKey::Rational:
        case AssumptionKey::Positive:
        case AssumptionKey::Negative:
        case AssumptionKey::Nonnegative:
        case AssumptionKey::Nonpositive:
        case AssumptionKey::Zero:
            return false;
    }
    return std::nullopt;
}

}  // namespace sympp
