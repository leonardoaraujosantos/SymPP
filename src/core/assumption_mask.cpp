#include <sympp/core/assumption_mask.hpp>

#include <cstddef>
#include <functional>

namespace sympp {

std::optional<bool> AssumptionMask::get(AssumptionKey k) const noexcept {
    switch (k) {
        case AssumptionKey::Real: return real;
        case AssumptionKey::Rational: return rational;
        case AssumptionKey::Integer: return integer;
        case AssumptionKey::Positive: return positive;
        case AssumptionKey::Negative: return negative;
        case AssumptionKey::Zero: return zero;
        case AssumptionKey::Finite: return finite;

        case AssumptionKey::Nonzero: {
            // nonzero = positive ∨ negative ∨ ¬zero (when known).
            if (zero == true) return false;
            if (zero == false) return true;
            if (positive == true || negative == true) return true;
            return std::nullopt;
        }
        case AssumptionKey::Nonnegative: {
            // x >= 0  ⇔  positive ∨ zero
            if (positive == true) return true;
            if (zero == true) return true;
            if (negative == true) return false;
            if (positive == false && zero == false) return false;
            return std::nullopt;
        }
        case AssumptionKey::Nonpositive: {
            if (negative == true) return true;
            if (zero == true) return true;
            if (positive == true) return false;
            if (negative == false && zero == false) return false;
            return std::nullopt;
        }
    }
    return std::nullopt;
}

void AssumptionMask::set(AssumptionKey k, bool value) noexcept {
    switch (k) {
        case AssumptionKey::Real: real = value; break;
        case AssumptionKey::Rational: rational = value; break;
        case AssumptionKey::Integer: integer = value; break;
        case AssumptionKey::Positive: positive = value; break;
        case AssumptionKey::Negative: negative = value; break;
        case AssumptionKey::Zero: zero = value; break;
        case AssumptionKey::Finite: finite = value; break;
        case AssumptionKey::Nonzero:
            // nonzero=true means zero=false; nonzero=false means zero=true.
            zero = !value;
            break;
        case AssumptionKey::Nonnegative:
            if (value) {
                negative = false;
            } else {
                negative = true;
                zero = false;
            }
            break;
        case AssumptionKey::Nonpositive:
            if (value) {
                positive = false;
            } else {
                positive = true;
                zero = false;
            }
            break;
    }
}

bool AssumptionMask::empty() const noexcept {
    return !real && !rational && !integer && !positive && !negative && !zero && !finite;
}

std::size_t AssumptionMask::hash() const noexcept {
    auto encode = [](std::optional<bool> v) noexcept -> std::size_t {
        if (!v) return 0;
        return *v ? 1 : 2;
    };
    std::size_t h = 0;
    auto mix = [&h](std::size_t v) noexcept {
        h ^= v + 0x9E37'79B9'7F4A'7C15ULL + (h << 6) + (h >> 2);
    };
    mix(encode(real));
    mix(encode(rational));
    mix(encode(integer));
    mix(encode(positive));
    mix(encode(negative));
    mix(encode(zero));
    mix(encode(finite));
    return h;
}

AssumptionMask close_assumptions(AssumptionMask m) noexcept {
    // Apply implications until fixpoint. The rules are simple enough that
    // two passes suffice; we run three for safety.
    for (int pass = 0; pass < 3; ++pass) {
        // positive => real, finite, nonzero (zero=false), negative=false
        if (m.positive == true) {
            if (!m.real) m.real = true;
            if (!m.finite) m.finite = true;
            if (!m.zero) m.zero = false;
            if (!m.negative) m.negative = false;
        }
        // negative => real, finite, nonzero, positive=false
        if (m.negative == true) {
            if (!m.real) m.real = true;
            if (!m.finite) m.finite = true;
            if (!m.zero) m.zero = false;
            if (!m.positive) m.positive = false;
        }
        // zero => real, integer, rational, finite, positive=false, negative=false
        if (m.zero == true) {
            if (!m.real) m.real = true;
            if (!m.integer) m.integer = true;
            if (!m.rational) m.rational = true;
            if (!m.finite) m.finite = true;
            if (!m.positive) m.positive = false;
            if (!m.negative) m.negative = false;
        }
        // integer => rational => real
        if (m.integer == true) {
            if (!m.rational) m.rational = true;
            if (!m.real) m.real = true;
        }
        if (m.rational == true) {
            if (!m.real) m.real = true;
        }
        // ¬real implies ¬rational, ¬integer, ¬positive, ¬negative, ¬zero
        if (m.real == false) {
            if (!m.rational) m.rational = false;
            if (!m.integer) m.integer = false;
            if (!m.positive) m.positive = false;
            if (!m.negative) m.negative = false;
            if (!m.zero) m.zero = false;
        }
        // ¬rational => ¬integer, ¬zero
        if (m.rational == false) {
            if (!m.integer) m.integer = false;
            if (!m.zero) m.zero = false;
        }
    }
    return m;
}

}  // namespace sympp
