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
        case AssumptionKey::Even: return even;
        case AssumptionKey::Odd: return odd;
        case AssumptionKey::Imaginary: return imaginary;
        case AssumptionKey::Complex: {
            // real ∨ imaginary ⇒ complex (a finite complex number).
            if (complex_) return complex_;
            if (real == true || imaginary == true) return true;
            return std::nullopt;
        }

        case AssumptionKey::Nonzero: {
            // nonzero = positive ∨ negative ∨ ¬zero (when known).
            if (zero == true) return false;
            if (zero == false) return true;
            if (positive == true || negative == true) return true;
            return std::nullopt;
        }
        case AssumptionKey::Nonnegative: {
            // A directly-declared nonnegative fact takes precedence; otherwise
            // x >= 0  ⇔  positive ∨ zero.
            if (nonnegative) return nonnegative;
            if (positive == true) return true;
            if (zero == true) return true;
            if (negative == true) return false;
            if (positive == false && zero == false) return false;
            return std::nullopt;
        }
        case AssumptionKey::Nonpositive: {
            if (nonpositive) return nonpositive;
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
        case AssumptionKey::Even: even = value; break;
        case AssumptionKey::Odd: odd = value; break;
        case AssumptionKey::Complex: complex_ = value; break;
        case AssumptionKey::Imaginary: imaginary = value; break;
        case AssumptionKey::Nonzero:
            // nonzero=true means zero=false; nonzero=false means zero=true.
            zero = !value;
            break;
        case AssumptionKey::Nonnegative:
            // Store as a primary fact; close_assumptions() derives the rest
            // (real, negative=false, and positive/zero once the other is known).
            nonnegative = value;
            break;
        case AssumptionKey::Nonpositive:
            nonpositive = value;
            break;
    }
}

bool AssumptionMask::empty() const noexcept {
    return !real && !rational && !integer && !positive && !negative && !zero
           && !nonnegative && !nonpositive && !finite && !even && !odd
           && !complex_ && !imaginary;
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
    mix(encode(nonnegative));
    mix(encode(nonpositive));
    mix(encode(finite));
    mix(encode(even));
    mix(encode(odd));
    mix(encode(complex_));
    mix(encode(imaginary));
    return h;
}

AssumptionMask close_assumptions(AssumptionMask m) noexcept {
    // Apply implications until fixpoint. The rules are simple enough that
    // three passes suffice; we run four for safety.
    for (int pass = 0; pass < 4; ++pass) {
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
        // zero => real, integer, rational, finite, positive=false, negative=false,
        // and 0 is even (not odd)
        if (m.zero == true) {
            if (!m.real) m.real = true;
            if (!m.integer) m.integer = true;
            if (!m.rational) m.rational = true;
            if (!m.finite) m.finite = true;
            if (!m.positive) m.positive = false;
            if (!m.negative) m.negative = false;
            if (!m.even) m.even = true;
            if (!m.odd) m.odd = false;
        }
        // nonnegative (x ≥ 0) / nonpositive (x ≤ 0), declarable in their own right.
        // Consistency with the sign primaries: positive ∨ zero ⇒ nonnegative,
        // negative ∨ zero ⇒ nonpositive, and the strict signs exclude the opposite.
        if (m.positive == true && !m.nonnegative) m.nonnegative = true;
        if (m.zero == true && !m.nonnegative) m.nonnegative = true;
        if (m.negative == true && !m.nonpositive) m.nonpositive = true;
        if (m.zero == true && !m.nonpositive) m.nonpositive = true;
        if (m.positive == true && !m.nonpositive) m.nonpositive = false;
        if (m.negative == true && !m.nonnegative) m.nonnegative = false;
        // nonnegative ⇒ real, finite, ¬negative; with the sign pinned it refines
        // to positive (x ≠ 0) or zero (x not positive).
        if (m.nonnegative == true) {
            if (!m.real) m.real = true;
            if (!m.finite) m.finite = true;
            if (!m.negative) m.negative = false;
            if (m.zero == false && !m.positive) m.positive = true;
            if (m.positive == false && !m.zero) m.zero = true;
            if (m.nonpositive == true && !m.zero) m.zero = true;  // ≥0 ∧ ≤0 ⇒ =0
        }
        if (m.nonnegative == false && m.real == true && !m.negative) {
            m.negative = true;  // real ∧ ¬(≥0) ⇒ < 0
        }
        if (m.nonpositive == true) {
            if (!m.real) m.real = true;
            if (!m.finite) m.finite = true;
            if (!m.positive) m.positive = false;
            if (m.zero == false && !m.negative) m.negative = true;
            if (m.negative == false && !m.zero) m.zero = true;
            if (m.nonnegative == true && !m.zero) m.zero = true;
        }
        if (m.nonpositive == false && m.real == true && !m.positive) {
            m.positive = true;  // real ∧ ¬(≤0) ⇒ > 0
        }
        // even => integer; odd => integer, nonzero; even/odd are mutually exclusive;
        // a known integer that is not even is odd (and vice versa).
        if (m.even == true) {
            if (!m.integer) m.integer = true;
            if (!m.odd) m.odd = false;
        }
        if (m.odd == true) {
            if (!m.integer) m.integer = true;
            if (!m.zero) m.zero = false;
            if (!m.even) m.even = false;
        }
        if (m.integer == true) {
            if (m.even == false && !m.odd) m.odd = true;
            if (m.odd == false && !m.even) m.even = true;
        }
        // ¬integer => ¬even, ¬odd
        if (m.integer == false) {
            if (!m.even) m.even = false;
            if (!m.odd) m.odd = false;
        }
        // integer => rational => real
        if (m.integer == true) {
            if (!m.rational) m.rational = true;
            if (!m.real) m.real = true;
        }
        if (m.rational == true) {
            if (!m.real) m.real = true;
        }
        // ¬real implies ¬rational, ¬integer, ¬positive, ¬negative, ¬zero,
        // ¬nonnegative, ¬nonpositive (the order predicates are defined on ℝ).
        if (m.real == false) {
            if (!m.rational) m.rational = false;
            if (!m.integer) m.integer = false;
            if (!m.positive) m.positive = false;
            if (!m.negative) m.negative = false;
            if (!m.zero) m.zero = false;
            if (!m.nonnegative) m.nonnegative = false;
            if (!m.nonpositive) m.nonpositive = false;
        }
        // ¬rational => ¬integer, ¬zero
        if (m.rational == false) {
            if (!m.integer) m.integer = false;
            if (!m.zero) m.zero = false;
        }
        // imaginary (a nonzero real multiple of i) => complex, ¬real, finite,
        // nonzero, and (since ¬real) ¬rational/¬integer/¬sign/¬parity.
        if (m.imaginary == true) {
            if (!m.complex_) m.complex_ = true;
            if (!m.real) m.real = false;
            if (!m.finite) m.finite = true;
            if (!m.zero) m.zero = false;
            if (!m.rational) m.rational = false;
            if (!m.integer) m.integer = false;
            if (!m.positive) m.positive = false;
            if (!m.negative) m.negative = false;
            if (!m.even) m.even = false;
            if (!m.odd) m.odd = false;
        }
        // real => complex, finite, ¬imaginary. A symbol declared `real` denotes a
        // finite real number; the unbounded values ±∞ are the separate Infinity
        // atoms, not real symbols (matching SymPy, where real ⇒ finite and oo is
        // only extended-real). Consistent with positive/negative/zero ⇒ finite.
        if (m.real == true) {
            if (!m.complex_) m.complex_ = true;
            if (!m.finite) m.finite = true;
            if (!m.imaginary) m.imaginary = false;
        }
        // zero => ¬imaginary (0 is real).
        if (m.zero == true) {
            if (!m.imaginary) m.imaginary = false;
        }
    }
    return m;
}

}  // namespace sympp
