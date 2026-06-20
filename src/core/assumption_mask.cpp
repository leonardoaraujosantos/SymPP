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
        case AssumptionKey::Prime: return prime;
        case AssumptionKey::Composite: return composite;
        case AssumptionKey::Irrational: return irrational;
        case AssumptionKey::Algebraic: return algebraic;
        case AssumptionKey::Transcendental: return transcendental;
        case AssumptionKey::ExtendedReal: return extended_real;
        case AssumptionKey::Infinite: return infinite;
        case AssumptionKey::Commutative: return commutative;
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
        case AssumptionKey::Prime: prime = value; break;
        case AssumptionKey::Composite: composite = value; break;
        case AssumptionKey::Irrational: irrational = value; break;
        case AssumptionKey::Algebraic: algebraic = value; break;
        case AssumptionKey::Transcendental: transcendental = value; break;
        case AssumptionKey::ExtendedReal: extended_real = value; break;
        case AssumptionKey::Infinite: infinite = value; break;
        case AssumptionKey::Commutative: commutative = value; break;
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
           && !complex_ && !imaginary && !prime && !composite && !irrational
           && !algebraic && !transcendental && !extended_real && !infinite
           && !commutative;
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
    mix(encode(prime));
    mix(encode(composite));
    mix(encode(irrational));
    mix(encode(algebraic));
    mix(encode(transcendental));
    mix(encode(extended_real));
    mix(encode(infinite));
    mix(encode(commutative));
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
        // ¬integer => ¬even, ¬odd, ¬prime, ¬composite
        if (m.integer == false) {
            if (!m.even) m.even = false;
            if (!m.odd) m.odd = false;
            if (!m.prime) m.prime = false;
            if (!m.composite) m.composite = false;
        }
        // prime => integer, positive (≥ 2), ¬composite. The positive/integer
        // closures above then cascade to real, finite, nonzero, nonnegative,
        // rational. Note a prime is NOT necessarily odd (2 is prime and even),
        // so no parity rule.
        if (m.prime == true) {
            if (!m.integer) m.integer = true;
            if (!m.positive) m.positive = true;
            if (!m.composite) m.composite = false;
        }
        // composite => integer, positive (≥ 4), ¬prime. Same cascade as prime;
        // again no parity rule (4 is even, 9 is odd). Mutually exclusive with
        // prime; the integer value 1 is neither prime nor composite.
        if (m.composite == true) {
            if (!m.integer) m.integer = true;
            if (!m.positive) m.positive = true;
            if (!m.prime) m.prime = false;
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
        // irrational ⟺ real ∧ ¬rational. Forward: irrational => real, ¬rational,
        // finite (the ¬rational rule above then gives ¬integer/¬zero, hence
        // nonzero; real gives complex/¬imaginary). Reverse / exclusions: a known
        // real non-rational is irrational; a rational or non-real value is not.
        if (m.irrational == true) {
            if (!m.real) m.real = true;
            if (!m.rational) m.rational = false;
            if (!m.finite) m.finite = true;
        }
        if (m.rational == true && !m.irrational) m.irrational = false;
        if (m.real == false && !m.irrational) m.irrational = false;
        if (m.real == true && m.rational == false && !m.irrational) {
            m.irrational = true;
        }
        // algebraic / transcendental — complex-domain predicates. rational ⇒
        // algebraic; algebraic ⇒ complex, finite, ¬transcendental; transcendental
        // ⇒ complex, finite, ¬algebraic, ¬rational (so ¬integer/¬zero via the
        // rational rule, and irrational once real is known). Within ℂ the two
        // partition: complex ∧ ¬algebraic ⇒ transcendental, and complex ∧
        // ¬transcendental ⇒ algebraic. Neither implies real (i is algebraic).
        if (m.rational == true && !m.algebraic) m.algebraic = true;
        if (m.algebraic == true) {
            if (!m.complex_) m.complex_ = true;
            if (!m.finite) m.finite = true;
            if (!m.transcendental) m.transcendental = false;
        }
        if (m.transcendental == true) {
            if (!m.complex_) m.complex_ = true;
            if (!m.finite) m.finite = true;
            if (!m.algebraic) m.algebraic = false;
            if (!m.rational) m.rational = false;
        }
        if (m.complex_ == true && m.algebraic == false && !m.transcendental) {
            m.transcendental = true;
        }
        if (m.complex_ == true && m.transcendental == false && !m.algebraic) {
            m.algebraic = true;
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
        // complex => finite. AssumptionKey::Complex denotes a *finite* complex
        // number (see the enum docs), matching SymPy's is_complex ⇒ is_finite;
        // the unbounded ∞ atoms answer Complex=false directly.
        if (m.complex_ == true && !m.finite) m.finite = true;
        // zero => ¬imaginary (0 is real).
        if (m.zero == true) {
            if (!m.imaginary) m.imaginary = false;
        }
        // extended_real: a point of ℝ ∪ {±∞}. real ⇒ extended_real; a nonzero
        // pure imaginary is off the real line. extended_real does NOT imply real,
        // finite or complex (±∞ are extended-real but neither), but it does
        // exclude imaginary.
        if (m.real == true && !m.extended_real) m.extended_real = true;
        if (m.imaginary == true && !m.extended_real) m.extended_real = false;
        if (m.extended_real == true && !m.imaginary) m.imaginary = false;
        // infinite ⟺ ¬finite. An infinite quantity is not a finite real/complex
        // number and is nonzero (zero is finite).
        if (m.infinite == true) {
            if (!m.finite) m.finite = false;
            if (!m.real) m.real = false;
            if (!m.complex_) m.complex_ = false;
            if (!m.zero) m.zero = false;
        }
        if (m.finite == true && !m.infinite) m.infinite = false;
        if (m.finite == false && !m.infinite) m.infinite = true;
    }
    return m;
}

namespace {
// std::optional<bool> shorthands for the consistency checks below.
[[nodiscard]] inline bool yes(std::optional<bool> v) noexcept { return v == true; }
[[nodiscard]] inline bool no(std::optional<bool> v) noexcept { return v == false; }
}  // namespace

bool assumptions_consistent(AssumptionMask m) noexcept {
    m = close_assumptions(m);

    // Mutually exclusive facts — both true is a contradiction.
    const bool clash =
        (yes(m.positive) && yes(m.negative)) ||
        (yes(m.positive) && yes(m.zero)) ||
        (yes(m.negative) && yes(m.zero)) ||
        (yes(m.even) && yes(m.odd)) ||
        (yes(m.prime) && yes(m.composite)) ||
        (yes(m.real) && yes(m.imaginary)) ||
        (yes(m.finite) && yes(m.infinite)) ||
        (yes(m.rational) && yes(m.irrational)) ||
        (yes(m.algebraic) && yes(m.transcendental)) ||
        (yes(m.nonnegative) && yes(m.negative)) ||
        (yes(m.nonpositive) && yes(m.positive)) ||
        (yes(m.zero) && yes(m.odd)) ||
        (yes(m.zero) && yes(m.imaginary)) ||
        (yes(m.zero) && yes(m.irrational)) ||
        (yes(m.zero) && yes(m.prime)) ||
        (yes(m.zero) && yes(m.composite));
    if (clash) return false;

    // A nonzero pure imaginary excludes every real-line / integer property.
    if (yes(m.imaginary) &&
        (yes(m.real) || yes(m.rational) || yes(m.integer) || yes(m.positive) ||
         yes(m.negative) || yes(m.even) || yes(m.odd) || yes(m.prime) ||
         yes(m.composite))) {
        return false;
    }
    // An infinite quantity is no finite real/complex number (±∞ keep only their
    // sign / extended-real / nonzero facts).
    if (yes(m.infinite) &&
        (yes(m.real) || yes(m.complex_) || yes(m.finite) || yes(m.zero) ||
         yes(m.rational) || yes(m.integer) || yes(m.even) || yes(m.odd) ||
         yes(m.prime) || yes(m.composite) || yes(m.algebraic) ||
         yes(m.transcendental) || yes(m.irrational) || yes(m.imaginary))) {
        return false;
    }
    // Order / tower properties presuppose membership.
    if (no(m.real) &&
        (yes(m.positive) || yes(m.negative) || yes(m.zero) || yes(m.rational) ||
         yes(m.integer) || yes(m.irrational))) {
        return false;
    }
    if (no(m.integer) &&
        (yes(m.even) || yes(m.odd) || yes(m.prime) || yes(m.composite))) {
        return false;
    }
    if (no(m.rational) && (yes(m.integer) || yes(m.zero))) return false;
    if (no(m.complex_) &&
        (yes(m.real) || yes(m.imaginary) || yes(m.integer) || yes(m.rational) ||
         yes(m.algebraic) || yes(m.transcendental))) {
        return false;
    }

    // Completeness laws — failing to be any branch of a total partition is a
    // contradiction (these are what make exhaustive disjunctions provable).
    if (yes(m.real) && no(m.positive) && no(m.negative) && no(m.zero)) return false;
    if (yes(m.real) && no(m.nonnegative) && no(m.nonpositive)) return false;
    if (yes(m.integer) && no(m.even) && no(m.odd)) return false;
    if (yes(m.complex_) && no(m.algebraic) && no(m.transcendental)) return false;

    return true;
}

}  // namespace sympp
