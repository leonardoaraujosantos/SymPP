#pragma once

// AssumptionMask — the set of facts a Symbol carries about itself.
//
// Construction:
//     auto x = symbol("x", AssumptionMask{}.set(AssumptionKey::Positive, true));
//     auto y = symbol("y", AssumptionMask{}.set_real(true).set_integer(true));
//
// On construction the mask is "closed" via close_assumptions(): any fact that
// follows by definition (positive => real, zero => integer, etc.) is filled
// in so query lookups don't need rule application at runtime.
//
// Reference: sympy/core/assumptions.py — assumption fact propagation table.

#include <cstddef>
#include <optional>

#include <sympp/core/api.hpp>
#include <sympp/core/assumption_key.hpp>

namespace sympp {

struct SYMPP_EXPORT AssumptionMask {
    std::optional<bool> real;
    std::optional<bool> rational;
    std::optional<bool> integer;
    std::optional<bool> positive;
    std::optional<bool> negative;
    std::optional<bool> zero;
    std::optional<bool> nonnegative;  // x >= 0 (declarable; not just positive ∨ zero)
    std::optional<bool> nonpositive;  // x <= 0
    std::optional<bool> finite;
    std::optional<bool> even;
    std::optional<bool> odd;
    std::optional<bool> complex_;
    std::optional<bool> imaginary;
    std::optional<bool> prime;
    std::optional<bool> composite;
    std::optional<bool> irrational;
    std::optional<bool> algebraic;
    std::optional<bool> transcendental;
    std::optional<bool> extended_real;
    std::optional<bool> infinite;
    std::optional<bool> extended_positive;
    std::optional<bool> extended_negative;
    std::optional<bool> extended_nonnegative;
    std::optional<bool> extended_nonpositive;
    std::optional<bool> hermitian;
    std::optional<bool> antihermitian;
    std::optional<bool> commutative;  // default-true; only stored when declared

    [[nodiscard]] std::optional<bool> get(AssumptionKey k) const noexcept;
    void set(AssumptionKey k, bool value) noexcept;

    // Builder helpers — return *this so they chain.
    AssumptionMask& set_real(bool v) noexcept { real = v; return *this; }
    AssumptionMask& set_rational(bool v) noexcept { rational = v; return *this; }
    AssumptionMask& set_integer(bool v) noexcept { integer = v; return *this; }
    AssumptionMask& set_positive(bool v) noexcept { positive = v; return *this; }
    AssumptionMask& set_negative(bool v) noexcept { negative = v; return *this; }
    AssumptionMask& set_zero(bool v) noexcept { zero = v; return *this; }
    AssumptionMask& set_nonnegative(bool v) noexcept { nonnegative = v; return *this; }
    AssumptionMask& set_nonpositive(bool v) noexcept { nonpositive = v; return *this; }
    AssumptionMask& set_nonzero(bool v) noexcept { zero = !v; return *this; }
    AssumptionMask& set_finite(bool v) noexcept { finite = v; return *this; }
    AssumptionMask& set_even(bool v) noexcept { even = v; return *this; }
    AssumptionMask& set_odd(bool v) noexcept { odd = v; return *this; }
    AssumptionMask& set_complex(bool v) noexcept { complex_ = v; return *this; }
    AssumptionMask& set_imaginary(bool v) noexcept { imaginary = v; return *this; }
    AssumptionMask& set_prime(bool v) noexcept { prime = v; return *this; }
    AssumptionMask& set_composite(bool v) noexcept { composite = v; return *this; }
    AssumptionMask& set_irrational(bool v) noexcept { irrational = v; return *this; }
    AssumptionMask& set_algebraic(bool v) noexcept { algebraic = v; return *this; }
    AssumptionMask& set_transcendental(bool v) noexcept {
        transcendental = v; return *this;
    }
    AssumptionMask& set_extended_real(bool v) noexcept {
        extended_real = v; return *this;
    }
    AssumptionMask& set_infinite(bool v) noexcept { infinite = v; return *this; }
    AssumptionMask& set_extended_positive(bool v) noexcept {
        extended_positive = v; return *this;
    }
    AssumptionMask& set_extended_negative(bool v) noexcept {
        extended_negative = v; return *this;
    }
    AssumptionMask& set_extended_nonnegative(bool v) noexcept {
        extended_nonnegative = v; return *this;
    }
    AssumptionMask& set_extended_nonpositive(bool v) noexcept {
        extended_nonpositive = v; return *this;
    }
    AssumptionMask& set_hermitian(bool v) noexcept { hermitian = v; return *this; }
    AssumptionMask& set_antihermitian(bool v) noexcept { antihermitian = v; return *this; }
    AssumptionMask& set_commutative(bool v) noexcept { commutative = v; return *this; }

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] std::size_t hash() const noexcept;
    [[nodiscard]] bool operator==(const AssumptionMask&) const noexcept = default;
};

// Compute the deductive closure of `m`: fill in all facts that follow from
// the given ones via the standard mathematical implications. If `m` carries
// inconsistent facts (e.g. positive=true and negative=true) the result has
// the conflicting fields set to their explicit input values; downstream code
// is responsible for treating the symbol as ill-formed.
[[nodiscard]] SYMPP_EXPORT AssumptionMask close_assumptions(AssumptionMask m) noexcept;

// True iff `m` (after closure) carries no contradictory pair of facts and
// satisfies the basic completeness laws of the ontology:
//   * a real number is positive, negative or zero (and is ≥0 or ≤0);
//   * an integer is even or odd;
//   * a finite complex number is algebraic or transcendental.
// A mask whose fields are all Unknown is trivially consistent. This is the
// satisfiability oracle the boolean `ask` solver uses to decide entailment by
// refutation (a proposition is entailed when asserting its negation yields an
// inconsistent mask). Only *genuine* contradictions return false, so a "true"
// here never over-claims.
[[nodiscard]] SYMPP_EXPORT bool assumptions_consistent(AssumptionMask m) noexcept;

}  // namespace sympp
