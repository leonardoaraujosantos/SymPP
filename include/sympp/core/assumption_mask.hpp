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
    std::optional<bool> finite;

    [[nodiscard]] std::optional<bool> get(AssumptionKey k) const noexcept;
    void set(AssumptionKey k, bool value) noexcept;

    // Builder helpers — return *this so they chain.
    AssumptionMask& set_real(bool v) noexcept { real = v; return *this; }
    AssumptionMask& set_rational(bool v) noexcept { rational = v; return *this; }
    AssumptionMask& set_integer(bool v) noexcept { integer = v; return *this; }
    AssumptionMask& set_positive(bool v) noexcept { positive = v; return *this; }
    AssumptionMask& set_negative(bool v) noexcept { negative = v; return *this; }
    AssumptionMask& set_zero(bool v) noexcept { zero = v; return *this; }
    AssumptionMask& set_finite(bool v) noexcept { finite = v; return *this; }

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

}  // namespace sympp
