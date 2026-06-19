#pragma once

// Internal helpers used by Add/Mul/Pow assumption rules. Header-only.

#include <cstddef>
#include <optional>

#include <sympp/core/assumption_key.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/fwd.hpp>

namespace sympp::detail {

template <typename Range>
[[nodiscard]] inline bool all_args_have(const Range& args,
                                        AssumptionKey k,
                                        bool target) noexcept {
    for (const auto& a : args) {
        auto v = direct_ask(a, k);
        if (!v.has_value() || *v != target) return false;
    }
    return true;
}

template <typename Range>
[[nodiscard]] inline bool any_arg_has(const Range& args,
                                      AssumptionKey k,
                                      bool target) noexcept {
    for (const auto& a : args) {
        auto v = direct_ask(a, k);
        if (v.has_value() && *v == target) return true;
    }
    return false;
}

// For a SUM: true iff every term is algebraic except exactly one transcendental
// term. Such a sum is transcendental (and not algebraic): if it were algebraic,
// subtracting the algebraic terms would make the transcendental one algebraic —
// a contradiction, since the algebraic numbers are closed under subtraction.
template <typename Range>
[[nodiscard]] inline bool sum_forces_transcendental(const Range& args) noexcept {
    int transcendental = 0;
    for (const auto& a : args) {
        if (direct_ask(a, AssumptionKey::Transcendental) == true) { ++transcendental; continue; }
        if (direct_ask(a, AssumptionKey::Algebraic) != true) return false;
    }
    return transcendental == 1;
}

// For a PRODUCT: true iff exactly one factor is transcendental and every other
// factor is a nonzero algebraic. The nonzero requirement matters because a zero
// algebraic factor would collapse the product to 0 (which is algebraic).
template <typename Range>
[[nodiscard]] inline bool product_forces_transcendental(const Range& args) noexcept {
    int transcendental = 0;
    for (const auto& a : args) {
        if (direct_ask(a, AssumptionKey::Transcendental) == true) { ++transcendental; continue; }
        if (direct_ask(a, AssumptionKey::Algebraic) != true) return false;
        if (direct_ask(a, AssumptionKey::Nonzero) != true) return false;
    }
    return transcendental == 1;
}

// Provably irrational at the node level: a directly-irrational value, or a real
// transcendental (every real transcendental is irrational, since the rationals
// are algebraic). The node-level `ask(Irrational)` is not itself derived from
// transcendental ∧ real, so spell that path out here.
[[nodiscard]] inline bool node_is_irrational(const Expr& a) noexcept {
    if (direct_ask(a, AssumptionKey::Irrational) == true) return true;
    return direct_ask(a, AssumptionKey::Transcendental) == true
           && direct_ask(a, AssumptionKey::Real) == true;
}

// For a SUM: true iff exactly one term is irrational and every other term is
// rational. The rationals are closed under subtraction, so subtracting the
// rational terms would make the irrational one rational — a contradiction.
// (rational + irrational = irrational; a sum of two irrationals may be rational,
// e.g. √2 + (−√2) = 0, so more than one irrational term leaves it Unknown.)
template <typename Range>
[[nodiscard]] inline bool sum_forces_irrational(const Range& args) noexcept {
    int irrational = 0;
    for (const auto& a : args) {
        if (node_is_irrational(a)) { ++irrational; continue; }
        if (direct_ask(a, AssumptionKey::Rational) != true) return false;
    }
    return irrational == 1;
}

// For a PRODUCT: true iff exactly one factor is irrational and every other factor
// is a nonzero rational. A nonzero rational times an irrational is irrational; a
// zero factor would collapse the product to 0 (rational), and two irrational
// factors may multiply to a rational (√2·√2 = 2), so both are excluded.
template <typename Range>
[[nodiscard]] inline bool product_forces_irrational(const Range& args) noexcept {
    int irrational = 0;
    for (const auto& a : args) {
        if (node_is_irrational(a)) { ++irrational; continue; }
        if (direct_ask(a, AssumptionKey::Rational) != true) return false;
        if (direct_ask(a, AssumptionKey::Nonzero) != true) return false;
    }
    return irrational == 1;
}

}  // namespace sympp::detail
