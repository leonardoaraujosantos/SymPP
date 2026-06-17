#pragma once

// Internal helpers used by Add/Mul/Pow assumption rules. Header-only.

#include <cstddef>
#include <optional>

#include <sympp/core/assumption_key.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/fwd.hpp>

namespace sympp::detail {

template <typename Range>
[[nodiscard]] inline bool all_args_have(const Range& args,
                                        AssumptionKey k,
                                        bool target) noexcept {
    for (const auto& a : args) {
        auto v = a->ask(k);
        if (!v.has_value() || *v != target) return false;
    }
    return true;
}

template <typename Range>
[[nodiscard]] inline bool any_arg_has(const Range& args,
                                      AssumptionKey k,
                                      bool target) noexcept {
    for (const auto& a : args) {
        auto v = a->ask(k);
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
        if (a->ask(AssumptionKey::Transcendental) == true) { ++transcendental; continue; }
        if (a->ask(AssumptionKey::Algebraic) != true) return false;
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
        if (a->ask(AssumptionKey::Transcendental) == true) { ++transcendental; continue; }
        if (a->ask(AssumptionKey::Algebraic) != true) return false;
        if (a->ask(AssumptionKey::Nonzero) != true) return false;
    }
    return transcendental == 1;
}

}  // namespace sympp::detail
