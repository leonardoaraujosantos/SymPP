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

}  // namespace sympp::detail
