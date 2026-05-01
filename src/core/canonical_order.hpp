#pragma once

// Internal helper — stable canonical ordering for Add/Mul args.
// Header-only so both add.cpp and mul.cpp can inline it.

#include <cstddef>

#include <sympp/core/basic.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp::detail {

// Returns true if `a` should sort before `b` in canonical Add/Mul args.
// Order: by type_id, then by hash, then by str representation. Deterministic
// and stable across runs. Numbers handled separately by callers.
[[nodiscard]] inline bool canonical_less(const Expr& a, const Expr& b) noexcept {
    auto ta = static_cast<std::size_t>(a->type_id());
    auto tb = static_cast<std::size_t>(b->type_id());
    if (ta != tb) return ta < tb;
    auto ha = a->hash();
    auto hb = b->hash();
    if (ha != hb) return ha < hb;
    // Hash collision — fall back to lexicographic string comparison.
    return a->str() < b->str();
}

}  // namespace sympp::detail
