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
// Order: by type_id, then by lexicographic str(), then by hash. The string
// comparison is the primary tiebreaker (rather than hash) so the canonical
// printed form is platform-independent — `x + y` reads `x + y` rather than
// whatever order the standard library's hash-of-string happens to pick on
// the host stdlib. Hash is the final tiebreaker for the rare case where
// two structurally distinct nodes share both type_id and str() (e.g.
// Symbol("x") with two different assumption masks).
//
// Reference: sympy/core/basic.py::_args_sort_key — same string-first idea.
[[nodiscard]] inline bool canonical_less(const Expr& a, const Expr& b) noexcept {
    auto ta = static_cast<std::size_t>(a->type_id());
    auto tb = static_cast<std::size_t>(b->type_id());
    if (ta != tb) return ta < tb;
    auto sa = a->str();
    auto sb = b->str();
    if (sa != sb) return sa < sb;
    return a->hash() < b->hash();
}

}  // namespace sympp::detail
