#pragma once

// Hash and equality functors over Expr, plus convenience map/set aliases.
// Use these whenever you need an unordered container keyed on Expr — they
// guarantee structural identity (two different shared_ptr but structurally
// equal Exprs map to the same bucket).

#include <cstddef>
#include <unordered_map>
#include <unordered_set>

#include <sympp/core/api.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

struct ExprHash {
    [[nodiscard]] std::size_t operator()(const Expr& e) const noexcept {
        return e ? e->hash() : 0;
    }
};

struct ExprEq {
    [[nodiscard]] bool operator()(const Expr& a, const Expr& b) const noexcept {
        return a == b;  // structural equality via Basic::equals
    }
};

template <typename V>
using ExprMap = std::unordered_map<Expr, V, ExprHash, ExprEq>;

using ExprSet = std::unordered_set<Expr, ExprHash, ExprEq>;

}  // namespace sympp
