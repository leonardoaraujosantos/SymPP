#pragma once

// Internal alias: detail-namespaced view of the public ExprMap.
// Kept separate so add.cpp / mul.cpp can still write `detail::ExprMap`.

#include <sympp/core/expr_collections.hpp>

namespace sympp::detail {

using sympp::ExprHash;
using sympp::ExprEq;

template <typename V>
using ExprMap = sympp::ExprMap<V>;

}  // namespace sympp::detail
