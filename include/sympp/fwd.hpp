#pragma once

// Forward declarations for SymPP public types.
// Include this header from your own headers to avoid pulling in implementations.

#include <memory>

namespace sympp {

class Basic;
class Symbol;

// Public alias: every user-facing expression is shared, immutable, polymorphic.
using Expr = std::shared_ptr<const Basic>;

}  // namespace sympp
