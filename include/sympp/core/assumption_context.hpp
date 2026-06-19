#pragma once

// assumption_context — a scoped (`assuming`) set of facts consulted by ask().
//
// Mirrors SymPy's `with assuming(Q.positive(x)): ...`. An `assuming` RAII scope
// augments the global assumption context with facts about a symbol (or any
// expression key) for the duration of its lifetime; every query routed through
// ask() — and therefore refine(), the limit/integration engines, etc. — sees
// those facts. Scopes nest with stack discipline.
//
//     {
//         assuming pos(x, AssumptionMask{}.set_positive(true));
//         refine(sqrt(pow(x, integer(2))));   // → x
//     }   // fact retracted here
//
// Reference: sympy/assumptions/assume.py — `assuming` / global_assumptions.

#include <optional>

#include <sympp/core/api.hpp>
#include <sympp/core/assumption_key.hpp>
#include <sympp/core/assumption_mask.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// First asserted value for predicate `k` about `e` across the active scopes
// (innermost first), or nullopt when no scope speaks to it.
[[nodiscard]] SYMPP_EXPORT std::optional<bool> assumption_context_get(
    const Expr& e, AssumptionKey k) noexcept;

// True when any `assuming` scope is currently active (lets ask() skip the
// lookup entirely in the common no-context case).
[[nodiscard]] SYMPP_EXPORT bool assumption_context_active() noexcept;

class SYMPP_EXPORT assuming {
public:
    assuming(const Expr& target, AssumptionMask facts);
    ~assuming();
    assuming(const assuming&) = delete;
    assuming& operator=(const assuming&) = delete;
    assuming(assuming&&) = delete;
    assuming& operator=(assuming&&) = delete;
};

}  // namespace sympp
