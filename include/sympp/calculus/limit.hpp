#pragma once

// limit(expr, var, target) — limit of expr as var -> target.
//
// Phase 6 minimal: substitute and simplify. Doesn't handle indeterminate
// forms (0/0, ∞·0, etc.) — Gruntz-quality limit reasoning ships in a
// follow-up phase.
//
// Reference: sympy/series/limits.py and sympy/series/gruntz.py

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr limit(const Expr& expr, const Expr& var,
                                        const Expr& target);

}  // namespace sympp
