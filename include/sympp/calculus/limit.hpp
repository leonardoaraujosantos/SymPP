#pragma once

// limit(expr, var, target) — limit of expr as var -> target.
//
// Strategy:
//   1. Direct substitution; if the value is finite and determinate,
//      return it.
//   2. Detect 0/0 indeterminate form by inspecting numerator and
//      denominator (after `together`) separately at the target.
//   3. Apply L'Hôpital iteratively up to 10 times: limit(f/g) =
//      limit(f'/g') when both → 0.
//
// Coverage: 0/0 forms (sin(x)/x, (1-cos(x))/x², log(1+x)/x, …) and
// the directly-substitutable cases.
//
// Not yet covered (deep-deferred): full Gruntz at infinity (∞/∞,
// 0·∞, ∞-∞), one-sided limits, limits requiring asymptotic
// expansion. These need an `Infinity` singleton and the
// asymptotic-comparison machinery from sympy/series/gruntz.py.
//
// Reference: sympy/series/limits.py and sympy/series/gruntz.py

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr limit(const Expr& expr, const Expr& var,
                                        const Expr& target);

// One-sided limit. `dir` selects the approach side: +1 from the right
// (var → target⁺), −1 from the left (var → target⁻), 0 the two-sided limit
// (identical to the 3-argument overload). One-sided queries resolve the cases a
// two-sided limit leaves undefined: a finite pole (1/x → +∞ as x→0⁺, −∞ as
// x→0⁻) and a sign/abs discontinuity (|x|/x → ±1). Note the 3-argument overload
// stays two-sided; SymPy's `limit` instead defaults to dir='+'.
[[nodiscard]] SYMPP_EXPORT Expr limit(const Expr& expr, const Expr& var,
                                        const Expr& target, int dir);

}  // namespace sympp
