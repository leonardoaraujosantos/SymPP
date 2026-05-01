#pragma once

// hyperexpand — rewrite recognized hypergeometric / Meijer-G calls into
// elementary closed forms when one is known.
//
// Recognized identities (all over the formal-series rather than the
// analytically-continued domain — assumes the principal branch is fine
// for the caller):
//
//   ₀F₀(z)              = exp(z)                        ← already at construction
//   ₁F₀(a; ; z)         = (1 − z)^(−a)                  ← already at construction
//   ₀F₁(; 1/2; −z²/4)   = cos(z)
//   ₀F₁(; 3/2; −z²/4)   = sin(z) / z
//   ₁F₁(1; 2; z)        = (e^z − 1) / z
//   ₂F₁(1, 1; 2; −z)    = log(1 + z) / z
//   ₂F₁(1, 1; 2; z)     = −log(1 − z) / z
//   ₂F₁(1/2, 1; 3/2; −z²) = arctan(z) / z
//   ₂F₁(1/2, 1/2; 3/2; z²) = arcsin(z) / z
//   ₂F₁(a, b; b; z)     = (1 − z)^(−a)                   ← via cancellation in factory
//
// Walks the tree recursively so embedded `hyper(...)` nodes inside
// Mul/Add/Pow/Function are rewritten too. Anything not matching the
// table is returned unchanged.
//
// Meijer-G expansion is currently a no-op (full Slater-theorem
// implementation is deferred-deep — it is the largest single piece of
// SymPy's `hyperexpand.py` we haven't ported).
//
// Reference: sympy/simplify/hyperexpand.py — `hyperexpand`,
// `_check_hyper`, the lookup tables `formulae` and `meijerg_formulae`.

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr hyperexpand(const Expr& e);

}  // namespace sympp
