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
//   ₀F₁(; 1/2; z)       = cosh(2√z)
//   ₀F₁(; 3/2; z)       = sinh(2√z) / (2√z)
//   ₁F₁(1; 2; z)        = (e^z − 1) / z
//   ₁F₁(1; 3/2; z)      = √π·e^z·erf(√z) / (2√z)
//   ₂F₁(1, 1; 2; z)     = −log(1 − z) / z
//   ₂F₁(1/2, 1; 3/2; z)   = atanh(√z) / √z   (→ arctan(y)/y at z = −y²)
//   ₂F₁(1/2, 1/2; 3/2; z) = asin(√z) / √z    (→ arcsinh(y)/y at z = −y²)
//   ₂F₁(a, b; b; z)     = (1 − z)^(−a)                   ← via cancellation in factory
//
// Walks the tree recursively so embedded `hyper(...)` nodes inside
// Mul/Add/Pow/Function are rewritten too. Anything not matching the
// table is returned unchanged.
//
// Meijer-G expansion uses Slater's theorem in the generic lower-parameter
// case (no two of b₁…b_m differ by an integer): the G-function is rewritten
// as Σ_k A_k·z^{b_k}·pF_{q−1}(…) and each hypergeometric is then expanded by
// the table above — so e.g. G^{1,1}_{1,1} → 1/(z+1), G^{2,0}_{0,2} →
// √π(cosh−sinh)(2√z). The confluent (integer-spaced) case and Mellin–Barnes
// definite integration remain staged work (OpenSpec: add-meijerg-slater-engine).
//
// Reference: sympy/simplify/hyperexpand.py — `hyperexpand`,
// `_check_hyper`, the lookup tables `formulae` and `meijerg_formulae`.

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr hyperexpand(const Expr& e);

}  // namespace sympp
