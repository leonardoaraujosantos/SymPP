#pragma once

// Mellin–Barnes definite integration via the Meijer-G function — the core of
// SymPy's `meijerint`. The Mellin transform of a Meijer-G is a ratio of Gamma
// functions, which makes ∫₀^∞ tractable in closed form:
//
//   M[G^{m,n}_{p,q}(η·x)](s) = ∫₀^∞ x^{s−1} G(η·x) dx
//     = η^{−s} · Πⱼ₌₁ᵐ Γ(b_j+s) · Πⱼ₌₁ⁿ Γ(1−a_j−s)
//              / ( Πⱼ₌ₘ₊₁^q Γ(1−b_j−s) · Πⱼ₌ₙ₊₁^p Γ(a_j+s) )
//
// and ∫₀^∞ G(η·x) dx is that transform evaluated at s = 1.
//
// Reference: sympy/integrals/meijerint.py, DLMF 16.17 / Meijer-G Mellin transform.

#include <optional>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// Mellin transform M[G(η·var)](s) of a Meijer-G whose argument is η·var (η free
// of var). Returns the Gamma-ratio formula above, or std::nullopt when `G` is
// not a Meijer-G in `var`.
[[nodiscard]] SYMPP_EXPORT std::optional<Expr> meijerg_mellin_transform(const Expr& G,
                                                                        const Expr& var,
                                                                        const Expr& s);

// Definite integral ∫₀^∞ G(η·var) dvar — the Mellin transform at s = 1. Returns
// std::nullopt when `G` is not a Meijer-G in `var` or the integral does not
// converge to a finite value (a Gamma pole in the formula).
[[nodiscard]] SYMPP_EXPORT std::optional<Expr> meijerg_integrate_0_inf(const Expr& G,
                                                                       const Expr& var);

}  // namespace sympp
