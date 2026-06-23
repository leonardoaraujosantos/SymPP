#pragma once

// product — closed-form symbolic products ∏_{var=lo}^{hi} expr.
//
// product(expr, var, lo, hi) computes the product when a closed form is
// recognized. Coverage:
//   * Single-term range (hi == lo): ∏ = expr(lo)
//   * Constant w.r.t. var: c → c^(hi − lo + 1)
//   * Multiplicative split: ∏ (f·g) = (∏ f)·(∏ g); constant factors → c^N
//   * Powers: ∏ base(k)^c = (∏ base(k))^c (constant c);
//             ∏ r^(c·k + d) = r^(c·Σk + d·N) for a var-free base r (geometric)
//   * Polynomial in var with closed-form roots: factor into linear pieces and
//     use the Gamma ratio ∏_{k=lo}^{hi} (k − r) = Γ(hi − r + 1)/Γ(lo − r).
//     Gives ∏ k = n! (factorial), ∏ (k+a) = Γ(n+a+1)/Γ(a+1), ∏ (k²−1), …
//
// Returns an unevaluated `Product(expr, var, lo, hi)` marker when no closed form
// is recognized, mirroring summation's `Sum(...)` — a product is never silently
// dropped to its bare term.
//
// Reference: sympy/concrete/products.py::Product / .doit

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr product(const Expr& expr, const Expr& var,
                                        const Expr& lo, const Expr& hi);

}  // namespace sympp
