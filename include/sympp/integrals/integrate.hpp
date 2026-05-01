#pragma once

// integrate(expr, var) — antiderivative.
// integrate(expr, var, lower, upper) — definite via Newton-Leibniz.
//
// Strategy chain (each falls through to the next on failure):
//
//   1. Linearity over Add — split sums into per-term integrals.
//   2. Elementary table on each term:
//        * Constant w.r.t. var          → c·var
//        * (a·var + b)^n                → (a·var + b)^(n+1) / (a·(n+1))
//          n == -1                      → log(a·var + b) / a
//        * sin(a·var + b)               → -cos(a·var + b) / a
//        * cos(a·var + b)               → sin(a·var + b) / a
//        * exp(a·var + b)               → exp(a·var + b) / a
//      Mul case: pull constant factors out, recurse via the public
//      integrate() so an Add inner gets linearity.
//   3. Trig reductions via identity rewrites:
//        sin²(u) → (1 - cos(2u))/2,  cos²(u) → (1 + cos(2u))/2
//        sin(p)·cos(q) → (sin(p+q) + sin(p-q))/2  (product-to-sum)
//   4. Rational-function path: together() into single fraction,
//      polynomial-divide for the improper part, apart() decomposition
//      of the proper part, recurse to integrate each piece.
//   5. Integration by parts:
//        * standalone log(affine)
//        * Mul of polynomial × {sin, cos, exp}(affine) — recursive
//   6. Numeric quadrature (vpaintegral) is available separately for
//      definite numeric values when symbolic methods don't succeed.
//
// Anything outside the chain returns the input wrapped in a placeholder
// "Integral" UndefinedFunction so downstream code can recognize the
// unevaluated case. Full Risch / Meijer G remains deferred.
//
// Reference: sympy/integrals/integrals.py and manualintegrate.py

#include <optional>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr integrate(const Expr& expr, const Expr& var);
[[nodiscard]] SYMPP_EXPORT Expr integrate(const Expr& expr, const Expr& var,
                                            const Expr& lower, const Expr& upper);

// manualintegrate — same pipeline as integrate(), but returns nullopt
// when the symbolic strategies couldn't find a closed form (instead of
// the Integral(_, _) marker that integrate() emits). Useful as a
// success/failure dispatch for callers that want to fall back to
// numeric quadrature.
[[nodiscard]] SYMPP_EXPORT std::optional<Expr> manualintegrate(
    const Expr& expr, const Expr& var);

}  // namespace sympp
