#pragma once

// simplify, together, collect, cancel.
//
// Phase 5 minimal: lightweight dispatcher built on top of the canonical
// form (Add term collection, Mul base collection from Phase 1f.2) plus
// expand from Phase 1i. The full SymPy simplify orchestrator chains a
// dozen sub-simplifiers (trigsimp, powsimp, radsimp, sqrtdenest, fu, etc.);
// we ship the high-leverage subset here.
//
// Reference: sympy/simplify/simplify.py::simplify

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// Canonical simplification: re-runs the Add / Mul / Pow factories so any
// substitution-driven simplifications collapse, then applies expand to
// flush nested distributions, then re-runs through the canonical form.
[[nodiscard]] SYMPP_EXPORT Expr simplify(const Expr& e);

// collect(expr, var) — group like-powers of `var` and combine their
// coefficients. Builds a polynomial in `var` and reconstitutes.
[[nodiscard]] SYMPP_EXPORT Expr collect(const Expr& e, const Expr& var);

// together(expr) is provided by sympp/polys/poly.hpp — it lives in the
// polynomial subsystem because it walks numerator/denominator structure.

// powsimp(expr) — combine same-exponent powers within Mul:
//   x^a * y^a → (x*y)^a   when assumptions justify (positive bases or
//   integer exponent).
//
// The complementary same-base combining (x^a * x^b → x^(a+b)) is already
// handled by Mul base collection in the canonical form.
//
// Reference: sympy/simplify/powsimp.py::powsimp
[[nodiscard]] SYMPP_EXPORT Expr powsimp(const Expr& e);

// expand_power_exp(expr) — apply x^(a+b) → x^a * x^b at every Pow whose
// exponent is an Add. Targeted variant of expand().
//
// Reference: sympy/core/expr.py expand(power_exp=True)
[[nodiscard]] SYMPP_EXPORT Expr expand_power_exp(const Expr& e);

// expand_power_base(expr) — apply (x*y)^n → x^n * y^n at every Pow whose
// base is a Mul, when the exponent is an integer (always safe) or all
// factors of the base are positive.
//
// Reference: sympy/core/expr.py expand(power_base=True)
[[nodiscard]] SYMPP_EXPORT Expr expand_power_base(const Expr& e);

// radsimp(expr) — rationalize denominators containing radicals. Handles
// the common case of a binomial denominator a + b*sqrt(c) by multiplying
// numerator and denominator by the conjugate a - b*sqrt(c). Larger
// denominators (more than one distinct radical) pass through unchanged.
//
// Reference: sympy/simplify/radsimp.py::radsimp
[[nodiscard]] SYMPP_EXPORT Expr radsimp(const Expr& e);

// sqrtdenest(expr) — denest sqrt(a + b*sqrt(c)) → sqrt(d) + sqrt(e) when
// the expression admits a rational denesting (i.e. when c² - 4*(a²-b²)/c
// — the discriminant of the auxiliary quadratic — is a perfect square).
//
// Reference: sympy/simplify/sqrtdenest.py::sqrtdenest
[[nodiscard]] SYMPP_EXPORT Expr sqrtdenest(const Expr& e);

// trigsimp(expr) — apply trigonometric identities (Pythagorean, ratio,
// reciprocal). Minimal Fu-style subset:
//   * sin²(x) + cos²(x) → 1   (with shared coefficient: a*sin²(x) + a*cos²(x) → a)
//   * tan(x) → sin(x)/cos(x), cot(x) → cos(x)/sin(x)
//   * sec(x) → 1/cos(x), csc(x) → 1/sin(x)
//   * 1 - sin²(x) → cos²(x), 1 - cos²(x) → sin²(x)
// Walks the tree at every Add and applies the rules to a fixed point.
//
// Reference: sympy/simplify/trigsimp.py, sympy/simplify/fu.py
[[nodiscard]] SYMPP_EXPORT Expr trigsimp(const Expr& e);

}  // namespace sympp
