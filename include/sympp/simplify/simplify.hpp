#pragma once

// Simplification orchestrator + per-strategy entry points.
//
// `simplify(e)` chains the following sub-simplifiers in order
// (each is a no-op when its pattern doesn't match):
//
//   1. re-canonicalize (Add term collection, Mul base collection)
//   2. expand          (distribute Mul over Add, expand integer Pow)
//   3. trigsimp        (sin² + cos² → 1, etc.)
//   4. powsimp         (combine same-exponent powers)
//   5. combsimp        (factorial(n+k)/factorial(n) → falling product)
//   6. gammasimp       (gamma(n+1)/gamma(n) → n)
//   7. radsimp         (rationalize binomial radical denominators)
//   8. sqrtdenest      (Borodin denesting on perfect-square discriminants)
//   9. re-canonicalize (final pass)
//
// Each step is also exposed as a free function so callers can invoke
// just the rule they want without paying for the others.
//
// Deferred-deep: hyperexpand (hypergeometric expansion) — this is the
// blocker for Meijer-G-driven integration as well; treated as one
// research-level subsystem to land jointly. See docs/04-roadmap.md.
//
// Reference: sympy/simplify/simplify.py::simplify and the supporting
// per-rule files in sympy/simplify/.

#include <utility>
#include <vector>

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

// combsimp(expr) — simplify combinatorial expressions. Minimal coverage:
//   * factorial(n + k) / factorial(n) → (n+k)*(n+k-1)*...*(n+1) for small
//     non-negative integer k (the "falling product").
//   * gamma(n + 1) / gamma(n) → n (analogous, via gamma's recurrence).
//
// Reference: sympy/simplify/combsimp.py
[[nodiscard]] SYMPP_EXPORT Expr combsimp(const Expr& e);

// gammasimp(expr) — apply gamma function identities. Minimal subset:
//   * gamma(n + 1) → n * gamma(n) for non-negative integer offsets
//   * gamma(1) → 1
// Combined with combsimp covers the common factorial/gamma manipulations.
//
// Reference: sympy/simplify/gammasimp.py
[[nodiscard]] SYMPP_EXPORT Expr gammasimp(const Expr& e);

// cse(expr) — common subexpression elimination. Walks the expression,
// finds subtrees that appear at least twice, replaces each with a fresh
// temporary symbol, and returns the substitution list alongside the
// simplified expression. Useful for code generation and any context
// that benefits from reifying shared subterms as variables.
//
// Reference: sympy/simplify/cse_main.py::cse
struct CSEResult {
    std::vector<std::pair<Expr, Expr>> substitutions;  // (temp, definition)
    Expr expr;
};
[[nodiscard]] SYMPP_EXPORT CSEResult cse(const Expr& e);

// nsimplify(expr) — heuristic numeric→exact-symbolic guesser. Try to
// recognize a Float / Rational as an exact symbolic constant from a
// short table: rational p/q (small denominators), pi, e, sqrt of small
// rationals, log of small integers. Returns the input unchanged when no
// clean exact form is found within tolerance.
//
// Reference: sympy/simplify/simplify.py::nsimplify
[[nodiscard]] SYMPP_EXPORT Expr nsimplify(const Expr& e, int dps = 15);

// trigsimp(expr) — Pythagorean-identity reduction over Add.
//
// At each Add node, group sin²/cos² terms by their argument and
// rewrite as
//
//     a·sin²(x) + b·cos²(x)  →  b + (a − b)·sin²(x)
//
// which collapses to `b` when a == b (the canonical
// sin² + cos² = 1 case). Walks the tree recursively so embedded
// occurrences inside Mul / Pow / Function arguments get reduced too.
//
// Not yet shipped: full Fu rule table (double-angle, half-angle,
// sum-to-product, cot/sec/csc reciprocal rewrites). The single
// Pythagorean rule covers the most common simplification need; the
// rest is deferred-deep.
//
// Reference: sympy/simplify/trigsimp.py, sympy/simplify/fu.py
[[nodiscard]] SYMPP_EXPORT Expr trigsimp(const Expr& e);

}  // namespace sympp
