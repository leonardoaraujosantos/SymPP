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

}  // namespace sympp
