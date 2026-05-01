#pragma once

// Asymptotes and inflection points.
//
//   * vertical_asymptotes(f, var) — roots of the simplified denominator
//     of f, i.e. points where f blows up. cancel() removes common
//     factors with the numerator first so removable singularities don't
//     show up.
//
//   * inflection_points(f, var) — roots of f''(var). This is a
//     superset of true inflection points (a sign change of f'' is also
//     required), but matches SymPy's coarse helper and is the textbook
//     starting point for inflection analysis.
//
// horizontal_asymptotes is intentionally absent for now: it requires
// limit-at-infinity machinery which depends on a proper Infinity
// singleton — that lands when we do Phase 6's full Gruntz upgrade.
//
// Reference: sympy/calculus/util.py — Polynomial inflection helpers,
//            sympy/calculus/singularities.py — singularities()

#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT std::vector<Expr> vertical_asymptotes(
    const Expr& f, const Expr& var);

[[nodiscard]] SYMPP_EXPORT std::vector<Expr> inflection_points(
    const Expr& f, const Expr& var);

}  // namespace sympp
