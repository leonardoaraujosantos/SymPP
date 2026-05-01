#pragma once

// Extrema and monotonicity helpers built on diff + solve.
//
// stationary_points(f, var) — roots of f'(var) = 0.
// minimum / maximum         — value of f at the critical point that
//                             minimizes / maximizes (numeric comparison).
// is_increasing / is_decreasing — query sign of f' globally; returns
//                             Unknown (nullopt) when undecidable.
//
// Reference: sympy/calculus/util.py — minimum, maximum, is_increasing,
//            stationary_points, etc.

#include <optional>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT std::vector<Expr> stationary_points(
    const Expr& f, const Expr& var);

[[nodiscard]] SYMPP_EXPORT Expr minimum(const Expr& f, const Expr& var);
[[nodiscard]] SYMPP_EXPORT Expr maximum(const Expr& f, const Expr& var);

[[nodiscard]] SYMPP_EXPORT std::optional<bool> is_increasing(
    const Expr& f, const Expr& var);
[[nodiscard]] SYMPP_EXPORT std::optional<bool> is_decreasing(
    const Expr& f, const Expr& var);

}  // namespace sympp
