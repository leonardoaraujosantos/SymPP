#pragma once

// series(expr, var, x0, n) — Taylor expansion of `expr` about `x0` up to
// degree `n - 1`. Computed via the standard Taylor formula:
//
//     sum_{k=0}^{n-1} (1/k!) * d^k expr / dx^k |_{x = x0} * (x - x0)^k
//
// Uses our diff() engine, so closed form is available for everything diff
// covers. Phase 6 minimal: returned as a plain Expr (no O(x^n) marker).
//
// Reference: sympy/series/series.py / sympy/core/expr.py::series

#include <cstddef>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr series(const Expr& expr, const Expr& var,
                                         const Expr& x0, std::size_t n);

}  // namespace sympp
