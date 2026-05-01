#pragma once

// expand — distribute multiplications over additions and expand integer
// powers of sums.
//
// Examples:
//     expand(x * (y + z))         -> x*y + x*z
//     expand((x + 1)**2)          -> x**2 + 2*x + 1
//     expand((x + 1) * (x + 2))   -> x**2 + 3*x + 2
//
// Cases NOT expanded by Phase 1i (defer to Phase 5 simplify/expand):
//     * Symbolic exponents:          (x+1)**y stays unevaluated
//     * Negative integer exponents:  (x+1)**(-2) stays unevaluated
//     * exp(a+b), log(a*b) — needs functions module (Phase 3)
//
// Reference: sympy/core/function.py::expand,
//            sympy/core/mul.py::_eval_expand_mul,
//            sympy/core/power.py::_eval_expand_power_base.

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr expand(const Expr& e);

}  // namespace sympp
