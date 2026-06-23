#pragma once

// rewrite — convert an expression between function families, the analogue of
// SymPy's `expr.rewrite(target)`.
//
//   rewrite(sin(x), RewriteTarget::Exp)       → (exp(I·x) − exp(−I·x)) / (2·I)
//   rewrite(tan(x), RewriteTarget::SinCos)    → sin(x) / cos(x)
//   rewrite(factorial(n), RewriteTarget::Gamma)    → gamma(n + 1)
//   rewrite(gamma(z), RewriteTarget::Factorial)    → factorial(z − 1)
//
// The transform is applied throughout the expression tree (inside sums,
// products, powers and other function arguments) and preserves value, so the
// result is always mathematically equal to the input.
//
// Reference: sympy/core/basic.py::Basic.rewrite.

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

enum class RewriteTarget {
    Exp,        // trig & hyperbolic → exponentials
    SinCos,     // tan/cot/sec/csc → sin & cos
    Gamma,      // factorial / binomial / (rising|falling) factorial → Γ
    Factorial,  // Γ / binomial / (rising|falling) factorial → factorial
};

[[nodiscard]] SYMPP_EXPORT Expr rewrite(const Expr& e, RewriteTarget target);

}  // namespace sympp
