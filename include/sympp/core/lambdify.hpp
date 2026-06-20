#pragma once

// lambdify — compile an expression into a fast native double-precision callable.
//
// SymPy's `lambdify` turns an expression into a Python function; SymPP's emits a
// `std::function<double(const std::vector<double>&)>` instead. The expression is
// walked once into a closure tree, so the returned callable can be invoked
// repeatedly with no further symbolic work — ideal for plotting, root-finding,
// quadrature, or any tight numeric loop.
//
//   auto x = symbol("x"), y = symbol("y");
//   auto f = lambdify({x, y}, x * x + sin(y));
//   f({3.0, 0.0});          // → 9.0
//
// Real-valued only. lambdify throws std::runtime_error at compile time if the
// body contains the imaginary unit, an unsupported (e.g. Bessel/Zeta) function,
// or a free symbol not listed in `vars`.
//
// Reference: sympy/utilities/lambdify.py (this is the interpreter backend, not
// the LLVM-JIT one).

#include <functional>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

using CompiledExpr = std::function<double(const std::vector<double>&)>;

[[nodiscard]] SYMPP_EXPORT CompiledExpr lambdify(const std::vector<Expr>& vars,
                                                 const Expr& body);

// Single-variable convenience.
[[nodiscard]] SYMPP_EXPORT CompiledExpr lambdify(const Expr& var, const Expr& body);

}  // namespace sympp
