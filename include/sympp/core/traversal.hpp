#pragma once

// Tree traversal — substitution, search, and symbol enumeration.
//
// References:
//   sympy/core/basic.py::Basic.subs / xreplace / has / free_symbols
//   MATLAB Symbolic Math Toolbox: subs() — see User's Guide §3.

#include <sympp/core/api.hpp>
#include <sympp/core/expr_collections.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// Replace every structurally-equal occurrence of `old_value` in `expr` with
// `new_value`, then re-canonicalize the result. For multiple substitutions,
// either chain calls or use `xreplace`.
//
// Reference: sympy/core/basic.py::Basic.subs
[[nodiscard]] SYMPP_EXPORT Expr subs(const Expr& expr,
                                     const Expr& old_value,
                                     const Expr& new_value);

// Apply a *simultaneous* substitution map to `expr`. Each subexpression that
// matches a key in `mapping` is replaced by its value; recursion does not
// re-process replaced subtrees, so the substitution is atomic. The result is
// re-canonicalized so e.g. xreplace(x+x, {x: 1}) returns 2 (not Add(1, 1)).
//
// Reference: sympy/core/basic.py::Basic.xreplace
[[nodiscard]] SYMPP_EXPORT Expr xreplace(const Expr& expr, const ExprMap<Expr>& mapping);

// True iff `target` appears anywhere in `expr` as a structurally-equal subtree.
//
// Reference: sympy/core/basic.py::Basic.has
[[nodiscard]] SYMPP_EXPORT bool has(const Expr& expr, const Expr& target);

// All Symbol leaves reachable from `expr`. Returns a structural-identity set,
// so `Symbol("x")` always coalesces.
//
// Reference: sympy/core/basic.py::Basic.free_symbols
[[nodiscard]] SYMPP_EXPORT ExprSet free_symbols(const Expr& expr);

}  // namespace sympp
