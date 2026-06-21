#pragma once

// LLVM ORC-JIT backend for lambdify — compiles an expression all the way to
// native machine code for maximum numeric throughput (the interpreter backend
// in <core/lambdify.hpp> is the portable default).
//
//   auto f = lambdify_jit({x}, sin(x) * x);   // native code, same signature
//   f({1.5});
//
// Available only when SymPP was built with LLVM (CMake option
// SYMPP_ENABLE_LLVM_JIT, auto-on when LLVM is found). When unavailable,
// lambdify_jit throws std::runtime_error; check lambdify_jit_available() first.
//
// Reference: sympy/utilities/lambdify.py — the JIT-compiled tier.

#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/lambdify.hpp>  // CompiledExpr
#include <sympp/fwd.hpp>

namespace sympp {

// True iff the LLVM ORC-JIT backend was compiled in.
[[nodiscard]] SYMPP_EXPORT bool lambdify_jit_available() noexcept;

// JIT-compile `body` over the ordered variable list into a native callable with
// the same contract as lambdify(). Throws std::runtime_error if the JIT backend
// is unavailable or the body uses an unsupported construct. The returned
// std::function owns the compiled code (kept alive for its lifetime).
[[nodiscard]] SYMPP_EXPORT CompiledExpr lambdify_jit(const std::vector<Expr>& vars,
                                                     const Expr& body);
[[nodiscard]] SYMPP_EXPORT CompiledExpr lambdify_jit(const Expr& var, const Expr& body);

}  // namespace sympp
