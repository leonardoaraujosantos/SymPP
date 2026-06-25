#pragma once

// autowrap — compile an expression to a native callable via the system C compiler.
//
// SymPy's `autowrap` generates source code for an expression, compiles it with a
// real toolchain, and binds the resulting binary back as a callable. SymPP's
// analogue emits C (reusing printing::c_function), invokes the system C compiler
// (`${CC:-cc} -O2 -shared -fPIC`) to produce a shared object, then dlopen/dlsym
// the function and wraps it as a CompiledExpr.
//
//   auto x = symbol("x");
//   if (autowrap_available()) {
//       auto f = autowrap(x, x * x + integer(1));   // x² + 1
//       f({3.0});                                    // → 10.0
//   }
//
// This differs from lambdify (closure interpreter) and lambdify_jit (LLVM) by
// going through generated C source and the platform compiler/linker — the same
// mechanism as SymPy's autowrap C backend.
//
// Real-valued only, with the same body restrictions as c_function/ccode. The
// returned callable keeps the dlopen handle alive for its whole lifetime; temp
// .c/.so files are removed once loaded. Throws std::runtime_error if no C
// compiler is available — callers should gate on autowrap_available().
//
// Reference: sympy/utilities/autowrap.py

#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/lambdify.hpp>  // CompiledExpr
#include <sympp/fwd.hpp>

namespace sympp::codegen {

// True iff a usable C compiler (from $CC, else cc/gcc/clang on PATH) is found.
[[nodiscard]] SYMPP_EXPORT bool autowrap_available() noexcept;

// Compile `body` (a function of `vars`) to native code and return a callable.
// Throws std::runtime_error if no C compiler is available or compilation fails.
[[nodiscard]] SYMPP_EXPORT CompiledExpr autowrap(const std::vector<Expr>& vars,
                                                 const Expr& body);

// Single-variable convenience.
[[nodiscard]] SYMPP_EXPORT CompiledExpr autowrap(const Expr& var,
                                                 const Expr& body);

}  // namespace sympp::codegen
