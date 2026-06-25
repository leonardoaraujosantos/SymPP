#pragma once

// Structured code-generation AST and emitters.
//
// Complements the string-only c_function / cxx_function helpers in
// printing.hpp by exposing a small, inspectable intermediate representation:
//
//   Assignment          — a single `lhs = rhs` statement (lhs is a Symbol)
//   CodeBlock           — a sequence of assignments plus a returned expression
//   FunctionDefinition  — a named function over a list of argument Symbols
//
// cse_codeblock(expr) builds a CodeBlock by common-subexpression elimination
// (reusing simplify::cse), turning each shared subterm into a temporary
// assignment `t0`, `t1`, ... and leaving the rewritten expression as the
// returned value.
//
// emit_c / emit_cxx render a complete C / C++ function, printing each
// temporary as `const double tN = <ccode/cxxcode>;` followed by `return ...;`.
//
// Reference: sympy/utilities/codegen.py, sympy/simplify/cse_main.py

#include <string>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp::codegen {

// A single `lhs = rhs;` statement. `lhs` is expected to be a Symbol Expr.
struct Assignment {
    Expr lhs;
    Expr rhs;
};

// A straight-line block: ordered assignments followed by a returned value.
struct CodeBlock {
    std::vector<Assignment> body;
    Expr ret;
};

// A named function: argument Symbols plus the body block.
struct FunctionDefinition {
    std::string name;
    std::vector<Expr> args;
    CodeBlock block;
};

// Build a CodeBlock from `expr` via common-subexpression elimination. Each
// shared subterm becomes a temporary assignment whose lhs is a Symbol named
// `t0`, `t1`, ...; `ret` holds the rewritten expression in terms of those
// temporaries.
[[nodiscard]] SYMPP_EXPORT CodeBlock cse_codeblock(const Expr& expr);

// Render a complete C function. Temporaries are emitted as
// `const double tN = <ccode(rhs)>;` and the block's return as
// `return <ccode(ret)>;`.
[[nodiscard]] SYMPP_EXPORT std::string emit_c(const FunctionDefinition& fn);

// Render a complete C++ function, using cxxcode for each rhs.
[[nodiscard]] SYMPP_EXPORT std::string emit_cxx(const FunctionDefinition& fn);

}  // namespace sympp::codegen
