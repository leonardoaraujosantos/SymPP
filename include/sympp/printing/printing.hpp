#pragma once

// Phase 13 — Code generation / printing.
//
// Five concrete printers:
//   * ccode(expr)      — C99 syntax with <math.h> functions (sin, cos, ...)
//   * cxxcode(expr)    — C++ with std::sin, std::pow, std::numbers::pi_v<double>
//   * fcode(expr)      — Fortran 90+ (** for power, intrinsic SIN/COS/EXP)
//   * latex(expr)      — TeX (\sin{x}, \frac{n}{d}, \sqrt{x})
//   * octave_code(expr) — MATLAB / Octave (element-wise operators)
//   * rust_code(expr)  — Rust (f64 literals, .powi/.powf, method-call funcs)
//   * julia_code(expr) — Julia (^ for power, pi, MathConstants.e)
//   * glsl_code(expr)  — GLSL shader (float literals, pow/sqrt, PI/E)
//   * mathml(expr)     — Presentation MathML (<math>, <msup>, <mfrac>, ...)
//   * pretty(expr)     — ASCII 2D pretty form
//   * dot(expr)        — Graphviz DOT digraph of the expression tree
//   * srepr(expr)      — SymPy constructor form (Add(...), Pow(...), ...)
//
// Plus higher-level function emission:
//   c_function(name, expr, args)   — wrap as a C function
//   cxx_function(...)              — wrap as a C++ function
//   fortran_function(...)          — wrap as a Fortran function
//   octave_function(...)           — wrap as an Octave function
//
// Reference: sympy/printing/{c, cxx, fortran, latex, octave, pycode,
//            codeprinter, precedence, pretty}.py

#include <string>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp::printing {

// --- Single-expression printers ---
[[nodiscard]] SYMPP_EXPORT std::string ccode(const Expr& e);
[[nodiscard]] SYMPP_EXPORT std::string cxxcode(const Expr& e);
[[nodiscard]] SYMPP_EXPORT std::string fcode(const Expr& e);
[[nodiscard]] SYMPP_EXPORT std::string latex(const Expr& e);
[[nodiscard]] SYMPP_EXPORT std::string octave_code(const Expr& e);
[[nodiscard]] SYMPP_EXPORT std::string rust_code(const Expr& e);
[[nodiscard]] SYMPP_EXPORT std::string julia_code(const Expr& e);
[[nodiscard]] SYMPP_EXPORT std::string glsl_code(const Expr& e);
[[nodiscard]] SYMPP_EXPORT std::string mathml(const Expr& e);
[[nodiscard]] SYMPP_EXPORT std::string pretty(const Expr& e);

// Graphviz DOT digraph of the expression tree (one node per subexpression,
// directed parent->child edges). Deterministic output.
[[nodiscard]] SYMPP_EXPORT std::string dot(const Expr& e);

// SymPy's srepr: the constructor form, e.g. srepr(x+1) yields
// "Add(Symbol('x'), Integer(1))".
[[nodiscard]] SYMPP_EXPORT std::string srepr(const Expr& e);

// --- Function emission ---
// args is the parameter list (each arg should be a Symbol Expr).
[[nodiscard]] SYMPP_EXPORT std::string c_function(
    const std::string& name, const Expr& body,
    const std::vector<Expr>& args);
[[nodiscard]] SYMPP_EXPORT std::string cxx_function(
    const std::string& name, const Expr& body,
    const std::vector<Expr>& args);
[[nodiscard]] SYMPP_EXPORT std::string fortran_function(
    const std::string& name, const Expr& body,
    const std::vector<Expr>& args);
[[nodiscard]] SYMPP_EXPORT std::string octave_function(
    const std::string& name, const Expr& body,
    const std::vector<Expr>& args);

}  // namespace sympp::printing
