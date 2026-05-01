#pragma once

// Phase 13 — Code generation / printing.
//
// Five concrete printers:
//   * ccode(expr)      — C99 syntax with <math.h> functions (sin, cos, ...)
//   * cxxcode(expr)    — C++ with std::sin, std::pow, std::numbers::pi_v<double>
//   * fcode(expr)      — Fortran 90+ (** for power, intrinsic SIN/COS/EXP)
//   * latex(expr)      — TeX (\sin{x}, \frac{n}{d}, \sqrt{x})
//   * octave_code(expr) — MATLAB / Octave (element-wise operators)
//   * pretty(expr)     — ASCII 2D pretty form
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
[[nodiscard]] SYMPP_EXPORT std::string pretty(const Expr& e);

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
