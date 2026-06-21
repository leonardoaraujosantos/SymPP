#pragma once

// Parse a LaTeX math string into an expression — the inverse of
// printing::latex(). Round-trips the dialect the LaTeX printer emits.
//
//   parse_latex("\\frac{x + 1}{y}")        // → (x + 1)/y
//   parse_latex("\\sqrt{x} + \\sin(x)^{2}") // → sqrt(x) + sin(x)**2
//
// Supports: + - · / ^ ^{…}, \frac, \sqrt (and \sqrt[n]), parentheses and
// {…}/\left(\right), the trig/hyperbolic/exp/log commands, \pi, Greek-letter
// commands (→ same-named symbol) and implicit multiplication. Throws
// std::runtime_error on a syntax error.
//
// Reference: sympy/parsing/latex.

#include <string>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

[[nodiscard]] SYMPP_EXPORT Expr parse_latex(const std::string& latex);

}  // namespace sympp
