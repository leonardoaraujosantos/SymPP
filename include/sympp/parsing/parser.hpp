#pragma once

// Phase 15 — Parser / serialization.
//
// parse(s) converts a math-syntax string into an Expr. Useful for
// reading expressions from non-C++ sources (config files, network
// protocols, REPLs) and for round-tripping with SymPy.
//
// Grammar (right-associative ** as in Python/SymPy):
//
//     expr      = add
//     add       = mul (('+' | '-') mul)*
//     mul       = pow_e (('*' | '/') pow_e)*
//     pow_e     = unary ('**' pow_e | '^' pow_e)?       # right-assoc
//     unary     = ('+' | '-')? atom
//     atom      = number
//               | identifier '(' args ')'
//               | identifier
//               | '(' expr ')'
//     args      = expr (',' expr)*
//
// Recognized identifiers (constants):
//     pi, E, I, True, False, oo (Infinity — not yet shipped, treated
//     as a Symbol)
//
// Recognized function names (matched against the SymPP factory set):
//     sin, cos, tan, asin, acos, atan, atan2,
//     sinh, cosh, tanh, asinh, acosh, atanh,
//     exp, log, sqrt, abs, sign, floor, ceiling,
//     factorial, binomial, gamma, loggamma,
//     erf, erfc, heaviside, dirac_delta
//
// Anything else: an UndefinedFunction (e.g. parse("f(x)") yields
// f(x) where f is a fresh function symbol).
//
// Throws sympp::parsing::ParseError on bad input.
//
// Reference: sympy/parsing/sympy_parser.py (algorithm only — we
// build a real recursive-descent parser, not Python's eval-based
// approach).

#include <stdexcept>
#include <string>
#include <string_view>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp::parsing {

class SYMPP_EXPORT ParseError : public std::runtime_error {
public:
    ParseError(std::string msg, std::size_t pos)
        : std::runtime_error(std::move(msg)), pos_(pos) {}
    [[nodiscard]] std::size_t position() const noexcept { return pos_; }
private:
    std::size_t pos_;
};

[[nodiscard]] SYMPP_EXPORT Expr parse(std::string_view source);

}  // namespace sympp::parsing
