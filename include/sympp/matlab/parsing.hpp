#pragma once

// MATLAB facade — string parsing.
//
// MATLAB users write `sym('cos(x)+1')` to mean "parse this as an
// expression". A naive sym(string_view) that called symbol(...)
// directly would create a single Symbol whose name is the literal
// string "cos(x)+1" — almost certainly a bug.
//
// Two entry points:
//
//   * str2sym(s)       — always parses, throws on bad input.
//   * sym(string_view) — sniffs `s` for operator / paren /
//                         whitespace characters; if any are present,
//                         routes through str2sym, otherwise falls
//                         through to the single-Symbol path.
//
// Backed by sympp::parsing::parse, which is a real recursive-descent
// parser supporting the full SymPP elementary-function set,
// constants (pi, E, I), operator precedence (** right-assoc, ^
// aliased), and undefined-function preservation.
//
// This header is included by matlab/matlab.hpp; it owns the
// definition of sym(string_view) and sym(long) for the facade.
//
// Reference: sympp/parsing/parser.hpp.

#include <string>
#include <string_view>

#include <sympp/core/api.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/fwd.hpp>
#include <sympp/parsing/parser.hpp>

namespace sympp::matlab {

// str2sym(s) — parse `s` as a math expression. Throws
// sympp::parsing::ParseError on bad input.
[[nodiscard]] inline Expr str2sym(std::string_view source) {
    return parsing::parse(source);
}

namespace detail {

// True iff `s` contains any character that suggests an expression
// rather than a bare identifier (operators, parens, comma,
// whitespace). All of these are illegal in a Symbol name, so the
// check is also a safety net against accidentally minting a Symbol
// named "x + 1".
[[nodiscard]] inline bool looks_like_expression(std::string_view s) noexcept {
    for (char c : s) {
        if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^'
            || c == '(' || c == ')' || c == ','
            || c == ' ' || c == '\t' || c == '\n') {
            return true;
        }
    }
    return false;
}

}  // namespace detail

// sym(name) — bare identifier becomes a Symbol; an expression string
// goes through the parser. Overload set with sym(long) below.
[[nodiscard]] inline Expr sym(std::string_view source) {
    if (detail::looks_like_expression(source)) {
        return str2sym(source);
    }
    return symbol(source);
}

// sym(value) — literal integer constant.
[[nodiscard]] inline Expr sym(long value) {
    return integer(value);
}

}  // namespace sympp::matlab
