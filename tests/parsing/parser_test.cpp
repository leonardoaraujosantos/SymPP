// Phase 15 — parser tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/parsing/parser.hpp>

using namespace sympp;
using namespace sympp::parsing;

// ----- Numbers --------------------------------------------------------------

TEST_CASE("parse: integer literal", "[15][parser]") {
    REQUIRE(parse("42") == integer(42));
    REQUIRE(parse("-7") == integer(-7));
    REQUIRE(parse("0") == integer(0));
}

TEST_CASE("parse: large integer literal", "[15][parser]") {
    REQUIRE(parse("18446744073709551616") == integer("18446744073709551616"));
}

TEST_CASE("parse: float literal", "[15][parser]") {
    auto e = parse("3.14");
    REQUIRE(e->type_id() == TypeId::Float);
}

TEST_CASE("parse: scientific notation", "[15][parser]") {
    auto e = parse("1.5e-3");
    REQUIRE(e->type_id() == TypeId::Float);
}

// ----- Symbols & constants --------------------------------------------------

TEST_CASE("parse: bare identifier becomes Symbol", "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("x") == x);
}

TEST_CASE("parse: pi recognized as constant", "[15][parser]") {
    REQUIRE(parse("pi") == S::Pi());
    REQUIRE(parse("Pi") == S::Pi());
}

TEST_CASE("parse: E recognized", "[15][parser]") {
    REQUIRE(parse("E") == S::E());
}

TEST_CASE("parse: I recognized as imaginary unit", "[15][parser]") {
    REQUIRE(parse("I") == S::I());
}

// ----- Arithmetic and precedence -------------------------------------------

TEST_CASE("parse: addition", "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("x + 2") == x + integer(2));
}

TEST_CASE("parse: subtraction", "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("x - 2") == x - integer(2));
}

TEST_CASE("parse: multiplication", "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("3 * x") == integer(3) * x);
}

TEST_CASE("parse: division", "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("x / 2") == x / integer(2));
}

TEST_CASE("parse: precedence — multiplication over addition",
          "[15][parser]") {
    auto x = symbol("x");
    auto y = symbol("y");
    // 1 + 2 * 3 → 1 + (2*3) = 7, not (1+2)*3 = 9
    REQUIRE(parse("1 + 2 * 3") == integer(7));
    // x + 2*y == (x) + (2*y), Add(x, 2y)
    REQUIRE(parse("x + 2*y") == x + integer(2) * y);
}

TEST_CASE("parse: precedence — power over multiplication",
          "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("2 * x**2") == integer(2) * pow(x, integer(2)));
}

TEST_CASE("parse: power is right-associative",
          "[15][parser]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    // x**y**z = x ** (y ** z), not (x ** y) ** z
    REQUIRE(parse("x ** y ** z") == pow(x, pow(y, z)));
}

TEST_CASE("parse: ^ is alias for **", "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("x^3") == pow(x, integer(3)));
}

TEST_CASE("parse: parentheses override precedence",
          "[15][parser]") {
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(parse("(x + y) * 2") == (x + y) * integer(2));
}

TEST_CASE("parse: unary minus", "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("-x") == mul(S::NegativeOne(), x));
    REQUIRE(parse("--x") == x);  // double negation
}

TEST_CASE("parse: unary plus is identity",
          "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("+x") == x);
}

// ----- Function calls --------------------------------------------------------

TEST_CASE("parse: sin(x)", "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("sin(x)") == sin(x));
}

TEST_CASE("parse: cos(x), tan(x), exp(x), log(x), sqrt(x)",
          "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("cos(x)") == cos(x));
    REQUIRE(parse("tan(x)") == tan(x));
    REQUIRE(parse("exp(x)") == exp(x));
    REQUIRE(parse("log(x)") == log(x));
    REQUIRE(parse("sqrt(x)") == sqrt(x));
}

TEST_CASE("parse: nested function calls", "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("sin(cos(x))") == sin(cos(x)));
}

TEST_CASE("parse: function with arithmetic argument",
          "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("sin(2*x + 1)") == sin(integer(2) * x + integer(1)));
}

TEST_CASE("parse: two-arg function atan2", "[15][parser]") {
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(parse("atan2(y, x)") == atan2(y, x));
}

TEST_CASE("parse: undefined function preserved",
          "[15][parser]") {
    auto x = symbol("x");
    auto e = parse("f(x)");
    REQUIRE(e->type_id() == TypeId::Function);
    REQUIRE(e->str() == "f(x)");
}

// ----- Real-world expressions ----------------------------------------------

TEST_CASE("parse: x^2 + sin(x) - 1", "[15][parser]") {
    auto x = symbol("x");
    auto e = parse("x^2 + sin(x) - 1");
    auto expected = pow(x, integer(2)) + sin(x) - integer(1);
    REQUIRE(e == expected);
}

TEST_CASE("parse: rational expression x/(y+1)",
          "[15][parser]") {
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(parse("x / (y + 1)") == x / (y + integer(1)));
}

TEST_CASE("parse: composite — exp(-x^2/2)",
          "[15][parser]") {
    auto x = symbol("x");
    auto e = parse("exp(-x^2/2)");
    auto expected = exp(mul(S::NegativeOne(), pow(x, integer(2))) / integer(2));
    REQUIRE(e == expected);
}

// ----- Whitespace handling --------------------------------------------------

TEST_CASE("parse: tolerates whitespace", "[15][parser]") {
    auto x = symbol("x");
    REQUIRE(parse("  x   +    2  ") == x + integer(2));
    REQUIRE(parse("\t\nsin( x )\n") == sin(x));
}

// ----- Errors ---------------------------------------------------------------

TEST_CASE("parse: unmatched paren throws",
          "[15][parser][error]") {
    REQUIRE_THROWS_AS(parse("(x + 2"), ParseError);
}

TEST_CASE("parse: trailing garbage throws",
          "[15][parser][error]") {
    REQUIRE_THROWS_AS(parse("x + 2 garbage"), ParseError);
}

TEST_CASE("parse: empty input throws",
          "[15][parser][error]") {
    REQUIRE_THROWS_AS(parse(""), ParseError);
}

TEST_CASE("parse: bad operator sequence throws",
          "[15][parser][error]") {
    REQUIRE_THROWS_AS(parse("x +"), ParseError);
}

TEST_CASE("parse: unknown character throws",
          "[15][parser][error]") {
    REQUIRE_THROWS_AS(parse("x @ y"), ParseError);
}

// ----- Round-trip with str() -----------------------------------------------

TEST_CASE("parse: str round-trip on x + 1",
          "[15][parser][roundtrip]") {
    auto x = symbol("x");
    auto e = x + integer(1);
    REQUIRE(parse(e->str()) == e);
}

TEST_CASE("parse: str round-trip on sin(x)",
          "[15][parser][roundtrip]") {
    auto x = symbol("x");
    auto e = sin(x);
    REQUIRE(parse(e->str()) == e);
}

TEST_CASE("parse: str round-trip on 2*x + sin(x*y)",
          "[15][parser][roundtrip]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = integer(2) * x + sin(x * y);
    REQUIRE(parse(e->str()) == e);
}
