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
#include <sympp/functions/special.hpp>
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

// ----- Canonical (SymPy) function-name parsing (regression, PARSE-1) ---------
// str() emits Abs/Heaviside/DiracDelta (SymPy's canonical capitalised names),
// but the parser only accepted the lowercase aliases — so parse(e->str())
// failed to round-trip and `Abs(-3)` (valid SymPy) stayed an undefined function.
TEST_CASE("parse: capitalised SymPy names (Abs/Heaviside/DiracDelta)",
          "[15][parser][regression]") {
    REQUIRE(parse("Abs(-3)") == integer(3));
    // Lowercase aliases still accepted.
    REQUIRE(parse("abs(-3)") == integer(3));
}

TEST_CASE("parse: round-trips Abs/Heaviside/DiracDelta str()",
          "[15][parser][roundtrip][regression]") {
    auto x = symbol("x");
    REQUIRE(parse(abs(x)->str()) == abs(x));
    REQUIRE(parse(heaviside(x)->str()) == heaviside(x));
    REQUIRE(parse(dirac_delta(x)->str()) == dirac_delta(x));
}

// ----- Min/Max canonical-name parsing (regression, PARSE-2) ------------------
// str() emits Min/Max (SymPy's names); the parser had no entry, so parse of a
// printed Min/Max produced an undefined function and broke round-trip.
TEST_CASE("parse: binary Min/Max round-trip and evaluate",
          "[15][parser][regression]") {
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(parse(min(x, y)->str()) == min(x, y));
    REQUIRE(parse(max(x, y)->str()) == max(x, y));
    REQUIRE(parse("Min(3, 5)") == integer(3));
    REQUIRE(parse("Max(3, 5)") == integer(5));
}

// ----- Variadic Min/Max parsing (regression, PARSE-3) ------------------------
// The parser only dispatched 2-argument Min/Max; 3+-argument calls fell through
// to an undefined function (Max(3, 7, 1) stayed unevaluated). Route any arity
// to the variadic min()/max(), which fold the numeric args.
TEST_CASE("parse: variadic Min/Max fold and round-trip",
          "[15][parser][regression]") {
    auto x = symbol("x");
    REQUIRE(parse("Max(3, 7, 1)") == integer(7));
    REQUIRE(parse("Min(3, 7, 1)") == integer(1));
    REQUIRE(parse("Max(1, 2, 3, 4)") == integer(4));
    REQUIRE(parse("Min(1, 2, 3, 4)") == integer(1));
    // Mixed symbolic + numeric: numeric args collapse to one extreme.
    REQUIRE(parse("Max(x, 3, 1)") == max({x, integer(3), integer(1)}));
    REQUIRE(parse(max({x, integer(3), integer(1)})->str())
            == max({x, integer(3), integer(1)}));
}
