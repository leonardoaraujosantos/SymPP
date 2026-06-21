// LaTeX parser — round-trips with printing::latex() and parses directly.

#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/parsing/latex_parser.hpp>
#include <sympp/printing/printing.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {
bool roundtrips(const Expr& e) {
    Expr back = parse_latex(printing::latex(e));
    return sympp::testing::Oracle::instance().equivalent(back->str(), e->str());
}
}  // namespace

TEST_CASE("LaTeX round-trip: parse(latex(e)) == e", "[latex][oracle]") {
    auto x = symbol("x"), y = symbol("y");
    std::vector<Expr> exprs = {
        x + integer(1),
        pow(x, integer(2)),
        x * y,
        (x + integer(1)) * pow(y, integer(-1)),
        sqrt(x),
        sin(x) + cos(y),
        pow(x, integer(2)) * rational(1, 2) + integer(1),
        pow(x, integer(3)) - integer(2) * x + integer(1),
        exp(x) * sin(y),
        tan(x + y),
    };
    for (const auto& e : exprs) REQUIRE(roundtrips(e));
}

TEST_CASE("LaTeX parser: direct parses", "[latex]") {
    auto x = symbol("x"), y = symbol("y");
    REQUIRE(parse_latex("x + 1") == x + integer(1));
    REQUIRE(parse_latex("\\frac{x + 1}{y}") == (x + integer(1)) * pow(y, integer(-1)));
    REQUIRE(parse_latex("\\sqrt{x}") == sqrt(x));
    REQUIRE(parse_latex("x^{2}") == pow(x, integer(2)));
    REQUIRE(parse_latex("2 x") == integer(2) * x);          // implicit multiplication
    REQUIRE(parse_latex("\\sin(x) + \\cos(y)") == sin(x) + cos(y));
    REQUIRE(parse_latex("\\pi") == S::Pi());
    REQUIRE(parse_latex("-x") == mul(integer(-1), x));
}

TEST_CASE("LaTeX parser: errors", "[latex]") {
    REQUIRE_THROWS_AS(parse_latex("x +"), std::runtime_error);
    REQUIRE_THROWS_AS(parse_latex("\\frac{1}"), std::runtime_error);  // missing 2nd group
}
