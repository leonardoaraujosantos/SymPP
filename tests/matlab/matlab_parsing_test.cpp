// MATLAB facade — string parsing tests.
//
// Validates `str2sym` / `sym(string_view)` route through the
// recursive-descent parser when the input contains expression syntax,
// and fall through to a single Symbol for bare identifiers.

#include <catch2/catch_test_macros.hpp>

#include <stdexcept>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/matlab/matlab.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

TEST_CASE("matlab::str2sym parses polynomial",
          "[15m][matlab][parsing][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = matlab::str2sym("x^2 + 2*x + 1");
    REQUIRE(oracle.equivalent(e->str(), "x**2 + 2*x + 1"));
}

TEST_CASE("matlab::str2sym parses transcendental product",
          "[15m][matlab][parsing][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = matlab::str2sym("sin(x)*exp(-x)");
    REQUIRE(oracle.equivalent(e->str(), "sin(x)*exp(-x)"));
}

TEST_CASE("matlab::str2sym parses parenthesized product",
          "[15m][matlab][parsing][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = matlab::str2sym("(x+1)*(x-1)");
    auto expanded = matlab::expand(e);
    REQUIRE(oracle.equivalent(expanded->str(), "x**2 - 1"));
}

TEST_CASE("matlab::sym bare identifier returns a Symbol",
          "[15m][matlab][parsing]") {
    REQUIRE(matlab::sym("alpha") == symbol("alpha"));
}

TEST_CASE("matlab::sym expression string routes through parser",
          "[15m][matlab][parsing][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = matlab::sym("2*x + 3");
    REQUIRE(oracle.equivalent(e->str(), "2*x + 3"));
}

TEST_CASE("matlab::sym integer literal makes Integer",
          "[15m][matlab][parsing]") {
    REQUIRE(matlab::sym(42L) == integer(42));
}

TEST_CASE("matlab::str2sym throws on malformed input",
          "[15m][matlab][parsing]") {
    REQUIRE_THROWS(matlab::str2sym("x + + + )"));
}

TEST_CASE("matlab::sym round-trips diff(parsed)",
          "[15m][matlab][parsing][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = matlab::syms("x");
    auto e = matlab::sym("x^3 - 2*x");
    auto d = matlab::diff(e, x);
    REQUIRE(oracle.equivalent(d->str(), "3*x**2 - 2"));
}
