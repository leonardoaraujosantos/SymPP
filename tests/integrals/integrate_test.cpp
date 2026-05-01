// integrate tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/integrals/integrate.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Polynomials ----------------------------------------------------------

TEST_CASE("integrate: constant", "[7][integrate]") {
    auto x = symbol("x");
    REQUIRE(integrate(integer(5), x) == integer(5) * x);
}

TEST_CASE("integrate: x", "[7][integrate]") {
    auto x = symbol("x");
    auto r = integrate(x, x);
    REQUIRE(r == mul(rational(1, 2), pow(x, integer(2))));
}

TEST_CASE("integrate: x^n", "[7][integrate][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(x, integer(3)), x);
    // x^4 / 4
    REQUIRE(oracle.equivalent(r->str(), "x**4 / 4"));
}

TEST_CASE("integrate: 1/x = log(x)", "[7][integrate]") {
    auto x = symbol("x");
    auto r = integrate(pow(x, S::NegativeOne()), x);
    REQUIRE(r == log(x));
}

TEST_CASE("integrate: linearity", "[7][integrate][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // ∫ (3*x^2 + 2*x + 1) dx = x^3 + x^2 + x
    auto e = integer(3) * pow(x, integer(2)) + integer(2) * x + integer(1);
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "x**3 + x**2 + x"));
}

// ----- exp / sin / cos ------------------------------------------------------

TEST_CASE("integrate: exp(x)", "[7][integrate]") {
    auto x = symbol("x");
    REQUIRE(integrate(exp(x), x) == exp(x));
}

TEST_CASE("integrate: sin(x) = -cos(x)", "[7][integrate]") {
    auto x = symbol("x");
    REQUIRE(integrate(sin(x), x) == mul(S::NegativeOne(), cos(x)));
}

TEST_CASE("integrate: cos(x) = sin(x)", "[7][integrate]") {
    auto x = symbol("x");
    REQUIRE(integrate(cos(x), x) == sin(x));
}

TEST_CASE("integrate: exp(2*x) = exp(2*x) / 2",
          "[7][integrate][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(exp(integer(2) * x), x);
    REQUIRE(oracle.equivalent(r->str(), "exp(2*x) / 2"));
}

TEST_CASE("integrate: sin(3*x) = -cos(3*x)/3",
          "[7][integrate][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(sin(integer(3) * x), x);
    REQUIRE(oracle.equivalent(r->str(), "-cos(3*x) / 3"));
}

// ----- Definite integration -------------------------------------------------

TEST_CASE("integrate: definite x from 0 to 1 = 1/2",
          "[7][integrate]") {
    auto x = symbol("x");
    auto r = integrate(x, x, integer(0), integer(1));
    REQUIRE(r == rational(1, 2));
}

TEST_CASE("integrate: definite x^2 from 0 to 3 = 9",
          "[7][integrate]") {
    auto x = symbol("x");
    auto r = integrate(pow(x, integer(2)), x, integer(0), integer(3));
    REQUIRE(r == integer(9));
}

TEST_CASE("integrate: definite cos(x) from 0 to pi = 0",
          "[7][integrate]") {
    auto x = symbol("x");
    auto r = integrate(cos(x), x, S::Zero(), S::Pi());
    REQUIRE(r == integer(0));
}

// ----- Round-trip property: diff(integrate(f, x), x) == f --------------------

TEST_CASE("integrate / diff round-trip on polynomial",
          "[7][integrate][diff]") {
    auto x = symbol("x");
    auto f = pow(x, integer(3)) + integer(2) * x + integer(7);
    auto F = integrate(f, x);
    auto fp = diff(F, x);
    REQUIRE(simplify(fp - f) == S::Zero());
}
