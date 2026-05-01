// Laplace transform tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/integrals/transforms.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

TEST_CASE("laplace: 1 -> 1/s", "[8][laplace]") {
    auto t = symbol("t");
    auto s = symbol("s");
    auto F = laplace_transform(integer(1), t, s);
    REQUIRE(F == pow(s, S::NegativeOne()));
}

TEST_CASE("laplace: t -> 1/s^2", "[8][laplace]") {
    auto t = symbol("t");
    auto s = symbol("s");
    auto F = laplace_transform(t, t, s);
    REQUIRE(F == pow(s, integer(-2)));
}

TEST_CASE("laplace: exp(a*t) -> 1/(s - a)", "[8][laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    auto a = symbol("a");
    auto F = laplace_transform(exp(a * t), t, s);
    REQUIRE(oracle.equivalent(F->str(), "1/(s - a)"));
}

TEST_CASE("laplace: sin(a*t) -> a/(s^2 + a^2)", "[8][laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    auto a = symbol("a");
    auto F = laplace_transform(sin(a * t), t, s);
    REQUIRE(oracle.equivalent(F->str(), "a/(s**2 + a**2)"));
}

TEST_CASE("laplace: cos(a*t) -> s/(s^2 + a^2)", "[8][laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    auto a = symbol("a");
    auto F = laplace_transform(cos(a * t), t, s);
    REQUIRE(oracle.equivalent(F->str(), "s/(s**2 + a**2)"));
}

TEST_CASE("laplace: linearity over Add", "[8][laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    // L{2 + 3*t + sin(t)} = 2/s + 3/s^2 + 1/(s^2 + 1)
    auto f = integer(2) + integer(3) * t + sin(t);
    auto F = laplace_transform(f, t, s);
    REQUIRE(oracle.equivalent(F->str(), "2/s + 3/s**2 + 1/(s**2 + 1)"));
}
