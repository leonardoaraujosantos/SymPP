// simplify / collect / together tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

TEST_CASE("simplify: trivial identity", "[5][simplify]") {
    auto x = symbol("x");
    auto e = x + x;
    REQUIRE(simplify(e) == integer(2) * x);
}

TEST_CASE("simplify: cancels x - x", "[5][simplify]") {
    auto x = symbol("x");
    REQUIRE(simplify(x - x) == S::Zero());
}

TEST_CASE("simplify: distributes and recombines", "[5][simplify][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // (x + 1)*(x - 1) - x^2 = -1
    auto e = (x + integer(1)) * (x - integer(1)) - pow(x, integer(2));
    auto s = simplify(e);
    REQUIRE(oracle.equivalent(s->str(), "-1"));
}

TEST_CASE("collect: groups powers of var", "[5][collect][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto a = symbol("a");
    auto b = symbol("b");
    // a*x^2 + b*x^2 + x + 3 -> (a + b)*x^2 + x + 3
    auto e = a * pow(x, integer(2)) + b * pow(x, integer(2)) + x + integer(3);
    auto c = collect(e, x);
    REQUIRE(oracle.equivalent(c->str(), "(a + b)*x**2 + x + 3"));
}

TEST_CASE("collect: linear in var", "[5][collect]") {
    auto x = symbol("x");
    auto a = symbol("a");
    auto b = symbol("b");
    // a*x + b*x + 7 -> (a+b)*x + 7
    auto e = a * x + b * x + integer(7);
    auto c = collect(e, x);
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(c->str(), "(a + b)*x + 7"));
}

TEST_CASE("simplify: matches SymPy on polynomial expansion",
          "[5][simplify][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x + integer(1), integer(3)) - pow(x - integer(1), integer(3));
    auto s = simplify(e);
    // Expected: 6*x^2 + 2
    REQUIRE(oracle.equivalent(s->str(), "6*x**2 + 2"));
}
