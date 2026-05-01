// UndefinedFunction tests — user-declared anonymous f(x) style.
//
// Reference: sympy/core/tests/test_function.py — `Function('f')(x)` semantics.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/undefined_function.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

TEST_CASE("UndefinedFunction: f(x) prints as 'f(x)'", "[3a][function]") {
    auto f = function_symbol("f");
    auto x = symbol("x");
    auto e = f(x);
    REQUIRE(e->str() == "f(x)");
    REQUIRE(e->type_id() == TypeId::Function);
}

TEST_CASE("UndefinedFunction: f(x) and f(x) are interned", "[3a][function]") {
    auto f = function_symbol("f");
    auto x = symbol("x");
    auto a = f(x);
    auto b = f(x);
    REQUIRE(a.get() == b.get());
}

TEST_CASE("UndefinedFunction: f(x) and g(x) compare unequal",
          "[3a][function]") {
    auto f = function_symbol("f");
    auto g = function_symbol("g");
    auto x = symbol("x");
    REQUIRE_FALSE(f(x) == g(x));
}

TEST_CASE("UndefinedFunction: multi-arg f(x, y)", "[3a][function]") {
    auto f = function_symbol("f");
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = f(x, y);
    REQUIRE(e->str() == "f(x, y)");
    REQUIRE(e->args().size() == 2);
}

TEST_CASE("UndefinedFunction: subs replaces inside f(x)", "[3a][function]") {
    auto f = function_symbol("f");
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = f(x);
    auto r = subs(e, x, y);
    REQUIRE(r == f(y));
    REQUIRE(r->str() == "f(y)");
}

TEST_CASE("UndefinedFunction: subs replaces f(x) entirely", "[3a][function]") {
    auto f = function_symbol("f");
    auto x = symbol("x");
    // Substitute the whole f(x) with a new value.
    auto target = f(x);
    auto e = target + integer(1);
    auto r = subs(e, target, integer(5));
    REQUIRE(r == integer(6));
}

TEST_CASE("UndefinedFunction: free_symbols sees through f(x, y)",
          "[3a][function]") {
    auto f = function_symbol("f");
    auto x = symbol("x");
    auto y = symbol("y");
    auto fs = free_symbols(f(x, y));
    REQUIRE(fs.size() == 2);
    REQUIRE(fs.count(x) == 1);
    REQUIRE(fs.count(y) == 1);
}

TEST_CASE("UndefinedFunction: appears inside compound expressions",
          "[3a][function]") {
    auto& oracle = Oracle::instance();
    auto f = function_symbol("f");
    auto x = symbol("x");
    auto e = integer(2) * f(x) + integer(3);
    REQUIRE(oracle.equivalent(e->str(), "2*f(x) + 3"));
}
