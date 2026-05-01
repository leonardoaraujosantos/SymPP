// Phase 0 Hello World: Symbol creation, equality, str(), oracle round-trip.
// Test cases ported from sympy/core/tests/test_symbol.py::test_Symbol.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/symbol.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

TEST_CASE("Symbol: name is preserved", "[0][symbol]") {
    auto x = symbol("x");
    REQUIRE(x->str() == "x");

    auto longname = symbol("alpha_beta");
    REQUIRE(longname->str() == "alpha_beta");
}

TEST_CASE("Symbol: structural equality by name", "[0][symbol]") {
    // Ported from sympy test_Symbol:
    //   x1 = Symbol("x"); x2 = Symbol("x")
    //   a  = Symbol("a")
    //   assert x1 == x2
    //   assert a != x1
    auto a = symbol("a");
    auto x1 = symbol("x");
    auto x2 = symbol("x");

    REQUIRE(x1 == x2);
    REQUIRE_FALSE(a == x1);
    REQUIRE_FALSE(a == x2);
}

TEST_CASE("Symbol: hash is stable for equal names", "[0][symbol]") {
    auto x1 = symbol("x");
    auto x2 = symbol("x");
    REQUIRE(x1->hash() == x2->hash());
}

TEST_CASE("Symbol: is atomic (no args)", "[0][symbol]") {
    auto x = symbol("x");
    REQUIRE(x->is_atomic());
    REQUIRE(x->args().empty());
}

TEST_CASE("Symbol: type tag is Symbol", "[0][symbol]") {
    auto x = symbol("x");
    REQUIRE(x->type_id() == TypeId::Symbol);
}

// ----- Oracle-validated tests --------------------------------------------------

TEST_CASE("Symbol: str() round-trips through SymPy", "[0][symbol][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");

    // Our str() output should parse back to a SymPy Symbol with the same name.
    REQUIRE(oracle.srepr(x->str()) == "Symbol('x')");
    REQUIRE(oracle.sympy_str(x->str()) == "x");
}

TEST_CASE("Symbol: equivalence with SymPy Symbol", "[0][symbol][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");

    // SymPy: Symbol('x') equals Symbol('x')
    REQUIRE(oracle.equivalent(x->str(), "x"));
    // ... and is not equal to a different symbol
    REQUIRE_FALSE(oracle.equivalent(x->str(), "y"));
}
