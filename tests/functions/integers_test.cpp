// floor / ceiling tests.
//
// Reference: sympy/functions/elementary/tests/test_integers.py
//   * test_floor, test_ceiling

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/functions/integers.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- floor -----------------------------------------------------------------

TEST_CASE("floor: integer is identity", "[3g][floor]") {
    REQUIRE(floor(integer(0)) == integer(0));
    REQUIRE(floor(integer(7)) == integer(7));
    REQUIRE(floor(integer(-3)) == integer(-3));
}

TEST_CASE("floor: rational rounds toward -inf", "[3g][floor]") {
    REQUIRE(floor(rational(7, 2)) == integer(3));     // 3.5 -> 3
    REQUIRE(floor(rational(-7, 2)) == integer(-4));   // -3.5 -> -4
    REQUIRE(floor(rational(1, 3)) == integer(0));     // 0.333.. -> 0
    REQUIRE(floor(rational(-1, 3)) == integer(-1));   // -0.333 -> -1
}

TEST_CASE("floor: float matches mpfr_floor", "[3g][floor]") {
    auto e = floor(float_value(2.7));
    REQUIRE(e->type_id() == TypeId::Float);
}

TEST_CASE("floor: integer-tagged symbol is identity",
          "[3g][floor][assumptions]") {
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    REQUIRE(floor(n) == n);
}

TEST_CASE("floor: stays unevaluated on a generic symbol", "[3g][floor]") {
    auto x = symbol("x");
    REQUIRE(floor(x)->str() == "floor(x)");
}

// ----- ceiling ---------------------------------------------------------------

TEST_CASE("ceiling: integer is identity", "[3g][ceiling]") {
    REQUIRE(ceiling(integer(7)) == integer(7));
    REQUIRE(ceiling(integer(-3)) == integer(-3));
}

TEST_CASE("ceiling: rational rounds toward +inf", "[3g][ceiling]") {
    REQUIRE(ceiling(rational(7, 2)) == integer(4));    // 3.5 -> 4
    REQUIRE(ceiling(rational(-7, 2)) == integer(-3));  // -3.5 -> -3
    REQUIRE(ceiling(rational(1, 3)) == integer(1));
    REQUIRE(ceiling(rational(-1, 3)) == integer(0));
}

TEST_CASE("ceiling: float matches mpfr_ceil", "[3g][ceiling]") {
    auto e = ceiling(float_value(2.3));
    REQUIRE(e->type_id() == TypeId::Float);
}

TEST_CASE("ceiling: stays unevaluated on a generic symbol", "[3g][ceiling]") {
    auto x = symbol("x");
    REQUIRE(ceiling(x)->str() == "ceiling(x)");
}

// ----- Substitution ----------------------------------------------------------

TEST_CASE("floor/ceiling: subs collapses on numeric substitution",
          "[3g][subs]") {
    auto x = symbol("x");
    REQUIRE(subs(floor(x), x, rational(11, 3)) == integer(3));
    REQUIRE(subs(ceiling(x), x, rational(11, 3)) == integer(4));
}

// ----- Oracle structural -----------------------------------------------------

TEST_CASE("floor / ceiling: structural and value match SymPy",
          "[3g][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(floor(x)->str(), "floor(x)"));
    REQUIRE(oracle.equivalent(ceiling(x)->str(), "ceiling(x)"));

    REQUIRE(oracle.equivalent(floor(rational(7, 2))->str(), "3"));
    REQUIRE(oracle.equivalent(ceiling(rational(7, 2))->str(), "4"));
    REQUIRE(oracle.equivalent(floor(rational(-7, 2))->str(), "-4"));
    REQUIRE(oracle.equivalent(ceiling(rational(-7, 2))->str(), "-3"));
}
