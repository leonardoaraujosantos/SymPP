// sqrt / Abs / sign tests.
//
// References:
//   sympy/functions/elementary/tests/test_miscellaneous.py — sqrt
//   sympy/functions/elementary/tests/test_complexes.py     — Abs, sign

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/functions/miscellaneous.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- sqrt ------------------------------------------------------------------

TEST_CASE("sqrt: identities", "[3d][sqrt]") {
    REQUIRE(sqrt(S::Zero()) == S::Zero());
    REQUIRE(sqrt(S::One()) == S::One());
}

TEST_CASE("sqrt: perfect squares fold to integers", "[3d][sqrt]") {
    // Reference: sympy/sqrt(4) == 2, sqrt(81) == 9.
    REQUIRE(sqrt(integer(4)) == integer(2));
    REQUIRE(sqrt(integer(81)) == integer(9));
    REQUIRE(sqrt(integer(100)) == integer(10));
}

TEST_CASE("sqrt: non-square integer stays as Pow(_, 1/2)", "[3d][sqrt]") {
    auto e = sqrt(integer(2));
    REQUIRE(e->type_id() == TypeId::Pow);
    REQUIRE(e->str() == "2**(1/2)");
}

TEST_CASE("sqrt: symbolic stays as Pow(_, 1/2)", "[3d][sqrt]") {
    auto x = symbol("x");
    auto e = sqrt(x);
    REQUIRE(e->type_id() == TypeId::Pow);
    REQUIRE(e->str() == "x**(1/2)");
}

TEST_CASE("sqrt: x^2 is positive — assumption-aware", "[3d][sqrt][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sqrt(pow(x, integer(2)));
    // No simplification yet (sqrt(x^2) = |x| comes via simplify, Phase 5);
    // for now it's just (x^2)^(1/2). Verify oracle agrees structurally.
    REQUIRE(oracle.equivalent(e->str(), "sqrt(x**2)"));
}

// ----- Abs -------------------------------------------------------------------

TEST_CASE("Abs: numeric reductions", "[3d][abs]") {
    REQUIRE(abs(integer(5)) == integer(5));
    REQUIRE(abs(integer(-5)) == integer(5));
    REQUIRE(abs(integer(0)) == integer(0));
    REQUIRE(abs(rational(-3, 4)) == rational(3, 4));
}

TEST_CASE("Abs: positive symbol → x", "[3d][abs][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    REQUIRE(abs(x) == x);
}

TEST_CASE("Abs: negative symbol → -x", "[3d][abs][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_negative(true));
    REQUIRE(abs(x) == mul(S::NegativeOne(), x));
}

TEST_CASE("Abs(-x) = Abs(x)", "[3d][abs]") {
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    REQUIRE(abs(neg) == abs(x));
}

TEST_CASE("Abs: stays unevaluated for unknown sign", "[3d][abs]") {
    auto x = symbol("x");
    auto e = abs(x);
    REQUIRE(e->type_id() == TypeId::Function);
    REQUIRE(e->str() == "Abs(x)");
}

TEST_CASE("Abs: predicates always real-nonneg", "[3d][abs][assumptions]") {
    auto x = symbol("x");
    auto e = abs(x);
    REQUIRE(is_real(e) == true);
    REQUIRE(is_nonnegative(e) == true);
    REQUIRE(is_negative(e) == false);
}

TEST_CASE("Abs: positive when arg known nonzero", "[3d][abs][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_positive(true));  // implies nonzero
    REQUIRE(is_positive(abs(x)) == true);
}

// ----- sign ------------------------------------------------------------------

TEST_CASE("sign: numeric reductions", "[3d][sign]") {
    REQUIRE(sign(integer(0)) == integer(0));
    REQUIRE(sign(integer(7)) == integer(1));
    REQUIRE(sign(integer(-3)) == integer(-1));
    REQUIRE(sign(rational(-1, 2)) == integer(-1));
}

TEST_CASE("sign: positive symbol → 1", "[3d][sign][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    REQUIRE(sign(x) == integer(1));
}

TEST_CASE("sign: negative symbol → -1", "[3d][sign][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_negative(true));
    REQUIRE(sign(x) == integer(-1));
}

TEST_CASE("sign(-x) = -sign(x)", "[3d][sign]") {
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    auto e = sign(neg);
    auto expected = mul(S::NegativeOne(), sign(x));
    REQUIRE(e == expected);
}

TEST_CASE("sign(Pi) = 1", "[3d][sign]") {
    REQUIRE(sign(S::Pi()) == integer(1));
}

// ----- Substitution ----------------------------------------------------------

TEST_CASE("Abs: subs propagates", "[3d][abs][subs]") {
    auto x = symbol("x");
    auto e = abs(x) + integer(1);
    REQUIRE(subs(e, x, integer(-7)) == integer(8));
}

TEST_CASE("sqrt: subs into sqrt(x) → exact integer", "[3d][sqrt][subs]") {
    auto x = symbol("x");
    auto e = sqrt(x);
    REQUIRE(subs(e, x, integer(9)) == integer(3));
}

// ----- Oracle structural ----------------------------------------------------

TEST_CASE("sqrt/abs/sign: structural output matches SymPy",
          "[3d][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(sqrt(integer(4))->str(), "2"));
    REQUIRE(oracle.equivalent(abs(integer(-5))->str(), "5"));
    REQUIRE(oracle.equivalent(sign(integer(-3))->str(), "-1"));
    REQUIRE(oracle.equivalent(sqrt(x)->str(), "sqrt(x)"));
    REQUIRE(oracle.equivalent(abs(x)->str(), "Abs(x)"));
    REQUIRE(oracle.equivalent(sign(x)->str(), "sign(x)"));
}
