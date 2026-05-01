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

// ============================================================================
// Phase 3h additions
// ============================================================================

// ----- re / im / conjugate / arg --------------------------------------------

TEST_CASE("re: real arg passes through", "[3h][re]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(re(x) == x);
    REQUIRE(re(integer(7)) == integer(7));
    REQUIRE(re(rational(-1, 2)) == rational(-1, 2));
}

TEST_CASE("im: real arg yields zero", "[3h][im]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(im(x) == S::Zero());
    REQUIRE(im(integer(7)) == S::Zero());
}

TEST_CASE("conjugate: real arg passes through; involution otherwise",
          "[3h][conjugate]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(conjugate(x) == x);

    auto y = symbol("y");  // unknown reality
    auto e = conjugate(conjugate(y));
    REQUIRE(e == y);
}

TEST_CASE("arg_: positive -> 0; negative -> pi", "[3h][arg]") {
    auto p = symbol("p", AssumptionMask{}.set_positive(true));
    auto n = symbol("n", AssumptionMask{}.set_negative(true));
    REQUIRE(arg_(p) == S::Zero());
    REQUIRE(arg_(n) == S::Pi());
    REQUIRE(arg_(integer(5)) == S::Zero());
}

TEST_CASE("complex parts: stay unevaluated for unknown arg", "[3h][complex]") {
    auto x = symbol("x");  // unknown reality
    REQUIRE(re(x)->type_id() == TypeId::Function);
    REQUIRE(im(x)->type_id() == TypeId::Function);
    REQUIRE(arg_(x)->type_id() == TypeId::Function);
}

TEST_CASE("re/im are real", "[3h][complex][assumptions]") {
    auto x = symbol("x");
    REQUIRE(is_real(re(x)) == true);
    REQUIRE(is_real(im(x)) == true);
}

// ----- Min, Max --------------------------------------------------------------

TEST_CASE("min/max: numeric inputs collapse", "[3h][minmax]") {
    REQUIRE(min(integer(2), integer(5)) == integer(2));
    REQUIRE(max(integer(2), integer(5)) == integer(5));
    REQUIRE(min({integer(7), integer(3), integer(11), integer(1)}) == integer(1));
    REQUIRE(max({integer(7), integer(3), integer(11), integer(1)}) == integer(11));
}

TEST_CASE("min/max: rational inputs", "[3h][minmax]") {
    REQUIRE(min(rational(1, 2), rational(1, 3)) == rational(1, 3));
    REQUIRE(max(rational(1, 2), rational(1, 3)) == rational(1, 2));
    REQUIRE(min(rational(-1, 2), integer(0)) == rational(-1, 2));
}

TEST_CASE("min/max: mixed numeric+symbolic keeps symbolics",
          "[3h][minmax]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = min({x, integer(2), integer(7), y});
    REQUIRE(e->type_id() == TypeId::Function);
    // After collection, expect 3 args: 2 (the numeric min), x, y.
    REQUIRE(e->args().size() == 3);
}

TEST_CASE("min/max: duplicate symbolics deduplicate",
          "[3h][minmax]") {
    auto x = symbol("x");
    auto e = max(x, x);
    REQUIRE(e == x);
}

TEST_CASE("min/max: real-arg propagation", "[3h][minmax][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    auto y = symbol("y", AssumptionMask{}.set_real(true));
    REQUIRE(is_real(min(x, y)) == true);
    REQUIRE(is_real(max(x, y)) == true);
}

TEST_CASE("min/max: oracle structural agreement", "[3h][minmax][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(min(integer(3), integer(7))->str(), "3"));
    REQUIRE(oracle.equivalent(max(integer(3), integer(7))->str(), "7"));
}
