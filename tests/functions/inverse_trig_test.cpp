// asin / acos / atan / atan2 tests.
//
// Reference: sympy/functions/elementary/tests/test_trigonometric.py
//   * test_asin   (asin(0)=0, asin(1)=pi/2, asin(-1)=-pi/2, odd identity)
//   * test_acos   (acos(0)=pi/2, acos(1)=0, acos(-1)=pi)
//   * test_atan   (atan(0)=0, atan(1)=pi/4, atan(-1)=-pi/4, odd identity)
//   * test_atan2  (atan2(0, x>0) = 0, atan2(y>0, 0) = pi/2, etc.)

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
#include <sympp/calculus/diff.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/parsing/parser.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- asin ------------------------------------------------------------------

TEST_CASE("asin: canonical values", "[3e][asin]") {
    REQUIRE(asin(S::Zero()) == S::Zero());
    REQUIRE(asin(S::One()) == mul(S::Half(), S::Pi()));
    REQUIRE(asin(S::NegativeOne())
            == mul(S::NegativeOne(), mul(S::Half(), S::Pi())));
}

TEST_CASE("asin: odd identity asin(-x) = -asin(x)", "[3e][asin]") {
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    auto e = asin(neg);
    auto expected = mul(S::NegativeOne(), asin(x));
    REQUIRE(e == expected);
}

TEST_CASE("asin: stays unevaluated on a generic symbol", "[3e][asin]") {
    auto x = symbol("x");
    REQUIRE(asin(x)->str() == "asin(x)");
}

// ----- Exact inverse values at special arguments (regression, TRIG-2) --------
// asin/acos/atan at 1/2, √2/2, √3/2, 1, √3, √3/3 used to stay unevaluated
// (only the trivial 0/±1 args folded). Reference: SymPy's exact-value table.
TEST_CASE("asin: exact special arguments", "[3e][asin][regression]") {
    REQUIRE(asin(S::Half()) == mul(rational(1, 6), S::Pi()));
    REQUIRE(asin(mul(S::NegativeOne(), S::Half()))
            == mul(rational(-1, 6), S::Pi()));
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(asin(pow(integer(2), rational(1, 2)) / integer(2))->str(),
                              "pi/4"));
}

TEST_CASE("acos: exact special arguments (via π/2 − asin)",
          "[3e][acos][regression]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(acos(S::Half())->str(), "pi/3"));
    REQUIRE(oracle.equivalent(acos(mul(S::NegativeOne(), S::Half()))->str(),
                              "2*pi/3"));
    REQUIRE(oracle.equivalent(
        acos(pow(integer(3), rational(1, 2)) / integer(2))->str(), "pi/6"));
}

TEST_CASE("atan: exact special arguments", "[3e][atan][regression]") {
    REQUIRE(atan(S::One()) == mul(rational(1, 4), S::Pi()));
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(atan(pow(integer(3), rational(1, 2)))->str(),
                              "pi/3"));
    REQUIRE(oracle.equivalent(
        atan(pow(integer(3), rational(1, 2)) / integer(3))->str(), "pi/6"));
}

TEST_CASE("acos: stays unevaluated on a non-special argument",
          "[3e][acos][regression]") {
    auto x = symbol("x");
    REQUIRE(acos(x)->type_id() == TypeId::Function);
    // 1/3 is not a special cosine value — must not fold.
    REQUIRE(acos(rational(1, 3))->type_id() == TypeId::Function);
}

// ----- f(f⁻¹(x)) = x compositions (regression, FUNC-1) -----------------------
// sin(asin(x)) = x and the cos/tan analogues. The reverse (asin(sin(x))) is NOT
// folded — it is only x on a restricted range, matching SymPy.
TEST_CASE("trig: f(inverse(x)) collapses to x", "[3e][regression]") {
    auto x = symbol("x");
    REQUIRE(sin(asin(x)) == x);
    REQUIRE(cos(acos(x)) == x);
    REQUIRE(tan(atan(x)) == x);
    // Reverse direction stays unevaluated.
    REQUIRE(asin(sin(x))->type_id() == TypeId::Function);
}

// ----- Cross-function inverse compositions (regression, FUNC-2) --------------
// cos(asin(x)) = √(1−x²) and the rest of the trig × inverse-trig table.
TEST_CASE("trig: cross-function inverse compositions",
          "[3e][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(cos(asin(x))->str(), "sqrt(1 - x**2)"));
    REQUIRE(oracle.equivalent(sin(acos(x))->str(), "sqrt(1 - x**2)"));
    REQUIRE(oracle.equivalent(tan(asin(x))->str(), "x/sqrt(1 - x**2)"));
    REQUIRE(oracle.equivalent(cos(atan(x))->str(), "1/sqrt(1 + x**2)"));
    REQUIRE(oracle.equivalent(sin(atan(x))->str(), "x/sqrt(1 + x**2)"));
    REQUIRE(oracle.equivalent(tan(acos(x))->str(), "sqrt(1 - x**2)/x"));
}

// ----- acos ------------------------------------------------------------------

TEST_CASE("acos: canonical values", "[3e][acos]") {
    REQUIRE(acos(S::Zero()) == mul(S::Half(), S::Pi()));
    REQUIRE(acos(S::One()) == S::Zero());
    REQUIRE(acos(S::NegativeOne()) == S::Pi());
}

TEST_CASE("acos: real arg implies real & nonnegative result",
          "[3e][acos][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(is_real(acos(x)) == true);
    REQUIRE(is_nonnegative(acos(x)) == true);
}

// ----- atan ------------------------------------------------------------------

TEST_CASE("atan: canonical values", "[3e][atan]") {
    REQUIRE(atan(S::Zero()) == S::Zero());
    REQUIRE(atan(S::One()) == mul(rational(1, 4), S::Pi()));
    REQUIRE(atan(S::NegativeOne()) == mul(rational(-1, 4), S::Pi()));
}

// Regression (FUNC-INF): atan(±oo) = ±pi/2.
TEST_CASE("atan: values at infinity", "[3e][atan][infinity][regression]") {
    REQUIRE(atan(S::Infinity()) == mul(S::Half(), S::Pi()));
    REQUIRE(atan(S::NegativeInfinity()) == mul(rational(-1, 2), S::Pi()));
}

TEST_CASE("atan: odd identity", "[3e][atan]") {
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    REQUIRE(atan(neg) == mul(S::NegativeOne(), atan(x)));
}

TEST_CASE("atan: real and finite for real arg",
          "[3e][atan][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(is_real(atan(x)) == true);
    REQUIRE(is_finite(atan(x)) == true);
}

// ----- atan2 -----------------------------------------------------------------

TEST_CASE("atan2: zero numerator with positive denominator -> 0",
          "[3e][atan2]") {
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    REQUIRE(atan2(S::Zero(), x) == S::Zero());
}

TEST_CASE("atan2: zero numerator with negative denominator -> pi",
          "[3e][atan2]") {
    auto x = symbol("x", AssumptionMask{}.set_negative(true));
    REQUIRE(atan2(S::Zero(), x) == S::Pi());
}

TEST_CASE("atan2: positive numerator with zero denominator -> pi/2",
          "[3e][atan2]") {
    auto y = symbol("y", AssumptionMask{}.set_positive(true));
    REQUIRE(atan2(y, S::Zero()) == mul(S::Half(), S::Pi()));
}

TEST_CASE("atan2(0, 0) stays unevaluated", "[3e][atan2]") {
    auto e = atan2(S::Zero(), S::Zero());
    REQUIRE(e->type_id() == TypeId::Function);
    REQUIRE(e->str() == "atan2(0, 0)");
}

// Regression (ATAN2-1): atan2 of nonzero numeric args reduced only on the axes;
// the general quadrant cases stayed unevaluated. atan2(y, x) = atan(y/x) with a
// quadrant correction now folds them (and atan folds the special values).
TEST_CASE("atan2: general quadrant reductions match SymPy",
          "[3e][atan2][oracle][regression]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(atan2(integer(1), integer(1))->str(), "pi/4"));
    REQUIRE(oracle.equivalent(atan2(integer(-1), integer(1))->str(), "-pi/4"));
    REQUIRE(oracle.equivalent(atan2(integer(1), integer(-1))->str(), "3*pi/4"));
    REQUIRE(
        oracle.equivalent(atan2(integer(-1), integer(-1))->str(), "-3*pi/4"));
    // First quadrant with an irrational denominator (atan rewrite is faithful
    // even where atan's own table does not fold the argument).
    REQUIRE(oracle.equivalent(atan2(integer(1), sqrt(integer(3)))->str(),
                              "pi/6"));
    // Non-special ratio: stays as atan(2), matching SymPy.
    REQUIRE(oracle.equivalent(atan2(integer(2), integer(1))->str(), "atan(2)"));
}

// ----- Numeric evalf ---------------------------------------------------------

TEST_CASE("asin/acos/atan: numeric Float arg evaluates",
          "[3e][evalf]") {
    REQUIRE(asin(float_value(0.5))->type_id() == TypeId::Float);
    REQUIRE(acos(float_value(0.5))->type_id() == TypeId::Float);
    REQUIRE(atan(float_value(1.0))->type_id() == TypeId::Float);
}

TEST_CASE("asin: high-precision matches SymPy",
          "[3e][asin][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = asin(float_value("0.5", 30));
    auto sympy_value = oracle.evalf("asin(Rational(1, 2))", 30);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

TEST_CASE("atan: high-precision matches SymPy",
          "[3e][atan][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = atan(float_value("1", 30));
    auto sympy_value = oracle.evalf("atan(Integer(1))", 30);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

// ----- Substitution ----------------------------------------------------------

TEST_CASE("asin/acos/atan: subs propagates", "[3e][subs]") {
    auto x = symbol("x");
    REQUIRE(subs(asin(x), x, S::One()) == mul(S::Half(), S::Pi()));
    REQUIRE(subs(acos(x), x, S::One()) == S::Zero());
    REQUIRE(subs(atan(x), x, S::One()) == mul(rational(1, 4), S::Pi()));
}

// ----- Oracle structural -----------------------------------------------------

TEST_CASE("inverse trig: structural output matches SymPy",
          "[3e][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(asin(x)->str(), "asin(x)"));
    REQUIRE(oracle.equivalent(acos(x)->str(), "acos(x)"));
    REQUIRE(oracle.equivalent(atan(x)->str(), "atan(x)"));
    REQUIRE(oracle.equivalent(asin(S::Zero())->str(), "0"));
    REQUIRE(oracle.equivalent(asin(S::One())->str(), "pi/2"));
    REQUIRE(oracle.equivalent(acos(S::NegativeOne())->str(), "pi"));
    REQUIRE(oracle.equivalent(atan(S::One())->str(), "pi/4"));
}

// ----- Inverse reciprocal trio: acot / asec / acsc (AINV-RECIP) --------------

TEST_CASE("acot: canonical values and odd parity", "[3e][acot]") {
    auto pi = S::Pi();
    REQUIRE(acot(S::Zero()) == mul(S::Half(), pi));      // acot(0) = pi/2
    REQUIRE(acot(S::One()) == mul(rational(1, 4), pi));   // acot(1) = pi/4
    REQUIRE(acot(S::NegativeOne()) == mul(rational(-1, 4), pi));  // odd
    REQUIRE(acot(pow(integer(3), rational(1, 2))) == mul(rational(1, 6), pi));
    REQUIRE(acot(S::Infinity()) == S::Zero());            // acot(oo) = 0
    auto x = symbol("x");
    REQUIRE(acot(integer(2))->type_id() == TypeId::Function);  // not a table angle
    REQUIRE(acot(x)->type_id() == TypeId::Function);
}

TEST_CASE("asec: canonical values, pole, no simple parity", "[3e][asec]") {
    auto pi = S::Pi();
    REQUIRE(asec(S::One()) == S::Zero());                 // asec(1) = 0
    REQUIRE(asec(integer(2)) == mul(rational(1, 3), pi));  // asec(2) = pi/3
    REQUIRE(asec(pow(integer(2), rational(1, 2))) == mul(rational(1, 4), pi));
    REQUIRE(asec(S::NegativeOne()) == pi);                // asec(-1) = pi
    REQUIRE(asec(S::Zero()) == S::ComplexInfinity());     // asec(0) = zoo
    auto x = symbol("x");
    REQUIRE(asec(x)->type_id() == TypeId::Function);
}

TEST_CASE("acsc: canonical values and odd parity", "[3e][acsc]") {
    auto pi = S::Pi();
    REQUIRE(acsc(S::One()) == mul(S::Half(), pi));        // acsc(1) = pi/2
    REQUIRE(acsc(integer(2)) == mul(rational(1, 6), pi));  // acsc(2) = pi/6
    REQUIRE(acsc(S::NegativeOne()) == mul(rational(-1, 2), pi));  // odd
    REQUIRE(acsc(S::Zero()) == S::ComplexInfinity());     // acsc(0) = zoo
    auto x = symbol("x");
    REQUIRE(acsc(x)->type_id() == TypeId::Function);
}

TEST_CASE("acot/asec/acsc: derivatives match SymPy", "[3e][acot][asec][acsc][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(diff(acot(x), x)->str(), "-1/(1 + x**2)"));
    REQUIRE(oracle.equivalent(diff(asec(x), x)->str(), "1/(x**2*sqrt(1 - 1/x**2))"));
    REQUIRE(oracle.equivalent(diff(acsc(x), x)->str(), "-1/(x**2*sqrt(1 - 1/x**2))"));
}

TEST_CASE("acot/asec/acsc: parse round-trip", "[3e][acot][asec][acsc][parser]") {
    auto x = symbol("x");
    REQUIRE(parsing::parse("acot(x)") == acot(x));
    REQUIRE(parsing::parse("asec(x)") == asec(x));
    REQUIRE(parsing::parse("acsc(x)") == acsc(x));
    REQUIRE(acot(x)->str() == "acot(x)");
    REQUIRE(asec(x)->str() == "asec(x)");
    REQUIRE(acsc(x)->str() == "acsc(x)");
}
