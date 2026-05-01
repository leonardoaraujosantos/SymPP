// sinh / cosh / tanh / asinh / acosh / atanh tests.
//
// Reference: sympy/functions/elementary/tests/test_hyperbolic.py
//   * test_sinh, test_cosh, test_tanh
//   * test_asinh, test_acosh, test_atanh

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/functions/hyperbolic.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- sinh / cosh / tanh ----------------------------------------------------

TEST_CASE("sinh: canonical values", "[3f][sinh]") {
    REQUIRE(sinh(S::Zero()) == S::Zero());
}

TEST_CASE("cosh: canonical values", "[3f][cosh]") {
    REQUIRE(cosh(S::Zero()) == S::One());
}

TEST_CASE("tanh: canonical values", "[3f][tanh]") {
    REQUIRE(tanh(S::Zero()) == S::Zero());
}

TEST_CASE("sinh: odd identity sinh(-x) = -sinh(x)", "[3f][sinh]") {
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    REQUIRE(sinh(neg) == mul(S::NegativeOne(), sinh(x)));
}

TEST_CASE("cosh: even identity cosh(-x) = cosh(x)", "[3f][cosh]") {
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    REQUIRE(cosh(neg) == cosh(x));
}

TEST_CASE("tanh: odd identity tanh(-x) = -tanh(x)", "[3f][tanh]") {
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    REQUIRE(tanh(neg) == mul(S::NegativeOne(), tanh(x)));
}

// ----- Assumptions -----------------------------------------------------------

TEST_CASE("sinh/tanh: real for real arg", "[3f][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(is_real(sinh(x)) == true);
    REQUIRE(is_real(tanh(x)) == true);
    REQUIRE(is_finite(sinh(x)) == true);
}

TEST_CASE("cosh: real arg yields positive cosh", "[3f][cosh][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(is_real(cosh(x)) == true);
    REQUIRE(is_positive(cosh(x)) == true);
    REQUIRE(is_nonzero(cosh(x)) == true);
}

// ----- Inverse hyperbolic ----------------------------------------------------

TEST_CASE("asinh: canonical and odd", "[3f][asinh]") {
    REQUIRE(asinh(S::Zero()) == S::Zero());
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    REQUIRE(asinh(neg) == mul(S::NegativeOne(), asinh(x)));
}

TEST_CASE("acosh(1) = 0", "[3f][acosh]") {
    REQUIRE(acosh(S::One()) == S::Zero());
}

TEST_CASE("atanh: canonical and odd", "[3f][atanh]") {
    REQUIRE(atanh(S::Zero()) == S::Zero());
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    REQUIRE(atanh(neg) == mul(S::NegativeOne(), atanh(x)));
}

// ----- Numeric evalf ---------------------------------------------------------

TEST_CASE("hyperbolic: numeric Float arg evaluates", "[3f][evalf]") {
    REQUIRE(sinh(float_value(0.5))->type_id() == TypeId::Float);
    REQUIRE(cosh(float_value(0.5))->type_id() == TypeId::Float);
    REQUIRE(tanh(float_value(0.5))->type_id() == TypeId::Float);
    REQUIRE(asinh(float_value(0.5))->type_id() == TypeId::Float);
    REQUIRE(acosh(float_value(2.0))->type_id() == TypeId::Float);
    REQUIRE(atanh(float_value(0.5))->type_id() == TypeId::Float);
}

TEST_CASE("sinh: high-precision matches SymPy",
          "[3f][sinh][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = sinh(float_value("1", 30));
    auto sympy_value = oracle.evalf("sinh(Integer(1))", 30);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

TEST_CASE("cosh: high-precision matches SymPy",
          "[3f][cosh][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = cosh(float_value("1", 30));
    auto sympy_value = oracle.evalf("cosh(Integer(1))", 30);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

TEST_CASE("tanh: high-precision matches SymPy",
          "[3f][tanh][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = tanh(float_value("0.5", 30));
    auto sympy_value = oracle.evalf("tanh(Rational(1, 2))", 30);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

// ----- Substitution ----------------------------------------------------------

TEST_CASE("hyperbolic: subs propagates", "[3f][subs]") {
    auto x = symbol("x");
    auto e = sinh(x) + cosh(x) + tanh(x);
    auto r = subs(e, x, S::Zero());
    // sinh(0)+cosh(0)+tanh(0) = 0+1+0 = 1
    REQUIRE(r == S::One());
}

// ----- Oracle structural -----------------------------------------------------

TEST_CASE("hyperbolic: structural output matches SymPy",
          "[3f][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(sinh(x)->str(), "sinh(x)"));
    REQUIRE(oracle.equivalent(cosh(x)->str(), "cosh(x)"));
    REQUIRE(oracle.equivalent(tanh(x)->str(), "tanh(x)"));
    REQUIRE(oracle.equivalent(asinh(x)->str(), "asinh(x)"));
    REQUIRE(oracle.equivalent(acosh(x)->str(), "acosh(x)"));
    REQUIRE(oracle.equivalent(atanh(x)->str(), "atanh(x)"));
}
