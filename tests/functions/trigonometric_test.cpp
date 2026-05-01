// sin / cos / tan tests.
//
// Reference: sympy/functions/elementary/tests/test_trigonometric.py
//   * test_sin   (line 52)
//   * test_cos   (line 309)
//   * test_tan
//
// Phase 3b ships the canonical-angle subset; multi-period reductions
// (sin(2*pi) = 0, sin(7*pi) = 0) come in a follow-up.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/functions/trigonometric.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Auto-eval -------------------------------------------------------------

TEST_CASE("sin: canonical angles", "[3b][trig][sin]") {
    REQUIRE(sin(S::Zero()) == S::Zero());
    REQUIRE(sin(S::Pi()) == S::Zero());
    REQUIRE(sin(mul(S::Half(), S::Pi())) == S::One());
}

TEST_CASE("cos: canonical angles", "[3b][trig][cos]") {
    REQUIRE(cos(S::Zero()) == S::One());
    REQUIRE(cos(S::Pi()) == S::NegativeOne());
    REQUIRE(cos(mul(S::Half(), S::Pi())) == S::Zero());
}

TEST_CASE("tan: canonical angles", "[3b][trig][tan]") {
    REQUIRE(tan(S::Zero()) == S::Zero());
    REQUIRE(tan(S::Pi()) == S::Zero());
}

TEST_CASE("sin: odd identity sin(-x) = -sin(x)", "[3b][trig][sin]") {
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    auto e = sin(neg);
    auto expected = mul(S::NegativeOne(), sin(x));
    REQUIRE(e == expected);
}

TEST_CASE("cos: even identity cos(-x) = cos(x)", "[3b][trig][cos]") {
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    REQUIRE(cos(neg) == cos(x));
}

TEST_CASE("tan: odd identity tan(-x) = -tan(x)", "[3b][trig][tan]") {
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    auto e = tan(neg);
    auto expected = mul(S::NegativeOne(), tan(x));
    REQUIRE(e == expected);
}

TEST_CASE("sin: stays unevaluated on a generic symbol", "[3b][trig][sin]") {
    auto x = symbol("x");
    auto e = sin(x);
    REQUIRE(e->type_id() == TypeId::Function);
    REQUIRE(e->str() == "sin(x)");
}

// ----- ask propagation -------------------------------------------------------

TEST_CASE("sin: real for real arg", "[3b][trig][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(is_real(sin(x)) == true);
    REQUIRE(is_finite(sin(x)) == true);
}

TEST_CASE("cos: real for real arg", "[3b][trig][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(is_real(cos(x)) == true);
    REQUIRE(is_finite(cos(x)) == true);
}

// ----- Numeric evalf ---------------------------------------------------------

TEST_CASE("sin: numeric Float arg evalf via mpfr_sin", "[3b][trig][evalf]") {
    auto half = float_value(0.5);
    auto e = sin(half);
    REQUIRE(e->type_id() == TypeId::Float);
}

// ----- Substitution roundtrip -----------------------------------------------

TEST_CASE("sin/cos: subs replaces inside the argument", "[3b][trig][subs]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = sin(x) + cos(x);
    auto r = subs(e, x, y);
    REQUIRE(r == sin(y) + cos(y));
}

TEST_CASE("sin: substituting x = 0 yields 0", "[3b][trig][subs]") {
    auto x = symbol("x");
    REQUIRE(subs(sin(x), x, S::Zero()) == S::Zero());
}

TEST_CASE("cos: substituting x = pi yields -1", "[3b][trig][subs]") {
    auto x = symbol("x");
    REQUIRE(subs(cos(x), x, S::Pi()) == S::NegativeOne());
}

// ----- Oracle-validated tests ------------------------------------------------

TEST_CASE("sin/cos/tan: structural output matches SymPy",
          "[3b][trig][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(sin(x)->str(), "sin(x)"));
    REQUIRE(oracle.equivalent(cos(x)->str(), "cos(x)"));
    REQUIRE(oracle.equivalent(tan(x)->str(), "tan(x)"));

    REQUIRE(oracle.equivalent(sin(S::Zero())->str(), "0"));
    REQUIRE(oracle.equivalent(cos(S::Zero())->str(), "1"));
    REQUIRE(oracle.equivalent(sin(S::Pi())->str(), "0"));
    REQUIRE(oracle.equivalent(cos(S::Pi())->str(), "-1"));
}

TEST_CASE("sin: numeric evalf matches SymPy", "[3b][trig][evalf][oracle]") {
    // The SymPy reference must use a Float constructed at the target precision
    // (or a Rational), otherwise sympify("sin(0.5)") locks in 53-bit precision
    // before evalf and the answer is just 53 bits printed at 30 dps with noise.
    auto& oracle = Oracle::instance();
    auto e = sin(float_value("0.5", 30));
    auto sympy_value = oracle.evalf("sin(Rational(1, 2))", 30);
    INFO("SymPP: " << e->str());
    INFO("SymPy: " << sympy_value);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

TEST_CASE("cos: numeric evalf matches SymPy", "[3b][trig][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    // 1.234567 = 1234567 / 1000000
    auto e = cos(float_value("1.234567", 30));
    auto sympy_value = oracle.evalf("cos(Rational(1234567, 1000000))", 30);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

TEST_CASE("trig: composite expressions match SymPy", "[3b][trig][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sin(x) + integer(2) * cos(x);
    REQUIRE(oracle.equivalent(e->str(), "sin(x) + 2*cos(x)"));
}
