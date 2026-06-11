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

// ----- Exact values at rational multiples of π (regression, TRIG-1) ----------
// sin/cos/tan at multiples of π/6, π/4, π/3 used to stay unevaluated (only
// 0, π/2, π were special-cased). Now evaluated via a π-coefficient table with
// full quadrant/period reduction. Reference: SymPy's exact-value table.
TEST_CASE("sin/cos: exact rational values (π/6, π/3)",
          "[3b][trig][regression]") {
    auto pi = S::Pi();
    REQUIRE(sin(mul(rational(1, 6), pi)) == S::Half());
    REQUIRE(cos(mul(rational(1, 3), pi)) == S::Half());
    REQUIRE(sin(mul(rational(7, 6), pi)) == mul(S::NegativeOne(), S::Half()));
}

TEST_CASE("tan: exact rational values (π/4, π/3, π/6)",
          "[3b][trig][tan][regression]") {
    auto pi = S::Pi();
    REQUIRE(tan(mul(rational(1, 4), pi)) == S::One());
    REQUIRE(tan(mul(rational(3, 4), pi)) == S::NegativeOne());
    REQUIRE(tan(mul(rational(1, 3), pi)) == pow(integer(3), rational(1, 2)));
}

TEST_CASE("sin/cos: radical exact values match SymPy",
          "[3b][trig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto pi = S::Pi();
    REQUIRE(oracle.equivalent(sin(mul(rational(1, 4), pi))->str(), "sqrt(2)/2"));
    REQUIRE(oracle.equivalent(cos(mul(rational(1, 6), pi))->str(), "sqrt(3)/2"));
    REQUIRE(oracle.equivalent(tan(mul(rational(1, 6), pi))->str(), "sqrt(3)/3"));
    REQUIRE(oracle.equivalent(cos(mul(rational(5, 6), pi))->str(), "-sqrt(3)/2"));
}

TEST_CASE("tan: pole at π/2 stays unevaluated",
          "[3b][trig][tan][regression]") {
    // cos = 0 there → tan is a pole; must not fold to a bogus value.
    auto r = tan(mul(S::Half(), S::Pi()));
    REQUIRE(r->type_id() == TypeId::Function);
}

TEST_CASE("sin: out-of-table angle (π/12) stays unevaluated",
          "[3b][trig][regression]") {
    // Denominator 12 is a nested radical — deliberately not in the table.
    auto r = sin(mul(rational(1, 12), S::Pi()));
    REQUIRE(r->type_id() == TypeId::Function);
}

// ----- Periodicity / π-shift argument reduction (regression, TRIG-3) ---------
// sin(x+kπ)=(−1)^k sin(x), cos(x+kπ)=(−1)^k cos(x), tan(x+kπ)=tan(x) for an
// integer k. Half-integer shifts (the co-function π/2 case) stay symbolic.
TEST_CASE("sin/cos: integer-π argument reduction", "[3b][trig][regression]") {
    auto x = symbol("x");
    auto pi = S::Pi();
    REQUIRE(sin(x + integer(2) * pi) == sin(x));
    REQUIRE(sin(x + pi) == mul(S::NegativeOne(), sin(x)));
    REQUIRE(sin(x + integer(3) * pi) == mul(S::NegativeOne(), sin(x)));
    REQUIRE(cos(x + pi) == mul(S::NegativeOne(), cos(x)));
    REQUIRE(cos(x + integer(4) * pi) == cos(x));
}

TEST_CASE("tan: period-π argument reduction", "[3b][trig][tan][regression]") {
    auto x = symbol("x");
    auto pi = S::Pi();
    REQUIRE(tan(x + pi) == tan(x));
    REQUIRE(tan(x + integer(5) * pi) == tan(x));
}

TEST_CASE("sin: multivariate π-shift reduces", "[3b][trig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto r = sin(x + y + S::Pi());
    REQUIRE(oracle.equivalent(r->str(), "-sin(x + y)"));
}

TEST_CASE("sin: half-integer π shift stays symbolic (co-function, out of scope)",
          "[3b][trig][regression]") {
    auto x = symbol("x");
    auto r = sin(x + mul(rational(1, 2), S::Pi()));
    REQUIRE(r->type_id() == TypeId::Function);  // not folded to cos(x) yet
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
