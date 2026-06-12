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
#include <sympp/calculus/diff.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/parsing/parser.hpp>

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

// ----- π/12 family exact values (regression, TRIG-PI12) ----------------------
// The 15° angles fold to clean sums of radicals (not nested): cos(π/12) =
// (√6+√2)/4, sin(π/12) = (√6−√2)/4, tan(π/12) = 2−√3. All 24 multiples of π/12
// reduce through cos_pi's period/reflection folds. π/8 (genuinely nested) stays
// symbolic.
TEST_CASE("sin/cos/tan: π/12 family exact values", "[3b][trig][regression][oracle]") {
    auto& oracle = Oracle::instance();
    auto pi = S::Pi();
    REQUIRE(oracle.equivalent(cos(mul(rational(1, 12), pi))->str(),
                              "sqrt(6)/4 + sqrt(2)/4"));
    REQUIRE(oracle.equivalent(sin(mul(rational(1, 12), pi))->str(),
                              "sqrt(6)/4 - sqrt(2)/4"));
    REQUIRE(oracle.equivalent(tan(mul(rational(1, 12), pi))->str(), "2 - sqrt(3)"));
    REQUIRE(oracle.equivalent(cos(mul(rational(5, 12), pi))->str(),
                              "sqrt(6)/4 - sqrt(2)/4"));
    REQUIRE(oracle.equivalent(tan(mul(rational(5, 12), pi))->str(), "2 + sqrt(3)"));
    // A higher multiple reduces by symmetry.
    REQUIRE(oracle.equivalent(cos(mul(rational(7, 12), pi))->str(),
                              "-sqrt(6)/4 + sqrt(2)/4"));
}

TEST_CASE("cos: π/8 (genuinely nested radical) stays unevaluated",
          "[3b][trig][regression]") {
    // cos(π/8) = √(2+√2)/2 — a nested radical, deliberately not in the table.
    auto r = cos(mul(rational(1, 8), S::Pi()));
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

// ----- Half-integer π co-function shifts (regression, TRIG-4) ----------------
// sin(x+(m/2)π) = ±cos(x), cos(x+(m/2)π) = ∓sin(x) for odd m.
TEST_CASE("sin/cos: half-integer π co-function shift", "[3b][trig][regression]") {
    auto x = symbol("x");
    auto pi = S::Pi();
    REQUIRE(sin(x + mul(rational(1, 2), pi)) == cos(x));
    REQUIRE(sin(x + mul(rational(-1, 2), pi)) == mul(S::NegativeOne(), cos(x)));
    REQUIRE(sin(x + mul(rational(3, 2), pi)) == mul(S::NegativeOne(), cos(x)));
    REQUIRE(cos(x + mul(rational(1, 2), pi)) == mul(S::NegativeOne(), sin(x)));
    REQUIRE(cos(x + mul(rational(-1, 2), pi)) == sin(x));
    REQUIRE(cos(x + mul(rational(3, 2), pi)) == sin(x));
}

TEST_CASE("tan: half-integer π shift stays symbolic (−cot, no cot type)",
          "[3b][trig][tan][regression]") {
    // tan(x+π/2) = −cot(x); SymPP has no cot, so it is left unevaluated.
    auto x = symbol("x");
    auto r = tan(x + mul(rational(1, 2), S::Pi()));
    REQUIRE(r->type_id() == TypeId::Function);
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

// ----- Reciprocal trio: cot, sec, csc ----------------------------------------

TEST_CASE("cot/sec/csc: canonical-angle exact values", "[3b][trig][reciprocal]") {
    auto pi = S::Pi();
    // cot(π/4) = 1, cot(3π/4) = -1, cot(π/2) = 0.
    REQUIRE(cot(mul(rational(1, 4), pi)) == S::One());
    REQUIRE(cot(mul(rational(3, 4), pi)) == S::NegativeOne());
    REQUIRE(cot(mul(S::Half(), pi)) == S::Zero());
    REQUIRE(cot(mul(rational(1, 6), pi)) == pow(integer(3), rational(1, 2)));
    // sec(0) = 1, sec(π/3) = 2, sec(2π/3) = -2.
    REQUIRE(sec(S::Zero()) == S::One());
    REQUIRE(sec(mul(rational(1, 3), pi)) == integer(2));
    REQUIRE(sec(mul(rational(2, 3), pi)) == integer(-2));
    REQUIRE(sec(S::Pi()) == S::NegativeOne());
    // csc(π/2) = 1, csc(π/6) = 2.
    REQUIRE(csc(mul(S::Half(), pi)) == S::One());
    REQUIRE(csc(mul(rational(1, 6), pi)) == integer(2));
}

TEST_CASE("cot/sec/csc: poles fold to ComplexInfinity (zoo)",
          "[3b][trig][reciprocal]") {
    auto pi = S::Pi();
    REQUIRE(cot(S::Zero()) == S::ComplexInfinity());  // sin = 0
    REQUIRE(cot(pi) == S::ComplexInfinity());
    REQUIRE(csc(S::Zero()) == S::ComplexInfinity());
    REQUIRE(sec(mul(S::Half(), pi)) == S::ComplexInfinity());  // cos = 0
}

TEST_CASE("cot/sec/csc: parity (cot/csc odd, sec even)",
          "[3b][trig][reciprocal]") {
    auto x = symbol("x");
    REQUIRE(cot(mul(S::NegativeOne(), x)) == mul(S::NegativeOne(), cot(x)));
    REQUIRE(csc(mul(S::NegativeOne(), x)) == mul(S::NegativeOne(), csc(x)));
    REQUIRE(sec(mul(S::NegativeOne(), x)) == sec(x));
}

TEST_CASE("cot/sec/csc: inverse-function compositions", "[3b][trig][reciprocal]") {
    auto x = symbol("x");
    REQUIRE(cot(atan(x)) == pow(x, S::NegativeOne()));   // cot(atan x) = 1/x
    REQUIRE(sec(acos(x)) == pow(x, S::NegativeOne()));   // sec(acos x) = 1/x
    REQUIRE(csc(asin(x)) == pow(x, S::NegativeOne()));   // csc(asin x) = 1/x
}

TEST_CASE("cot/sec/csc: real for real argument", "[3b][trig][reciprocal][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(is_real(cot(x)) == true);
    REQUIRE(is_real(sec(x)) == true);
    REQUIRE(is_real(csc(x)) == true);
}

TEST_CASE("cot/sec/csc: parse round-trip", "[3b][trig][reciprocal][parser]") {
    auto x = symbol("x");
    REQUIRE(parsing::parse("cot(x)") == cot(x));
    REQUIRE(parsing::parse("sec(x)") == sec(x));
    REQUIRE(parsing::parse("csc(x)") == csc(x));
    // str() emits the canonical spelling, so parse(e->str()) == e.
    REQUIRE(cot(x)->str() == "cot(x)");
    REQUIRE(sec(x)->str() == "sec(x)");
    REQUIRE(csc(x)->str() == "csc(x)");
}

TEST_CASE("cot/sec/csc: derivatives match SymPy", "[3b][trig][reciprocal][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // d/dx cot = -csc², d/dx sec = sec·tan, d/dx csc = -csc·cot.
    REQUIRE(oracle.equivalent(diff(cot(x), x)->str(), "-csc(x)**2"));
    REQUIRE(oracle.equivalent(diff(sec(x), x)->str(), "sec(x)*tan(x)"));
    REQUIRE(oracle.equivalent(diff(csc(x), x)->str(), "-csc(x)*cot(x)"));
}

TEST_CASE("cot/sec/csc: exact values match SymPy (radical forms)",
          "[3b][trig][reciprocal][oracle]") {
    auto& oracle = Oracle::instance();
    auto pi = S::Pi();
    REQUIRE(oracle.equivalent(cot(mul(rational(1, 3), pi))->str(), "sqrt(3)/3"));
    REQUIRE(oracle.equivalent(sec(mul(rational(1, 4), pi))->str(), "sqrt(2)"));
    REQUIRE(oracle.equivalent(sec(mul(rational(1, 6), pi))->str(), "2*sqrt(3)/3"));
    REQUIRE(oracle.equivalent(csc(mul(rational(1, 4), pi))->str(), "sqrt(2)"));
}

TEST_CASE("cot/sec/csc: numeric Float arg evalf matches SymPy",
          "[3b][trig][reciprocal][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = cot(float_value("1.5", 30));
    REQUIRE(e->type_id() == TypeId::Float);
    auto sympy_value = oracle.evalf("cot(Rational(3, 2))", 30);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

TEST_CASE("sin/cos/tan at an integer multiple of pi (ASSUME-7)",
          "[3b][trig][assumptions][regression]") {
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    auto x = symbol("x");  // generic
    // sin(kπ)=0, tan(kπ)=0, cos(kπ)=(-1)^k for an integer k.
    REQUIRE(sin(n * S::Pi()) == S::Zero());
    REQUIRE(tan(n * S::Pi()) == S::Zero());
    REQUIRE(cos(n * S::Pi()) == pow(S::NegativeOne(), n));
    REQUIRE(sin(integer(2) * n * S::Pi()) == S::Zero());          // 2n integer
    REQUIRE(cos((n + integer(1)) * S::Pi())
            == pow(S::NegativeOne(), n + integer(1)));
    // A non-integer (generic) coefficient stays unevaluated.
    REQUIRE(sin(x * S::Pi())->str() == "sin(x*pi)");
}

TEST_CASE("sin/cos at an odd half-integer multiple of pi (ASSUME-9)",
          "[3b][trig][assumptions][regression]") {
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    // (2n+1)·π/2 is an odd half-integer multiple: cos = 0, sin = (-1)^n.
    Expr arg = (integer(2) * n + integer(1)) * S::Pi() / integer(2);
    REQUIRE(cos(arg) == S::Zero());
    REQUIRE(sin(arg) == pow(S::NegativeOne(), n));
    // Literal multiples still resolve precisely.
    REQUIRE(cos(S::Pi() / integer(2)) == S::Zero());
    REQUIRE(sin(integer(3) * S::Pi() / integer(2)) == S::NegativeOne());
}

TEST_CASE("cot/sec/csc at integer / half-integer multiples of pi (ASSUME-10)",
          "[3b][trig][reciprocal][assumptions][regression]") {
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    Expr half = (integer(2) * n + integer(1)) * S::Pi() / integer(2);
    // Integer k·π: cot/csc are poles (sin=0); sec = (-1)^k (1/cos).
    REQUIRE(cot(n * S::Pi()) == S::ComplexInfinity());
    REQUIRE(csc(n * S::Pi()) == S::ComplexInfinity());
    REQUIRE(sec(n * S::Pi()) == pow(S::NegativeOne(), n));
    // Odd half-integer: sec is a pole (cos=0); cot=0; csc = (-1)^n (1/sin).
    REQUIRE(sec(half) == S::ComplexInfinity());
    REQUIRE(cot(half) == S::Zero());
    REQUIRE(csc(half) == pow(S::NegativeOne(), n));
}
