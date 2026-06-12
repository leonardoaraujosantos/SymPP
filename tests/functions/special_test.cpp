// erf / erfc / Heaviside / DiracDelta tests.
//
// References:
//   sympy/functions/special/tests/test_error_functions.py
//   sympy/functions/special/tests/test_delta_functions.py

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
#include <sympp/functions/special.hpp>
#include <sympp/parsing/parser.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- erf -------------------------------------------------------------------

TEST_CASE("erf: erf(0) = 0", "[3j][erf]") {
    REQUIRE(erf(S::Zero()) == S::Zero());
}

TEST_CASE("erf: odd identity erf(-x) = -erf(x)", "[3j][erf]") {
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    REQUIRE(erf(neg) == mul(S::NegativeOne(), erf(x)));
}

TEST_CASE("erf: numeric Float matches mpfr_erf", "[3j][erf]") {
    auto e = erf(float_value(1.0));
    REQUIRE(e->type_id() == TypeId::Float);
}

TEST_CASE("erf: high-precision matches SymPy",
          "[3j][erf][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = erf(float_value("1", 30));
    auto sympy_value = oracle.evalf("erf(Integer(1))", 30);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

TEST_CASE("erf: real arg propagates", "[3j][erf][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(is_real(erf(x)) == true);
    REQUIRE(is_finite(erf(x)) == true);
}

// ----- erfc ------------------------------------------------------------------

TEST_CASE("erfc: erfc(0) = 1", "[3j][erfc]") {
    REQUIRE(erfc(S::Zero()) == S::One());
}

// Regression (FUNC-INF): erf/erfc at ±oo.
TEST_CASE("erf/erfc: values at infinity", "[3j][erf][erfc][infinity][regression]") {
    REQUIRE(erf(S::Infinity()) == S::One());
    REQUIRE(erf(S::NegativeInfinity()) == S::NegativeOne());
    REQUIRE(erfc(S::Infinity()) == S::Zero());
    REQUIRE(erfc(S::NegativeInfinity()) == integer(2));
}

TEST_CASE("erfc: numeric matches mpfr_erfc", "[3j][erfc]") {
    auto e = erfc(float_value(1.0));
    REQUIRE(e->type_id() == TypeId::Float);
}

// ----- Heaviside -------------------------------------------------------------

TEST_CASE("Heaviside(0) = 1/2 (SymPy default)", "[3j][heaviside]") {
    REQUIRE(heaviside(S::Zero()) == S::Half());
}

TEST_CASE("Heaviside(positive) = 1", "[3j][heaviside]") {
    auto p = symbol("p", AssumptionMask{}.set_positive(true));
    REQUIRE(heaviside(p) == S::One());
    REQUIRE(heaviside(integer(5)) == S::One());
}

TEST_CASE("Heaviside(negative) = 0", "[3j][heaviside]") {
    auto n = symbol("n", AssumptionMask{}.set_negative(true));
    REQUIRE(heaviside(n) == S::Zero());
    REQUIRE(heaviside(integer(-3)) == S::Zero());
}

TEST_CASE("Heaviside: stays unevaluated for unknown sign",
          "[3j][heaviside]") {
    auto x = symbol("x");
    auto e = heaviside(x);
    REQUIRE(e->type_id() == TypeId::Function);
    REQUIRE(e->str() == "Heaviside(x)");
}

// ----- DiracDelta ------------------------------------------------------------

TEST_CASE("DiracDelta(nonzero) = 0", "[3j][dirac]") {
    auto p = symbol("p", AssumptionMask{}.set_positive(true));
    REQUIRE(dirac_delta(p) == S::Zero());
    REQUIRE(dirac_delta(integer(5)) == S::Zero());
    REQUIRE(dirac_delta(integer(-5)) == S::Zero());
}

TEST_CASE("DiracDelta(0) stays unevaluated", "[3j][dirac]") {
    auto e = dirac_delta(S::Zero());
    REQUIRE(e->type_id() == TypeId::Function);
    REQUIRE(e->str() == "DiracDelta(0)");
}

TEST_CASE("DiracDelta(generic) stays unevaluated", "[3j][dirac]") {
    auto x = symbol("x");
    auto e = dirac_delta(x);
    REQUIRE(e->str() == "DiracDelta(x)");
}

// ----- Substitution ----------------------------------------------------------

TEST_CASE("erf/erfc/Heaviside: subs propagates", "[3j][subs]") {
    auto x = symbol("x");
    REQUIRE(subs(erf(x), x, S::Zero()) == S::Zero());
    REQUIRE(subs(erfc(x), x, S::Zero()) == S::One());
    REQUIRE(subs(heaviside(x), x, integer(7)) == S::One());
    REQUIRE(subs(heaviside(x), x, integer(-7)) == S::Zero());
}

// ----- Oracle structural -----------------------------------------------------

TEST_CASE("special functions: structural output matches SymPy",
          "[3j][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(erf(x)->str(), "erf(x)"));
    REQUIRE(oracle.equivalent(erfc(x)->str(), "erfc(x)"));
    REQUIRE(oracle.equivalent(heaviside(x)->str(), "Heaviside(x)"));
    REQUIRE(oracle.equivalent(dirac_delta(x)->str(), "DiracDelta(x)"));
}

// ----- Derivatives (regression, DIFF-2) --------------------------------------
// erf/erfc/Heaviside had no diff_arg override, so diff() returned 0 — the same
// class of bug as DIFF-1 (Abs). Now: erf' = 2·exp(−x²)/√π, erfc' = −that,
// Heaviside' = DiracDelta(x).
TEST_CASE("erf/erfc: derivative is the Gaussian", "[3j][erf][diff][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(diff(erf(x), x)->str(), "2*exp(-x**2)/sqrt(pi)"));
    REQUIRE(oracle.equivalent(diff(erfc(x), x)->str(), "-2*exp(-x**2)/sqrt(pi)"));
    // Chain rule still applies.
    REQUIRE(oracle.equivalent(diff(erf(integer(2) * x), x)->str(),
                              "4*exp(-4*x**2)/sqrt(pi)"));
}

TEST_CASE("Heaviside: derivative is DiracDelta", "[3j][diff][regression]") {
    auto x = symbol("x");
    auto d = diff(heaviside(x), x);
    REQUIRE(d->type_id() == TypeId::Function);
    REQUIRE(d->str() == "DiracDelta(x)");
}

// ----- Riemann zeta (ZETA-FN) ------------------------------------------------
// zeta(0)=-1/2, zeta(1)=zoo (pole), even positives zeta(2n)=r*pi^(2n),
// negative integers rational (trivial zeros at even), odd positives symbolic.
TEST_CASE("zeta: special values", "[3h][zeta][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(zeta(integer(2))->str(), "pi**2/6"));
    REQUIRE(oracle.equivalent(zeta(integer(4))->str(), "pi**4/90"));
    REQUIRE(oracle.equivalent(zeta(integer(12))->str(), "691*pi**12/638512875"));
    REQUIRE(zeta(integer(0)) == rational(-1, 2));
    REQUIRE(zeta(integer(1)) == S::ComplexInfinity());     // pole
    REQUIRE(zeta(integer(-1)) == rational(-1, 12));
    REQUIRE(zeta(integer(-2)) == integer(0));               // trivial zero
    REQUIRE(zeta(integer(-3)) == rational(1, 120));
}

TEST_CASE("zeta: odd positive and symbolic stay unevaluated",
          "[3h][zeta][parser]") {
    auto s = symbol("s");
    REQUIRE(zeta(integer(3))->type_id() == TypeId::Function);
    REQUIRE(zeta(s)->type_id() == TypeId::Function);
    REQUIRE(parsing::parse("zeta(s)") == zeta(s));
    REQUIRE(zeta(s)->str() == "zeta(s)");
}
