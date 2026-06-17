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
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/calculus/diff.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/parsing/parser.hpp>

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

// Regression (FUNC-INF): hyperbolic functions at ±oo.
TEST_CASE("hyperbolic: values at infinity", "[3f][infinity][regression]") {
    auto oo = S::Infinity();
    auto noo = S::NegativeInfinity();
    REQUIRE(sinh(oo) == oo);
    REQUIRE(sinh(noo) == noo);
    REQUIRE(cosh(oo) == oo);
    REQUIRE(cosh(noo) == oo);
    REQUIRE(tanh(oo) == S::One());
    REQUIRE(tanh(noo) == S::NegativeOne());
    REQUIRE(asinh(oo) == oo);
    REQUIRE(asinh(noo) == noo);
    REQUIRE(acosh(oo) == oo);
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

// ----- f(f⁻¹(x)) = x compositions (regression, FUNC-1) -----------------------
TEST_CASE("hyperbolic: f(inverse(x)) collapses to x", "[3f][regression]") {
    auto x = symbol("x");
    REQUIRE(sinh(asinh(x)) == x);
    REQUIRE(cosh(acosh(x)) == x);
    REQUIRE(tanh(atanh(x)) == x);
    // Reverse direction stays unevaluated.
    REQUIRE(asinh(sinh(x))->type_id() == TypeId::Function);
}

// ----- Cross-function inverse compositions (regression, FUNC-3) --------------
// cosh(asinh(x)) = √(x²+1), etc. — the hyperbolic analogue of FUNC-2.
TEST_CASE("hyperbolic: cross-function inverse compositions",
          "[3f][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(cosh(asinh(x))->str(), "sqrt(x**2 + 1)"));
    REQUIRE(oracle.equivalent(tanh(asinh(x))->str(), "x/sqrt(x**2 + 1)"));
    REQUIRE(oracle.equivalent(cosh(atanh(x))->str(), "1/sqrt(1 - x**2)"));
    REQUIRE(oracle.equivalent(sinh(atanh(x))->str(), "x/sqrt(1 - x**2)"));
    REQUIRE(oracle.equivalent(sinh(acosh(x))->str(), "sqrt(x - 1)*sqrt(x + 1)"));
    REQUIRE(oracle.equivalent(tanh(acosh(x))->str(),
                              "sqrt(x - 1)*sqrt(x + 1)/x"));
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

// ACOSH-IMAG-1: acosh(x) = i·acos(x) for a real x ∈ [−1, 1] whose acos has a
// closed form — acosh(0)=iπ/2, acosh(½)=iπ/3, acosh(−1)=iπ. A rational with no
// nice acos value (acosh(⅓)) or |x|>1 (acosh(2)) is left symbolic, as SymPy does.
TEST_CASE("acosh: imaginary values on [−1, 1] (ACOSH-IMAG-1)",
          "[3f][acosh][oracle][regression]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(acosh(S::Zero())->str(), "I*pi/2"));
    REQUIRE(oracle.equivalent(acosh(rational(1, 2))->str(), "I*pi/3"));
    REQUIRE(oracle.equivalent(acosh(rational(-1, 2))->str(), "2*I*pi/3"));
    REQUIRE(oracle.equivalent(acosh(integer(-1))->str(), "I*pi"));
    // No over-reach: no closed form, or outside [−1, 1], stays an Acosh node.
    REQUIRE(acosh(rational(1, 3))->type_id() == TypeId::Function);
    REQUIRE(acosh(integer(2))->type_id() == TypeId::Function);
    REQUIRE(acosh(integer(-2))->type_id() == TypeId::Function);
    // Inverse composition still collapses.
    auto x = symbol("x");
    REQUIRE(cosh(acosh(x)) == x);
}

TEST_CASE("atanh: canonical and odd", "[3f][atanh]") {
    REQUIRE(atanh(S::Zero()) == S::Zero());
    auto x = symbol("x");
    auto neg = mul(S::NegativeOne(), x);
    REQUIRE(atanh(neg) == mul(S::NegativeOne(), atanh(x)));
}

// ATANH-POLE-1: atanh(±1) = ±∞ (the real pole, since ½·log((1+x)/(1−x)) → ±∞),
// and acoth(±1) = ±∞ follows via acoth = atanh of the reciprocal. Matches SymPy.
TEST_CASE("atanh/acoth: poles at ±1 (ATANH-POLE-1)", "[3f][atanh][acoth][regression]") {
    REQUIRE(atanh(S::One()) == S::Infinity());
    REQUIRE(atanh(integer(-1)) == S::NegativeInfinity());
    REQUIRE(acoth(S::One()) == S::Infinity());
    REQUIRE(acoth(integer(-1)) == S::NegativeInfinity());
    // No over-reach: interior and exterior arguments stay symbolic.
    REQUIRE(atanh(rational(1, 2))->type_id() == TypeId::Function);
    REQUIRE(atanh(integer(2))->type_id() == TypeId::Function);
    REQUIRE(acoth(integer(2))->type_id() == TypeId::Function);
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

// ----- Reciprocal trio: coth, sech, csch -------------------------------------

TEST_CASE("coth/sech/csch: canonical values and poles", "[3f][reciprocal]") {
    // sech(0) = 1; coth(0) = csch(0) = zoo (poles).
    REQUIRE(sech(S::Zero()) == S::One());
    REQUIRE(coth(S::Zero()) == S::ComplexInfinity());
    REQUIRE(csch(S::Zero()) == S::ComplexInfinity());
    // At ±oo: coth → ±1, sech → 0, csch → 0.
    REQUIRE(coth(S::Infinity()) == S::One());
    REQUIRE(coth(S::NegativeInfinity()) == S::NegativeOne());
    REQUIRE(sech(S::Infinity()) == S::Zero());
    REQUIRE(csch(S::Infinity()) == S::Zero());
}

TEST_CASE("coth/sech/csch: parity (coth/csch odd, sech even)", "[3f][reciprocal]") {
    auto x = symbol("x");
    REQUIRE(coth(mul(S::NegativeOne(), x)) == mul(S::NegativeOne(), coth(x)));
    REQUIRE(csch(mul(S::NegativeOne(), x)) == mul(S::NegativeOne(), csch(x)));
    REQUIRE(sech(mul(S::NegativeOne(), x)) == sech(x));
}

TEST_CASE("coth/sech/csch: inverse-function compositions", "[3f][reciprocal]") {
    auto x = symbol("x");
    REQUIRE(coth(atanh(x)) == pow(x, integer(-1)));  // coth(atanh x) = 1/x
    REQUIRE(sech(acosh(x)) == pow(x, integer(-1)));  // sech(acosh x) = 1/x
    REQUIRE(csch(asinh(x)) == pow(x, integer(-1)));  // csch(asinh x) = 1/x
}

TEST_CASE("coth/sech/csch: real for real argument",
          "[3f][reciprocal][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(is_real(coth(x)) == true);
    REQUIRE(is_real(sech(x)) == true);
    REQUIRE(is_real(csch(x)) == true);
}

TEST_CASE("coth/sech/csch: parse round-trip", "[3f][reciprocal][parser]") {
    auto x = symbol("x");
    REQUIRE(parsing::parse("coth(x)") == coth(x));
    REQUIRE(parsing::parse("sech(x)") == sech(x));
    REQUIRE(parsing::parse("csch(x)") == csch(x));
    REQUIRE(coth(x)->str() == "coth(x)");
    REQUIRE(sech(x)->str() == "sech(x)");
    REQUIRE(csch(x)->str() == "csch(x)");
}

TEST_CASE("coth/sech/csch: derivatives match SymPy", "[3f][reciprocal][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // coth' = -csch², sech' = -sech·tanh, csch' = -csch·coth.
    REQUIRE(oracle.equivalent(diff(coth(x), x)->str(), "-csch(x)**2"));
    REQUIRE(oracle.equivalent(diff(sech(x), x)->str(), "-sech(x)*tanh(x)"));
    REQUIRE(oracle.equivalent(diff(csch(x), x)->str(), "-csch(x)*coth(x)"));
}

TEST_CASE("coth/sech/csch: numeric Float arg evalf matches SymPy",
          "[3f][reciprocal][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = coth(float_value("1.5", 30));
    REQUIRE(e->type_id() == TypeId::Float);
    REQUIRE(oracle.equivalent(e->str(), oracle.evalf("coth(Rational(3, 2))", 30)));
}

// ----- Inverse reciprocal trio: acoth, asech, acsch (HYP-AINV-RECIP) ---------

TEST_CASE("acoth/asech/acsch: canonical values and poles", "[3f][reciprocal]") {
    // acoth(±oo) = 0; asech(0) = oo, asech(1) = 0; acsch(0) = zoo, acsch(oo) = 0.
    REQUIRE(acoth(S::Infinity()) == S::Zero());
    REQUIRE(acoth(S::NegativeInfinity()) == S::Zero());
    REQUIRE(asech(S::Zero()) == S::Infinity());
    REQUIRE(asech(S::One()) == S::Zero());
    REQUIRE(acsch(S::Zero()) == S::ComplexInfinity());
    REQUIRE(acsch(S::Infinity()) == S::Zero());
    // Non-special integer args stay unevaluated (SymPy keeps them too / gives
    // complex or log forms SymPP doesn't model).
    auto x = symbol("x");
    REQUIRE(acoth(integer(2))->type_id() == TypeId::Function);
    REQUIRE(acoth(x)->type_id() == TypeId::Function);
    REQUIRE(asech(x)->type_id() == TypeId::Function);
    REQUIRE(acsch(x)->type_id() == TypeId::Function);
}

TEST_CASE("acoth/acsch: odd parity (asech has none)", "[3f][reciprocal]") {
    auto x = symbol("x");
    REQUIRE(acoth(mul(S::NegativeOne(), x)) == mul(S::NegativeOne(), acoth(x)));
    REQUIRE(acsch(mul(S::NegativeOne(), x)) == mul(S::NegativeOne(), acsch(x)));
}

TEST_CASE("acoth/asech/acsch: derivatives match SymPy",
          "[3f][reciprocal][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(diff(acoth(x), x)->str(), "1/(1 - x**2)"));
    REQUIRE(oracle.equivalent(diff(asech(x), x)->str(), "-1/(x*sqrt(1 - x**2))"));
    REQUIRE(oracle.equivalent(diff(acsch(x), x)->str(),
                              "-1/(x**2*sqrt(1 + 1/x**2))"));
}

TEST_CASE("acoth/asech/acsch: parse round-trip", "[3f][reciprocal][parser]") {
    auto x = symbol("x");
    REQUIRE(parsing::parse("acoth(x)") == acoth(x));
    REQUIRE(parsing::parse("asech(x)") == asech(x));
    REQUIRE(parsing::parse("acsch(x)") == acsch(x));
    REQUIRE(acoth(x)->str() == "acoth(x)");
    REQUIRE(asech(x)->str() == "asech(x)");
    REQUIRE(acsch(x)->str() == "acsch(x)");
}
