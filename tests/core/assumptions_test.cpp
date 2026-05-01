// Phase 2a — assumption foundation tests.
//
// Reference inputs from sympy/core/tests/test_assumptions.py:
//   * test_symbol_unset   (line 21)   — Symbol with explicit assumptions
//   * test_zero           (line 30)   — Integer(0) properties
//   * test_one            (line 56)   — Integer(1) properties
//   * test_negativeone    (line 82)   — Integer(-1) properties
//   * test_pi / test_E    (line 274/299) — NumberSymbol properties
//
// Subset: real/positive/negative/zero/nonzero/nonneg/nonpos/integer/rational/finite.
// Even/odd/prime/composite/imaginary/algebraic/transcendental/etc. are deferred.

#include <optional>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number_symbol.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/refine.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

// ----- AssumptionMask closure -----------------------------------------------

TEST_CASE("AssumptionMask: positive => real, finite, nonzero", "[2a][assumptions]") {
    auto m = close_assumptions(AssumptionMask{}.set_positive(true));
    REQUIRE(m.real == true);
    REQUIRE(m.finite == true);
    REQUIRE(m.zero == false);
    REQUIRE(m.negative == false);
}

TEST_CASE("AssumptionMask: integer => rational => real", "[2a][assumptions]") {
    auto m = close_assumptions(AssumptionMask{}.set_integer(true));
    REQUIRE(m.rational == true);
    REQUIRE(m.real == true);
}

TEST_CASE("AssumptionMask: zero implies integer + nonpositive + nonnegative",
          "[2a][assumptions]") {
    auto m = close_assumptions(AssumptionMask{}.set_zero(true));
    REQUIRE(m.real == true);
    REQUIRE(m.integer == true);
    REQUIRE(m.rational == true);
    REQUIRE(m.positive == false);
    REQUIRE(m.negative == false);
    REQUIRE(m.get(AssumptionKey::Nonpositive) == true);
    REQUIRE(m.get(AssumptionKey::Nonnegative) == true);
}

TEST_CASE("AssumptionMask: nonzero is derived from positive | negative",
          "[2a][assumptions]") {
    auto pos = close_assumptions(AssumptionMask{}.set_positive(true));
    REQUIRE(pos.get(AssumptionKey::Nonzero) == true);

    auto neg = close_assumptions(AssumptionMask{}.set_negative(true));
    REQUIRE(neg.get(AssumptionKey::Nonzero) == true);
}

// ----- Symbol with assumptions ----------------------------------------------

TEST_CASE("Symbol: bare symbol has all-unknown assumptions", "[2a][assumptions]") {
    auto x = symbol("x");
    REQUIRE_FALSE(is_real(x).has_value());
    REQUIRE_FALSE(is_positive(x).has_value());
    REQUIRE_FALSE(is_integer(x).has_value());
}

TEST_CASE("Symbol: positive=true propagates to real and nonzero",
          "[2a][assumptions]") {
    // Reference: test_symbol_unset
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    REQUIRE(is_positive(x) == true);
    REQUIRE(is_real(x) == true);
    REQUIRE(is_negative(x) == false);
    REQUIRE(is_zero(x) == false);
    REQUIRE(is_nonzero(x) == true);
    REQUIRE(is_nonnegative(x) == true);
    REQUIRE(is_finite(x) == true);
}

TEST_CASE("Symbol: integer=true propagates to rational and real",
          "[2a][assumptions]") {
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    REQUIRE(is_integer(n) == true);
    REQUIRE(is_rational(n) == true);
    REQUIRE(is_real(n) == true);
}

TEST_CASE("Symbol: real=false propagates to false on rational/integer",
          "[2a][assumptions]") {
    auto z = symbol("z", AssumptionMask{}.set_real(false));
    REQUIRE(is_real(z) == false);
    REQUIRE(is_rational(z) == false);
    REQUIRE(is_integer(z) == false);
    REQUIRE(is_positive(z) == false);
}

TEST_CASE("Symbol: same name, different assumptions are distinct",
          "[2a][assumptions]") {
    auto a = symbol("x");
    auto b = symbol("x", AssumptionMask{}.set_positive(true));
    REQUIRE_FALSE(a == b);
    REQUIRE(a.get() != b.get());  // hash-cons keeps them distinct
}

TEST_CASE("Symbol: same name, same assumptions are interned to one Expr",
          "[2a][assumptions]") {
    auto a = symbol("x", AssumptionMask{}.set_positive(true));
    auto b = symbol("x", AssumptionMask{}.set_positive(true));
    REQUIRE(a.get() == b.get());
}

// ----- Number predicates via ask --------------------------------------------

TEST_CASE("Integer(0): all standard predicates", "[2a][assumptions]") {
    // Reference: test_zero
    auto z = integer(0);
    REQUIRE(is_real(z) == true);
    REQUIRE(is_rational(z) == true);
    REQUIRE(is_integer(z) == true);
    REQUIRE(is_zero(z) == true);
    REQUIRE(is_nonzero(z) == false);
    REQUIRE(is_positive(z) == false);
    REQUIRE(is_negative(z) == false);
    REQUIRE(is_nonnegative(z) == true);
    REQUIRE(is_nonpositive(z) == true);
    REQUIRE(is_finite(z) == true);
}

TEST_CASE("Integer(1): all standard predicates", "[2a][assumptions]") {
    // Reference: test_one
    auto z = integer(1);
    REQUIRE(is_positive(z) == true);
    REQUIRE(is_negative(z) == false);
    REQUIRE(is_zero(z) == false);
    REQUIRE(is_nonzero(z) == true);
    REQUIRE(is_nonnegative(z) == true);
    REQUIRE(is_nonpositive(z) == false);
    REQUIRE(is_integer(z) == true);
    REQUIRE(is_real(z) == true);
}

TEST_CASE("Integer(-1): all standard predicates", "[2a][assumptions]") {
    // Reference: test_negativeone
    auto z = integer(-1);
    REQUIRE(is_positive(z) == false);
    REQUIRE(is_negative(z) == true);
    REQUIRE(is_zero(z) == false);
    REQUIRE(is_nonpositive(z) == true);
    REQUIRE(is_nonnegative(z) == false);
}

TEST_CASE("Rational(1/2): positive non-integer real", "[2a][assumptions]") {
    auto r = rational(1, 2);
    REQUIRE(is_real(r) == true);
    REQUIRE(is_rational(r) == true);
    REQUIRE(is_integer(r) == false);
    REQUIRE(is_positive(r) == true);
    REQUIRE(is_nonzero(r) == true);
}

TEST_CASE("Rational(-3/4): negative", "[2a][assumptions]") {
    auto r = rational(-3, 4);
    REQUIRE(is_real(r) == true);
    REQUIRE(is_negative(r) == true);
    REQUIRE(is_positive(r) == false);
    REQUIRE(is_integer(r) == false);
}

TEST_CASE("Float(1.5): real, positive, but rational/integer unknown",
          "[2a][assumptions]") {
    // SymPy convention: Float assumptions for rational/integer return None.
    auto f = float_value(1.5);
    REQUIRE(is_real(f) == true);
    REQUIRE(is_positive(f) == true);
    REQUIRE_FALSE(is_rational(f).has_value());
    REQUIRE_FALSE(is_integer(f).has_value());
    REQUIRE(is_finite(f) == true);
}

// ----- NumberSymbol predicates ----------------------------------------------

TEST_CASE("S::Pi: positive real irrational finite", "[2a][assumptions]") {
    // Reference: test_pi
    auto pi = S::Pi();
    REQUIRE(is_real(pi) == true);
    REQUIRE(is_positive(pi) == true);
    REQUIRE(is_negative(pi) == false);
    REQUIRE(is_zero(pi) == false);
    REQUIRE(is_rational(pi) == false);
    REQUIRE(is_integer(pi) == false);
    REQUIRE(is_finite(pi) == true);
}

TEST_CASE("S::E: positive real irrational finite", "[2a][assumptions]") {
    // Reference: test_E
    auto e = S::E();
    REQUIRE(is_real(e) == true);
    REQUIRE(is_positive(e) == true);
    REQUIRE(is_rational(e) == false);
    REQUIRE(is_integer(e) == false);
}

TEST_CASE("S::EulerGamma and S::Catalan: positive real",
          "[2a][assumptions]") {
    REQUIRE(is_positive(S::EulerGamma()) == true);
    REQUIRE(is_real(S::EulerGamma()) == true);
    REQUIRE(is_positive(S::Catalan()) == true);
    REQUIRE(is_real(S::Catalan()) == true);
}

// ----- Add propagation -------------------------------------------------------

TEST_CASE("Add: sum of positives is positive", "[2b][assumptions][add]") {
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    auto y = symbol("y", AssumptionMask{}.set_positive(true));
    auto e = x + y;
    REQUIRE(is_positive(e) == true);
    REQUIRE(is_real(e) == true);
    REQUIRE(is_finite(e) == true);
    REQUIRE(is_nonzero(e) == true);
}

TEST_CASE("Add: sum of negatives is negative", "[2b][assumptions][add]") {
    auto x = symbol("x", AssumptionMask{}.set_negative(true));
    auto y = symbol("y", AssumptionMask{}.set_negative(true));
    auto e = x + y;
    REQUIRE(is_negative(e) == true);
    REQUIRE(is_real(e) == true);
}

TEST_CASE("Add: integer + integer = integer", "[2b][assumptions][add]") {
    auto m = symbol("m", AssumptionMask{}.set_integer(true));
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    auto e = m + n;
    REQUIRE(is_integer(e) == true);
    REQUIRE(is_rational(e) == true);
    REQUIRE(is_real(e) == true);
}

TEST_CASE("Add: positive + 1 is positive", "[2b][assumptions][add]") {
    // Reference: SymPy (Symbol('x', positive=True) + 1).is_positive == True
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    auto e = x + integer(1);
    REQUIRE(is_positive(e) == true);
}

TEST_CASE("Add: real + unknown stays Unknown", "[2b][assumptions][add]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    auto y = symbol("y");  // unknown
    auto e = x + y;
    REQUIRE_FALSE(is_real(e).has_value());
}

// ----- Mul propagation -------------------------------------------------------

TEST_CASE("Mul: product of positives is positive", "[2b][assumptions][mul]") {
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    auto y = symbol("y", AssumptionMask{}.set_positive(true));
    auto e = x * y;
    REQUIRE(is_positive(e) == true);
    REQUIRE(is_real(e) == true);
    REQUIRE(is_nonzero(e) == true);
}

TEST_CASE("Mul: positive * negative is negative", "[2b][assumptions][mul]") {
    auto p = symbol("p", AssumptionMask{}.set_positive(true));
    auto n = symbol("n", AssumptionMask{}.set_negative(true));
    auto e = p * n;
    REQUIRE(is_negative(e) == true);
}

TEST_CASE("Mul: negative * negative is positive", "[2b][assumptions][mul]") {
    auto a = symbol("a", AssumptionMask{}.set_negative(true));
    auto b = symbol("b", AssumptionMask{}.set_negative(true));
    auto e = a * b;
    REQUIRE(is_positive(e) == true);
}

TEST_CASE("Mul: any-zero factor zeroes the product", "[2b][assumptions][mul]") {
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    auto e = x * integer(0);
    REQUIRE(is_zero(e) == true);
}

TEST_CASE("Mul: product of integers is integer", "[2b][assumptions][mul]") {
    auto m = symbol("m", AssumptionMask{}.set_integer(true));
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    auto e = m * n;
    REQUIRE(is_integer(e) == true);
}

// ----- Pow propagation -------------------------------------------------------

TEST_CASE("Pow: positive base raised to anything is positive",
          "[2b][assumptions][pow]") {
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    auto y = symbol("y", AssumptionMask{}.set_real(true));
    auto e = pow(x, y);
    REQUIRE(is_positive(e) == true);
    REQUIRE(is_real(e) == true);
    REQUIRE(is_nonzero(e) == true);
}

TEST_CASE("Pow: integer base ^ nonneg integer exp is integer",
          "[2b][assumptions][pow]") {
    auto m = symbol("m", AssumptionMask{}.set_integer(true));
    auto n = symbol("n", AssumptionMask{}.set_integer(true).set_positive(true));
    auto e = pow(m, n);
    REQUIRE(is_integer(e) == true);
    REQUIRE(is_rational(e) == true);
}

TEST_CASE("Pow: rational base ^ integer exp is rational",
          "[2b][assumptions][pow]") {
    // "nonzero" is encoded as zero=false on the mask.
    auto q = symbol("q", AssumptionMask{}.set_rational(true).set_zero(false));
    // Note: positive isn't required — even negative integer exponent on a
    // nonzero rational is rational.
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    auto e = pow(q, n);
    REQUIRE(is_rational(e) == true);
}

// ----- refine() (Phase 2c) ---------------------------------------------------

TEST_CASE("refine: (x^a)^b -> x^(a*b) when x is positive",
          "[2c][refine]") {
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    auto a = symbol("a");
    auto b = symbol("b");
    auto e = pow(pow(x, a), b);
    auto r = refine(e);
    // Expected canonical: pow(x, a*b)
    REQUIRE(r->type_id() == TypeId::Pow);
    auto& p = static_cast<const Pow&>(*r);
    REQUIRE(p.base() == x);
    REQUIRE(p.exp() == mul(a, b));
}

TEST_CASE("refine: (x^2)^3 collapses to x^6 even without assumptions",
          "[2c][refine]") {
    // Both exponents are integers, so the rewrite is unconditionally valid.
    auto x = symbol("x");
    auto e = pow(pow(x, integer(2)), integer(3));
    auto r = refine(e);
    REQUIRE(r->type_id() == TypeId::Pow);
    auto& p = static_cast<const Pow&>(*r);
    REQUIRE(p.base() == x);
    REQUIRE(p.exp() == integer(6));
}

TEST_CASE("refine: (x^a)^b stays unevaluated without assumptions on x",
          "[2c][refine]") {
    auto x = symbol("x");
    auto a = symbol("a");
    auto b = symbol("b");
    auto e = pow(pow(x, a), b);
    auto r = refine(e);
    // Without x>0 and with non-integer exponents, refine cannot rewrite.
    REQUIRE(r == e);
}

TEST_CASE("refine: descends through Add and Mul", "[2c][refine]") {
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    auto y = symbol("y");
    auto a = symbol("a");
    auto b = symbol("b");

    // (x^a)^b + y -> x^(a*b) + y
    auto e = pow(pow(x, a), b) + y;
    auto r = refine(e);
    REQUIRE(r != e);
}

TEST_CASE("refine: leaves atomic and Number unchanged",
          "[2c][refine]") {
    REQUIRE(refine(symbol("x")) == symbol("x"));
    REQUIRE(refine(integer(7)) == integer(7));
    REQUIRE(refine(S::Pi()) == S::Pi());
}

TEST_CASE("refine: matches SymPy on assumption-gated rewrite",
          "[2c][refine][oracle]") {
    // SymPy: refine((x**a)**b, Q.positive(x)) == x**(a*b)
    // We've encoded x as a positive Symbol so the assumption is implicit.
    // Verify our rewritten output is structurally equal to SymPy's.
    auto& oracle = sympp::testing::Oracle::instance();
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    auto a = symbol("a");
    auto b = symbol("b");
    auto r = refine(pow(pow(x, a), b));
    REQUIRE(oracle.equivalent(r->str(), "x**(a*b)"));
}
