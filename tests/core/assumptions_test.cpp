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

#include <sympp/core/assumption_context.hpp>
#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number_symbol.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/refine.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/simplify/simplify.hpp>

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

TEST_CASE("AssumptionMask: even/odd closure (ASSUME-11)", "[2a][assumptions]") {
    // even ⇒ integer ⇒ rational ⇒ real; even excludes odd.
    auto e = close_assumptions(AssumptionMask{}.set_even(true));
    REQUIRE(e.get(AssumptionKey::Integer) == true);
    REQUIRE(e.get(AssumptionKey::Real) == true);
    REQUIRE(e.get(AssumptionKey::Odd) == false);
    // odd ⇒ integer, nonzero (zero=false); odd excludes even.
    auto o = close_assumptions(AssumptionMask{}.set_odd(true));
    REQUIRE(o.get(AssumptionKey::Integer) == true);
    REQUIRE(o.get(AssumptionKey::Nonzero) == true);
    REQUIRE(o.get(AssumptionKey::Even) == false);
    // zero ⇒ even (not odd).
    auto z = close_assumptions(AssumptionMask{}.set_zero(true));
    REQUIRE(z.get(AssumptionKey::Even) == true);
    REQUIRE(z.get(AssumptionKey::Odd) == false);
    // ¬integer ⇒ ¬even, ¬odd.
    auto ni = close_assumptions(AssumptionMask{}.set_integer(false));
    REQUIRE(ni.get(AssumptionKey::Even) == false);
    REQUIRE(ni.get(AssumptionKey::Odd) == false);
}

TEST_CASE("Symbol/Integer: even & odd predicates and consumers (ASSUME-11)",
          "[2b][assumptions][regression]") {
    auto e = symbol("e", AssumptionMask{}.set_even(true));
    auto o = symbol("o", AssumptionMask{}.set_odd(true));
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    REQUIRE(is_even(e) == true);
    REQUIRE(is_integer(e) == true);
    REQUIRE(is_odd(o) == true);
    REQUIRE(is_nonzero(o) == true);
    // A bare integer has unknown parity.
    REQUIRE(!is_even(n).has_value());
    // Integer literals answer parity directly.
    REQUIRE(is_even(integer(4)) == true);
    REQUIRE(is_odd(integer(7)) == true);
    // Declared parity drives the (-1)^k consumer (ASSUME-8); trig consumers are
    // covered in the trigonometric tests.
    REQUIRE(pow(S::NegativeOne(), e) == S::One());
    REQUIRE(pow(S::NegativeOne(), o) == S::NegativeOne());
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

TEST_CASE("Mul/Add/Pow: parity inference at the ask level (ASSUME-12)",
          "[2b][assumptions][regression]") {
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    auto m = symbol("m", AssumptionMask{}.set_integer(true));
    auto e = symbol("e", AssumptionMask{}.set_even(true));
    auto o = symbol("o", AssumptionMask{}.set_odd(true));
    // Product: even iff a factor is even; odd iff all factors odd.
    REQUIRE(is_even(integer(2) * n) == true);
    REQUIRE(is_even(o * e) == true);
    REQUIRE(is_odd(o * o) == true);          // o² (a Pow) — base parity preserved
    REQUIRE(is_even(o * o) == false);
    REQUIRE(!is_even(n * m).has_value());    // unknown parity factors
    // Sum: fold parities (XOR of odd-term count).
    REQUIRE(is_odd(integer(2) * n + integer(1)) == true);
    REQUIRE(is_even(o + o) == true);
    REQUIRE(is_odd(e + o) == true);
    // Pow: base^(positive int) preserves base parity.
    REQUIRE(is_even(pow(e, integer(3))) == true);
    REQUIRE(is_odd(pow(o, integer(5))) == true);
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

TEST_CASE("Pow: (-1)^k folds by structural parity (ASSUME-8)",
          "[2b][assumptions][pow][regression]") {
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    auto g = symbol("g");  // not known integer
    auto m1 = S::NegativeOne();
    // 2n even ⇒ 1; 2n+1, 4n+3 odd ⇒ -1; 2n+2 even ⇒ 1.
    REQUIRE(pow(m1, integer(2) * n) == S::One());
    REQUIRE(pow(m1, integer(2) * n + integer(1)) == m1);
    REQUIRE(pow(m1, integer(2) * n + integer(2)) == S::One());
    REQUIRE(pow(m1, integer(4) * n + integer(3)) == m1);
    // Unknown parity (bare integer) and non-integer coefficient stay.
    REQUIRE(pow(m1, n)->str() == "(-1)**n");
    REQUIRE(pow(m1, integer(2) * g)->str() == "(-1)**(2*g)");
}

TEST_CASE("Pow: integer base ^ nonneg integer exp is integer",
          "[2b][assumptions][pow]") {
    auto m = symbol("m", AssumptionMask{}.set_integer(true));
    auto n = symbol("n", AssumptionMask{}.set_integer(true).set_positive(true));
    auto e = pow(m, n);
    REQUIRE(is_integer(e) == true);
    REQUIRE(is_rational(e) == true);
}

TEST_CASE("Pow: real base ^ even integer exp is nonnegative (ASSUME-2)",
          "[2b][assumptions][pow][regression]") {
    auto r = symbol("r", AssumptionMask{}.set_real(true));
    // x² ≥ 0 for real x, but not provably > 0 (could be 0) and not nonpositive.
    REQUIRE(is_nonnegative(pow(r, integer(2))) == true);
    REQUIRE(!is_positive(pow(r, integer(2))).has_value());
    REQUIRE(is_nonnegative(pow(r, integer(4))) == true);
    // Cascades through Add: x² + 1 > 0 for real x.
    REQUIRE(is_positive(pow(r, integer(2)) + integer(1)) == true);
    // Real nonzero (zero=false) base ⇒ even power is strictly positive.
    auto rn = symbol("r", AssumptionMask{}.set_real(true).set_zero(false));
    REQUIRE(is_positive(pow(rn, integer(2))) == true);
    // Odd exponent: sign unknown. Generic (maybe complex) base: nothing inferred.
    REQUIRE(!is_nonnegative(pow(r, integer(3))).has_value());
    auto z = symbol("z");
    REQUIRE(!is_nonnegative(pow(z, integer(2))).has_value());
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

// ASSUMING-CONTEXT-1: a scoped `assuming(...)` context augments ask() globally —
// the fact derives related predicates, propagates through queries, unlocks
// refine, and retracts at scope exit. Mirrors SymPy's `assuming(Q...(x))`.
TEST_CASE("assuming: scoped relational facts drive ask and refine (ASSUMING-CONTEXT-1)",
          "[2c][refine][assuming][oracle]") {
    auto x = symbol("x");  // no built-in assumptions

    // Outside any scope the symbol is unconstrained.
    REQUIRE(!is_positive(x).has_value());
    REQUIRE(refine(abs(x)) == abs(x));
    REQUIRE(refine(sqrt(pow(x, integer(2)))) == sqrt(pow(x, integer(2))));

    {
        assuming pos(x, AssumptionMask{}.set_positive(true));
        // The fact and its consequences are visible.
        REQUIRE(is_positive(x) == true);
        REQUIRE(is_real(x) == true);      // positive ⇒ real
        REQUIRE(is_nonzero(x) == true);   // positive ⇒ nonzero
        // refine is unlocked: |x| → x, √(x²) → x.
        REQUIRE(refine(abs(x)) == x);
        REQUIRE(refine(sqrt(pow(x, integer(2)))) == x);
    }

    // Retracted at scope exit.
    REQUIRE(!is_positive(x).has_value());
    REQUIRE(refine(abs(x)) == abs(x));

    // A negative scope flips the abs sign (matches SymPy refine(Abs(x),Q.negative)).
    {
        assuming neg(x, AssumptionMask{}.set_negative(true));
        REQUIRE(refine(abs(x)) == mul(S::NegativeOne(), x));
        REQUIRE(is_real(x) == true);
    }

    // Scopes nest and apply to independent symbols.
    auto y = symbol("y");
    {
        assuming px(x, AssumptionMask{}.set_positive(true));
        assuming ny(y, AssumptionMask{}.set_negative(true));
        REQUIRE(refine(abs(x)) == x);
        REQUIRE(refine(abs(y)) == mul(S::NegativeOne(), y));
    }
    REQUIRE(!is_positive(x).has_value());
}

// REFINE-SQRT-SQUARE-1: √(g²) → |g| for a real/imaginary base (the principal
// root of a square is the absolute value, not g). A known sign then collapses it
// to ±g. SymPy leaves it as √(g²) without a real assumption, and does not
// generalize to other even roots like (x⁴)^(1/4).
TEST_CASE("refine: sqrt(x^2) -> |x| gated on reality (REFINE-SQRT-SQUARE-1)",
          "[2c][refine][assuming][oracle]") {
    auto& oracle = sympp::testing::Oracle::instance();
    auto x = symbol("x");
    Expr sx2 = sqrt(pow(x, integer(2)));

    // No reality fact → unchanged (matches SymPy).
    REQUIRE(refine(sx2) == sx2);

    {
        assuming re(x, AssumptionMask{}.set_real(true));
        REQUIRE(refine(sx2) == abs(x));  // → |x|
        // The general even root is not rewritten (matches SymPy).
        Expr q4 = pow(pow(x, integer(4)), rational(1, 4));
        REQUIRE(refine(q4) == q4);
    }
    {
        assuming neg(x, AssumptionMask{}.set_negative(true));
        REQUIRE(refine(sx2) == mul(S::NegativeOne(), x));  // → -x
    }
    {
        assuming pos(x, AssumptionMask{}.set_positive(true));
        REQUIRE(refine(sx2) == x);  // → x
    }
    // Oracle cross-check of the real-assumption form.
    {
        assuming re(x, AssumptionMask{}.set_real(true));
        REQUIRE(oracle.equivalent(refine(sx2)->str(), "Abs(x)"));
    }
}

// REFINE-ABS-MUL-1: |∏fᵢ| = ∏|fᵢ| — absolute value distributes over a product,
// so each factor's known sign applies independently. The distributed form is
// returned only when a factor actually collapses (|x·y| with no facts stays).
TEST_CASE("refine: abs distributes over a product (REFINE-ABS-MUL-1)",
          "[2c][refine][assuming][oracle]") {
    auto& oracle = sympp::testing::Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    Expr axy = abs(mul(x, y));

    // No sign facts → unchanged.
    REQUIRE(refine(axy) == axy);

    {  // Only x known positive → x·|y|.
        assuming px(x, AssumptionMask{}.set_positive(true));
        REQUIRE(oracle.equivalent(refine(axy)->str(), "x*Abs(y)"));
    }
    {  // Both positive → x·y.
        assuming px(x, AssumptionMask{}.set_positive(true));
        assuming py(y, AssumptionMask{}.set_positive(true));
        REQUIRE(refine(axy) == mul(x, y));
    }
    {  // x>0, y<0 → -x·y.
        assuming px(x, AssumptionMask{}.set_positive(true));
        assuming ny(y, AssumptionMask{}.set_negative(true));
        REQUIRE(refine(axy) == mul(S::NegativeOne(), mul(x, y)));
    }
}

// ASSUMING-PROPAGATE-1: a scoped fact reaches compound expressions — Mul and Add
// ask() query their children through the context-aware path, so x>0 ∧ y>0 makes
// x·y and x+y positive, etc. With no scope active behavior is unchanged.
TEST_CASE("assuming: facts propagate through products and sums (ASSUMING-PROPAGATE-1)",
          "[2c][assuming]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto prod = mul(x, y);
    auto sum = add(x, y);

    // Unconstrained outside any scope.
    REQUIRE(!is_positive(prod).has_value());
    REQUIRE(!is_real(sum).has_value());

    {
        assuming px(x, AssumptionMask{}.set_positive(true));
        assuming py(y, AssumptionMask{}.set_positive(true));
        REQUIRE(is_positive(prod) == true);    // (+)(+) > 0
        REQUIRE(is_real(prod) == true);        // positive ⇒ real
        REQUIRE(is_nonnegative(prod) == true);
        REQUIRE(is_positive(sum) == true);     // (+)+(+) > 0
    }
    {
        assuming nx(x, AssumptionMask{}.set_negative(true));
        assuming py(y, AssumptionMask{}.set_positive(true));
        REQUIRE(is_negative(prod) == true);    // (−)(+) < 0
    }
    {
        assuming rx(x, AssumptionMask{}.set_real(true));
        assuming ry(y, AssumptionMask{}.set_real(true));
        REQUIRE(is_real(prod) == true);
        REQUIRE(is_real(sum) == true);
    }

    // Powers propagate too (Pow::ask routed through the context-aware path).
    {
        assuming rx(x, AssumptionMask{}.set_real(true));
        REQUIRE(is_real(pow(x, integer(3))) == true);        // real^int is real
        REQUIRE(is_nonnegative(pow(x, integer(2))) == true);  // real² ≥ 0
    }
    {
        assuming px(x, AssumptionMask{}.set_positive(true));
        REQUIRE(is_positive(pow(x, integer(3))) == true);     // (+)^int > 0
        REQUIRE(is_positive(sqrt(x)) == true);                // √(+) > 0
    }

    // Retracted.
    REQUIRE(!is_positive(prod).has_value());
    REQUIRE(!is_real(sum).has_value());
    REQUIRE(!is_positive(pow(x, integer(3))).has_value());
}

// ASSUME-IMAG-1: Imaginary / Complex predicates with arithmetic propagation.
// imaginary ⇒ complex ∧ ¬real ∧ finite; I·(real≠0) is imaginary; I·I and
// imaginary² are real; sums of imaginaries are imaginary. Matches SymPy.
TEST_CASE("assumptions: imaginary / complex predicates (ASSUME-IMAG-1)",
          "[2][assumptions][imaginary][regression]") {
    using sympp::is_imaginary;
    using sympp::is_complex;
    using sympp::is_real;
    auto xr = symbol("xr", AssumptionMask{}.set_real(true).set_zero(false));
    auto yr = symbol("yr", AssumptionMask{}.set_real(true).set_zero(false));
    auto xi = symbol("xi", AssumptionMask{}.set_imaginary(true));
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};

    // I itself.
    REQUIRE(is_imaginary(S::I()) == T);
    REQUIRE(is_real(S::I()) == F);
    REQUIRE(is_complex(S::I()) == T);
    // A declared-imaginary symbol: imaginary ⇒ complex ∧ ¬real ∧ finite.
    REQUIRE(is_imaginary(xi) == T);
    REQUIRE(is_real(xi) == F);
    REQUIRE(sympp::is_finite(xi) == T);
    // Real numbers / symbols are complex but not imaginary; 0 is real.
    REQUIRE(is_imaginary(integer(2)) == F);
    REQUIRE(is_complex(integer(2)) == T);
    REQUIRE(is_imaginary(xr) == F);
    REQUIRE(is_imaginary(integer(0)) == F);
    REQUIRE(is_real(integer(0)) == T);
    // Mul propagation: I·(real ≠ 0) is imaginary and not real.
    REQUIRE(is_imaginary(S::I() * xr) == T);
    REQUIRE(is_real(S::I() * xr) == F);
    REQUIRE(is_imaginary(integer(2) * S::I()) == T);
    // Even number of imaginary factors → real: I·I, imaginary².
    REQUIRE(is_real(S::I() * S::I()) == T);          // folds to −1
    REQUIRE(is_imaginary(pow(xi, integer(2))) == F);
    REQUIRE(is_real(pow(xi, integer(2))) == T);
    REQUIRE(is_imaginary(pow(xi, integer(3))) == T);  // odd power stays imaginary
    // real · real is not imaginary.
    REQUIRE(is_imaginary(xr * yr) == F);
    // Add: imaginary + imaginary is imaginary; real + imaginary is not.
    REQUIRE(is_imaginary(S::I() * xr + S::I() * yr) == T);
    REQUIRE(is_imaginary(xr + S::I() * yr) == F);
    REQUIRE(is_complex(xr + S::I() * yr) == T);
}

// ASSUME-NONNEG-1: nonnegative (x ≥ 0) and nonpositive (x ≤ 0) are declarable
// primary facts, not just positive∨zero. Previously set(Nonnegative,true) only
// recorded negative=false, losing the fact (is_nonnegative → None). Now they close
// to real/finite, exclude the opposite strict sign, and refine to the strict sign
// once the zero case is ruled in or out. Matches SymPy's nonnegative/nonpositive.
TEST_CASE("AssumptionMask: declarable nonnegative / nonpositive (ASSUME-NONNEG-1)",
          "[2a][assumptions]") {
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};

    auto x = symbol("x", AssumptionMask{}.set_nonnegative(true));
    REQUIRE(is_nonnegative(x) == T);
    REQUIRE(is_real(x) == T);
    REQUIRE(is_negative(x) == F);
    REQUIRE(is_positive(x) == U);   // could be 0
    REQUIRE(is_zero(x) == U);

    auto w = symbol("w", AssumptionMask{}.set_nonpositive(true));
    REQUIRE(is_nonpositive(w) == T);
    REQUIRE(is_real(w) == T);
    REQUIRE(is_positive(w) == F);

    // Refinement: nonnegative ∧ nonzero ⇒ positive; nonnegative ∧ ¬positive ⇒ zero.
    auto y = symbol("y", AssumptionMask{}.set_nonnegative(true).set_nonzero(true));
    REQUIRE(is_positive(y) == T);

    // Consistency with the primaries: positive ⇒ nonnegative, negative ⇒ nonpositive.
    auto p = symbol("p", AssumptionMask{}.set_positive(true));
    REQUIRE(is_nonnegative(p) == T);
    REQUIRE(is_nonpositive(p) == F);

    // Order predicates are defined on ℝ: imaginary ⇒ ¬nonnegative, ¬nonpositive.
    auto i = symbol("i", AssumptionMask{}.set_imaginary(true));
    REQUIRE(is_nonnegative(i) == F);
    REQUIRE(is_nonpositive(i) == F);
}

// ASSUME-NONNEG-2: the declared nonnegativity flows into simplification —
// √(x²) → x and |x| → x for x ≥ 0 (and the nonpositive analogues).
TEST_CASE("simplify: uses declared nonnegative/nonpositive (ASSUME-NONNEG-2)",
          "[2a][assumptions][simplify]") {
    auto x = symbol("x", AssumptionMask{}.set_nonnegative(true));
    auto w = symbol("w", AssumptionMask{}.set_nonpositive(true));
    REQUIRE(simplify(pow(pow(x, integer(2)), rational(1, 2))) == x);
    REQUIRE(simplify(abs(x)) == x);
    REQUIRE(simplify(abs(w)) == mul(S::NegativeOne(), w));
}

// ASSUME-MULSIGN-1: sign direction of a product propagates through nonnegative /
// nonpositive factors (not just strict signs). positive·nonnegative → nonnegative,
// nonpositive·nonpositive → nonnegative, positive·nonpositive → nonpositive.
// Conservative: the unprovable opposite (which would need ruling out a zero factor)
// stays unknown. Matches SymPy.
TEST_CASE("Mul: nonnegative/nonpositive sign propagation (ASSUME-MULSIGN-1)",
          "[2a][assumptions]") {
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};
    auto a = symbol("a", AssumptionMask{}.set_positive(true));
    auto c = symbol("c", AssumptionMask{}.set_negative(true));
    auto nn = symbol("nn", AssumptionMask{}.set_nonnegative(true));
    auto np = symbol("np", AssumptionMask{}.set_nonpositive(true));
    REQUIRE(is_nonnegative(a * nn) == T);          // >0 · ≥0
    REQUIRE(is_nonnegative(c * c) == T);           // <0 · <0
    REQUIRE(is_nonnegative(np * np) == T);         // ≤0 · ≤0
    REQUIRE(is_nonpositive(a * np) == T);          // >0 · ≤0
    REQUIRE(is_nonpositive(a * c) == T);           // >0 · <0
    REQUIRE(is_nonnegative(a * c) == F);           // strictly < 0
    REQUIRE(is_nonnegative(a * np) == U);          // could be 0 or negative
}

// ASSUME-REALFINITE-1: a symbol declared real denotes a finite real number, so
// real ⇒ finite (the unbounded ±∞ are the separate Infinity atoms). Matches SymPy
// and is consistent with positive/negative/zero ⇒ finite. Closes exp(real) → finite.
TEST_CASE("AssumptionMask: real => finite (ASSUME-REALFINITE-1)",
          "[2a][assumptions]") {
    const auto T = std::optional<bool>{true};
    auto r = symbol("r", AssumptionMask{}.set_real(true));
    REQUIRE(is_finite(r) == T);
    REQUIRE(is_finite(exp(r)) == T);
    REQUIRE(is_nonzero(exp(r)) == T);
}

// PRIME-1: the `prime` predicate. Concrete integers answer via GMP primality
// (n < 2 → False); a declared-prime symbol closes to integer ∧ positive ∧
// nonzero ∧ real ∧ rational, but NOT odd (2 is prime and even). A non-integer
// is never prime. Matches SymPy's Q.prime / Symbol(..., prime=True).
TEST_CASE("AssumptionMask: prime predicate (PRIME-1)", "[2a][assumptions]") {
    using sympp::is_prime;
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};

    // Concrete integers.
    REQUIRE(is_prime(integer(2)) == T);
    REQUIRE(is_prime(integer(3)) == T);
    REQUIRE(is_prime(integer(11)) == T);
    REQUIRE(is_prime(integer(1)) == F);
    REQUIRE(is_prime(integer(4)) == F);
    REQUIRE(is_prime(integer(8)) == F);
    REQUIRE(is_prime(integer(0)) == F);
    REQUIRE(is_prime(integer(-3)) == F);

    // Declared-prime symbol: prime ⇒ integer, positive, nonzero, real, rational.
    auto p = symbol("p", AssumptionMask{}.set_prime(true));
    REQUIRE(is_prime(p) == T);
    REQUIRE(is_integer(p) == T);
    REQUIRE(is_positive(p) == T);
    REQUIRE(is_nonzero(p) == T);
    REQUIRE(is_real(p) == T);
    REQUIRE(is_rational(p) == T);
    REQUIRE(is_even(p) == U);  // 2 is prime and even → parity is not pinned
    REQUIRE(is_odd(p) == U);

    // Non-integers are never prime; an integer-unknown symbol stays Unknown.
    auto q = symbol("q", AssumptionMask{}.set_rational(true).set_integer(false));
    REQUIRE(is_prime(q) == F);
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(is_prime(x) == U);
    REQUIRE(is_prime(S::I()) == F);
}

// COMPOSITE-1: the `composite` predicate (integer > 1 that is not prime, so
// ≥ 4). Concrete integers answer via GMP primality (n < 4 → False); a
// declared-composite symbol closes to integer ∧ positive ∧ nonzero ∧ real ∧
// ¬prime, but NOT a fixed parity (4 even, 9 odd). Prime and composite are
// mutually exclusive; 1 is neither. Matches SymPy's Q.composite.
TEST_CASE("AssumptionMask: composite predicate (COMPOSITE-1)", "[2a][assumptions]") {
    using sympp::is_composite;
    using sympp::is_prime;
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};

    // Concrete integers: composite iff ≥ 4 and not prime.
    REQUIRE(is_composite(integer(4)) == T);
    REQUIRE(is_composite(integer(6)) == T);
    REQUIRE(is_composite(integer(9)) == T);
    REQUIRE(is_composite(integer(2)) == F);   // prime
    REQUIRE(is_composite(integer(3)) == F);   // prime
    REQUIRE(is_composite(integer(5)) == F);   // prime
    REQUIRE(is_composite(integer(1)) == F);   // neither prime nor composite
    REQUIRE(is_composite(integer(0)) == F);
    REQUIRE(is_composite(integer(-4)) == F);  // negatives are not composite

    // Declared-composite symbol: composite ⇒ integer, positive, nonzero, real,
    // ¬prime; parity stays Unknown.
    auto c = symbol("c", AssumptionMask{}.set_composite(true));
    REQUIRE(is_composite(c) == T);
    REQUIRE(is_integer(c) == T);
    REQUIRE(is_positive(c) == T);
    REQUIRE(is_nonzero(c) == T);
    REQUIRE(is_real(c) == T);
    REQUIRE(is_prime(c) == F);   // mutually exclusive
    REQUIRE(is_even(c) == U);
    REQUIRE(is_odd(c) == U);

    // A declared-prime symbol is not composite; non-integers never composite.
    auto p = symbol("p", AssumptionMask{}.set_prime(true));
    REQUIRE(is_composite(p) == F);
    auto q = symbol("q", AssumptionMask{}.set_rational(true).set_integer(false));
    REQUIRE(is_composite(q) == F);
    auto i = symbol("i", AssumptionMask{}.set_integer(true).set_positive(true));
    REQUIRE(is_composite(i) == U);  // a generic positive integer is unknown
    REQUIRE(is_composite(S::I()) == F);
}

// IRRATIONAL-1: the `irrational` predicate, defined as real ∧ ¬rational. The
// proven-irrational constants π and e answer True; integers/rationals answer
// False; I and the infinities are not real so answer False. EulerGamma and
// Catalan have unknown rationality (open problems) so every rationality-
// dependent predicate is Unknown for them — this also corrected their previous
// (wrong) rational=False answer. A declared-irrational symbol closes to real ∧
// finite ∧ nonzero ∧ ¬rational ∧ ¬integer ∧ ¬parity ∧ complex ∧ ¬imaginary.
// Matches SymPy's Q.irrational / Symbol(..., irrational=True).
TEST_CASE("AssumptionMask: irrational predicate (IRRATIONAL-1)", "[2a][assumptions]") {
    using sympp::is_irrational;
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};

    // Concrete rationals are not irrational; π and e are.
    REQUIRE(is_irrational(integer(3)) == F);
    REQUIRE(is_irrational(integer(0)) == F);
    REQUIRE(is_irrational(rational(3, 2)) == F);
    REQUIRE(is_irrational(S::Pi()) == T);
    REQUIRE(is_irrational(S::E()) == T);
    REQUIRE(is_irrational(S::I()) == F);          // not real
    REQUIRE(is_irrational(S::Infinity()) == F);   // not a real number

    // EulerGamma / Catalan: rationality is an open problem ⇒ Unknown (and the
    // dependent predicates are Unknown too, not False).
    REQUIRE(is_irrational(S::EulerGamma()) == U);
    REQUIRE(is_rational(S::EulerGamma()) == U);
    REQUIRE(is_integer(S::EulerGamma()) == U);
    REQUIRE(is_irrational(S::Catalan()) == U);
    REQUIRE(is_real(S::EulerGamma()) == T);       // still a positive real
    REQUIRE(is_positive(S::EulerGamma()) == T);

    // Declared-irrational symbol: irrational ⇒ real ∧ ¬rational ∧ finite ∧
    // nonzero ∧ ¬integer ∧ ¬parity ∧ complex ∧ ¬imaginary.
    auto q = symbol("q", AssumptionMask{}.set_irrational(true));
    REQUIRE(is_irrational(q) == T);
    REQUIRE(is_real(q) == T);
    REQUIRE(is_rational(q) == F);
    REQUIRE(is_integer(q) == F);
    REQUIRE(is_even(q) == F);
    REQUIRE(is_odd(q) == F);
    REQUIRE(is_nonzero(q) == T);
    REQUIRE(is_finite(q) == T);
    REQUIRE(is_complex(q) == T);
    REQUIRE(is_imaginary(q) == F);

    // Biconditional: real ∧ ¬rational ⇒ irrational; a plain real symbol is
    // Unknown; rational / integer symbols are not irrational.
    auto x = symbol("x", AssumptionMask{}.set_real(true).set_rational(false));
    REQUIRE(is_irrational(x) == T);
    auto r = symbol("r", AssumptionMask{}.set_real(true));
    REQUIRE(is_irrational(r) == U);
    auto ra = symbol("ra", AssumptionMask{}.set_rational(true));
    REQUIRE(is_irrational(ra) == F);
    auto in = symbol("in", AssumptionMask{}.set_integer(true));
    REQUIRE(is_irrational(in) == F);
}

// IRRATIONAL-ROOT-1: a surviving radical of a positive rational base — √2,
// 2^(1/3), 6^(1/2) — is a real algebraic irrational. The Pow factory reduces
// every perfect rational root (√4 → 2, √(9/4) → 3/2), so any Pow that persists
// with a positive rational base ≠ 1 and a non-integer rational exponent is not
// rational. Previously these answered Unknown for irrational / rational.
// Matches SymPy (sqrt(2).is_irrational == True).
TEST_CASE("AssumptionMask: irrational algebraic roots (IRRATIONAL-ROOT-1)",
          "[2a][assumptions]") {
    using sympp::is_algebraic;
    using sympp::is_irrational;
    using sympp::is_rational;
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};

    auto root = [](long b, long p, long q) {
        return pow(integer(b), rational(p, q));
    };
    // √2, 2^(1/3), 6^(1/2): irrational, not rational, algebraic, real. (A
    // reciprocal/scaled root like 1/√2 = √2/2 is a Mul and needs Mul-level
    // irrationality propagation, which is left for a follow-up.)
    for (const Expr& e : {root(2, 1, 2), root(3, 1, 3), root(6, 1, 2)}) {
        REQUIRE(is_irrational(e) == T);
        REQUIRE(is_rational(e) == F);
        REQUIRE(is_integer(e) != T);  // not provably an integer
        REQUIRE(is_algebraic(e) == T);
        REQUIRE(is_real(e) == T);
    }
    // A perfect root reduces to an integer/rational and stays rational.
    REQUIRE(is_rational(root(4, 1, 2)) == T);   // √4 → 2
    REQUIRE(is_irrational(root(4, 1, 2)) == F);
    REQUIRE(is_rational(pow(rational(9, 4), rational(1, 2))) == T);  // → 3/2
    // A negative base gives an imaginary value, not an irrational real.
    REQUIRE(is_irrational(root(-2, 1, 2)) != T);
    REQUIRE(is_real(root(-2, 1, 2)) == F);
}

// IRRATIONAL-PROP-1: irrationality propagates through a Mul/Add against rational
// terms — a nonzero rational times exactly one irrational factor is irrational
// (2π, √2/2, 3√2), and a rational sum plus exactly one irrational term is
// irrational (π + 1, √2 − 3). A real transcendental (π, e) counts as the
// irrational factor. Two irrational operands leave it Unknown (√2 + √3, since the
// sum could in principle be rational). Matches SymPy.
TEST_CASE("AssumptionMask: Mul/Add irrationality propagation (IRRATIONAL-PROP-1)",
          "[2a][assumptions]") {
    using sympp::is_irrational;
    using sympp::is_rational;
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};
    auto s2 = pow(integer(2), rational(1, 2));
    auto s3 = pow(integer(3), rational(1, 2));

    // nonzero rational × one irrational → irrational.
    for (const Expr& e : {mul(integer(2), S::Pi()),
                          mul(rational(1, 2), S::Pi()),
                          mul(rational(2, 3), S::Pi()), mul(integer(3), s2),
                          pow(integer(2), rational(-1, 2))}) {  // 1/√2 = √2/2
        REQUIRE(is_irrational(e) == T);
        REQUIRE(is_rational(e) == F);
    }
    // rational + one irrational → irrational.
    for (const Expr& e : {add(S::Pi(), integer(1)), add(S::Pi(), integer(-3)),
                          add(s2, integer(1)), add(mul(integer(2), S::E()),
                                                   integer(5))}) {
        REQUIRE(is_irrational(e) == T);
        REQUIRE(is_rational(e) == F);
    }
    // Two irrational operands: Unknown (the combination may be rational).
    REQUIRE(is_irrational(add(s2, s3)) == U);
    REQUIRE(is_irrational(mul(S::Pi(), S::E())) == U);
    // EulerGamma's rationality is open, so a sum/product with it stays Unknown.
    REQUIRE(is_irrational(add(S::EulerGamma(), integer(1))) == U);
}

// ALGTRANS-1: the `algebraic` / `transcendental` pair. algebraic = root of a
// rational polynomial (implies complex, finite; rational ⇒ algebraic; i is
// algebraic but not real); transcendental ⟺ complex ∧ ¬algebraic (implies
// complex, finite, ¬algebraic, ¬rational; a real transcendental is irrational).
// π and e are transcendental; EulerGamma/Catalan are unknown. Matches SymPy.
TEST_CASE("AssumptionMask: algebraic / transcendental (ALGTRANS-1)", "[2a][assumptions]") {
    using sympp::is_algebraic;
    using sympp::is_transcendental;
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};

    // Concrete numbers.
    REQUIRE(is_algebraic(integer(3)) == T);
    REQUIRE(is_transcendental(integer(3)) == F);
    REQUIRE(is_algebraic(rational(3, 2)) == T);
    REQUIRE(is_algebraic(S::Pi()) == F);
    REQUIRE(is_transcendental(S::Pi()) == T);
    REQUIRE(is_transcendental(S::E()) == T);
    REQUIRE(is_algebraic(S::I()) == T);            // i is a root of x² + 1
    REQUIRE(is_transcendental(S::I()) == F);
    REQUIRE(is_algebraic(S::Infinity()) == F);     // not a finite complex number
    REQUIRE(is_transcendental(S::Infinity()) == F);

    // EulerGamma / Catalan: neither status is known.
    REQUIRE(is_algebraic(S::EulerGamma()) == U);
    REQUIRE(is_transcendental(S::EulerGamma()) == U);
    REQUIRE(is_algebraic(S::Catalan()) == U);

    // Declared-algebraic symbol: complex ∧ finite ∧ ¬transcendental, but real /
    // rational / integer / imaginary all stay Unknown (i is algebraic).
    auto a = symbol("a", AssumptionMask{}.set_algebraic(true));
    REQUIRE(is_algebraic(a) == T);
    REQUIRE(is_transcendental(a) == F);
    REQUIRE(is_complex(a) == T);
    REQUIRE(is_finite(a) == T);
    REQUIRE(is_real(a) == U);
    REQUIRE(is_rational(a) == U);
    REQUIRE(is_imaginary(a) == U);

    // Declared-transcendental symbol: complex ∧ finite ∧ ¬algebraic ∧ ¬rational
    // ∧ ¬integer; real stays Unknown, so irrational is Unknown.
    auto tr = symbol("tr", AssumptionMask{}.set_transcendental(true));
    REQUIRE(is_transcendental(tr) == T);
    REQUIRE(is_algebraic(tr) == F);
    REQUIRE(is_complex(tr) == T);
    REQUIRE(is_finite(tr) == T);
    REQUIRE(is_rational(tr) == F);
    REQUIRE(is_integer(tr) == F);
    REQUIRE(is_real(tr) == U);
    REQUIRE(is_irrational(tr) == U);

    // A real transcendental is irrational; rational/integer ⇒ algebraic; an
    // imaginary symbol's algebraicity is unknown.
    auto rt = symbol("rt", AssumptionMask{}.set_transcendental(true).set_real(true));
    REQUIRE(is_irrational(rt) == T);
    REQUIRE(is_algebraic(symbol("r", AssumptionMask{}.set_rational(true))) == T);
    REQUIRE(is_algebraic(symbol("i", AssumptionMask{}.set_integer(true))) == T);
    REQUIRE(is_algebraic(symbol("im", AssumptionMask{}.set_imaginary(true))) == U);
}

// EXTREAL-INF-1: the `extended_real` / `infinite` predicates. extended_real is a
// point of ℝ ∪ {±∞} (real ⇒ extended_real, but ±∞ are extended-real yet not
// real); it excludes imaginary and does not imply finite/complex. infinite ⟺
// ¬finite, and implies ¬real/¬complex/¬zero. The ∞ atoms answer directly (oo/-oo
// extended-real, zoo not); a declared-complex symbol is now finite (so not
// infinite), matching the enum's "finite complex number" contract. Matches SymPy.
TEST_CASE("AssumptionMask: extended_real / infinite (EXTREAL-INF-1)", "[2a][assumptions]") {
    using sympp::is_extended_real;
    using sympp::is_infinite;
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};

    // Finite reals: extended-real, not infinite.
    REQUIRE(is_extended_real(integer(3)) == T);
    REQUIRE(is_infinite(integer(3)) == F);
    REQUIRE(is_extended_real(S::Pi()) == T);
    // I is off the real line and finite.
    REQUIRE(is_extended_real(S::I()) == F);
    REQUIRE(is_infinite(S::I()) == F);
    // The infinity atoms.
    REQUIRE(is_extended_real(S::Infinity()) == T);
    REQUIRE(is_infinite(S::Infinity()) == T);
    REQUIRE(is_finite(S::Infinity()) == F);
    REQUIRE(is_extended_real(S::NegativeInfinity()) == T);
    REQUIRE(is_infinite(S::NegativeInfinity()) == T);
    REQUIRE(is_extended_real(S::ComplexInfinity()) == F);  // zoo is off the line
    REQUIRE(is_infinite(S::ComplexInfinity()) == T);

    // real ⇒ extended_real ∧ finite ∧ ¬infinite.
    auto re = symbol("re", AssumptionMask{}.set_real(true));
    REQUIRE(is_extended_real(re) == T);
    REQUIRE(is_infinite(re) == F);

    // Declared-infinite symbol: ¬finite ∧ ¬real ∧ ¬complex ∧ ¬zero; extended_real
    // stays Unknown (could be ±∞ or zoo).
    auto inf = symbol("inf", AssumptionMask{}.set_infinite(true));
    REQUIRE(is_infinite(inf) == T);
    REQUIRE(is_finite(inf) == F);
    REQUIRE(is_real(inf) == F);
    REQUIRE(is_complex(inf) == F);
    REQUIRE(is_extended_real(inf) == U);

    // Declared-extended_real symbol: only ¬imaginary follows; real/finite/complex
    // stay Unknown (it might be ±∞).
    auto er = symbol("er", AssumptionMask{}.set_extended_real(true));
    REQUIRE(is_extended_real(er) == T);
    REQUIRE(is_imaginary(er) == F);
    REQUIRE(is_real(er) == U);
    REQUIRE(is_finite(er) == U);
    REQUIRE(is_infinite(er) == U);

    // imaginary ⇒ ¬extended_real; a complex symbol is finite ⇒ not infinite.
    REQUIRE(is_extended_real(symbol("im", AssumptionMask{}.set_imaginary(true))) == F);
    REQUIRE(is_infinite(symbol("cx", AssumptionMask{}.set_complex(true))) == F);
    REQUIRE(is_finite(symbol("cx2", AssumptionMask{}.set_complex(true))) == T);
}

// ALGCLOSURE-1: the algebraic numbers form a field, so Add and Mul propagate
// algebraicity. A sum/product of algebraic args is algebraic; a sum of
// algebraics plus exactly one transcendental term is transcendental (and not
// algebraic); a product of nonzero algebraics with exactly one transcendental
// factor is transcendental. Previously Add/Mul returned Unknown for these.
// Matches SymPy.
TEST_CASE("Add/Mul: algebraic field closure (ALGCLOSURE-1)",
          "[2b][assumptions][algebraic]") {
    using sympp::is_algebraic;
    using sympp::is_transcendental;
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};

    auto a = symbol("a", AssumptionMask{}.set_algebraic(true));
    auto b = symbol("b", AssumptionMask{}.set_algebraic(true));
    auto tr = symbol("tr", AssumptionMask{}.set_transcendental(true));
    auto tu = symbol("tu", AssumptionMask{}.set_transcendental(true));
    auto anz = symbol("anz", AssumptionMask{}.set_algebraic(true).set_nonzero(true));

    // Closure under + and ×, including the algebraic atom I and integers.
    REQUIRE(is_algebraic(a + b) == T);
    REQUIRE(is_transcendental(a + b) == F);
    REQUIRE(is_algebraic(a * b) == T);
    REQUIRE(is_algebraic(a * S::I()) == T);          // i is algebraic
    REQUIRE(is_algebraic(a + b + S::I()) == T);
    REQUIRE(is_algebraic(integer(2) + integer(3) * a) == T);

    // Algebraics + exactly one transcendental ⇒ transcendental.
    REQUIRE(is_algebraic(a + tr) == F);
    REQUIRE(is_transcendental(a + tr) == T);

    // Product: needs the algebraic factor nonzero (0·t = 0 is algebraic).
    REQUIRE(is_algebraic(a * tr) == U);              // a might be zero
    REQUIRE(is_transcendental(a * tr) == U);
    REQUIRE(is_algebraic(anz * tr) == F);            // nonzero ⇒ transcendental
    REQUIRE(is_transcendental(anz * tr) == T);

    // Two transcendentals can cancel — Unknown.
    REQUIRE(is_algebraic(tr + tu) == U);

    // A plain symbolic sum stays Unknown.
    REQUIRE(is_algebraic(symbol("x") + symbol("y")) == U);
}

// ALGCLOSURE-POW-1: a rational power of an algebraic number is algebraic
// (radicals: sqrt(2) = 2^(1/2)). The 0^negative case is guarded by requiring a
// nonnegative exponent or a nonzero base. A transcendental base to a nonzero
// rational power is transcendental. Previously Pow returned Unknown for these.
// Matches SymPy.
TEST_CASE("Pow: algebraic of rational power (ALGCLOSURE-POW-1)",
          "[2b][assumptions][algebraic]") {
    using sympp::is_algebraic;
    using sympp::is_transcendental;
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};

    auto a = symbol("a", AssumptionMask{}.set_algebraic(true));
    auto anz = symbol("anz", AssumptionMask{}.set_algebraic(true).set_nonzero(true));
    auto tr = symbol("tr", AssumptionMask{}.set_transcendental(true));
    auto n = symbol("n", AssumptionMask{}.set_integer(true));

    // Radicals and rational powers of algebraic numbers.
    REQUIRE(is_algebraic(sqrt(integer(2))) == T);          // 2^(1/2)
    REQUIRE(is_algebraic(pow(integer(2), rational(1, 3))) == T);
    REQUIRE(is_algebraic(pow(a, rational(3, 2))) == T);
    REQUIRE(is_algebraic(pow(a, integer(2))) == T);

    // 0^negative guard: a might be zero ⇒ Unknown; a nonzero base ⇒ algebraic.
    REQUIRE(is_algebraic(pow(a, rational(-1, 2))) == U);
    REQUIRE(is_algebraic(pow(anz, rational(-1, 2))) == T);

    // A symbolic integer exponent could be negative with a possibly-zero base.
    REQUIRE(is_algebraic(pow(a, n)) == U);

    // Transcendental base, nonzero rational exponent ⇒ transcendental.
    REQUIRE(is_algebraic(pow(tr, integer(2))) == F);
    REQUIRE(is_transcendental(pow(tr, integer(2))) == T);
    REQUIRE(is_transcendental(pow(tr, rational(1, 2))) == T);

    // Algebraic base, transcendental exponent ⇒ Unknown (Gelfond out of scope).
    REQUIRE(is_algebraic(pow(integer(2), tr)) == U);

    // Composition: (1 + sqrt(2))^3 is algebraic via Add-closure then Pow.
    REQUIRE(is_algebraic(pow(integer(1) + sqrt(integer(2)), integer(3))) == T);
}

// FUNC-SIGN-1: exp is never zero, and the odd strictly-increasing hyperbolics
// sinh/tanh preserve the argument's sign. Previously these were Unknown. Matches
// SymPy (exp(x).is_zero == False, sinh(p>0).is_positive == True, …).
TEST_CASE("AssumptionMask: exp-nonzero and hyperbolic sign propagation (FUNC-SIGN-1)",
          "[2a][assumptions]") {
    using sympp::is_negative;
    using sympp::is_nonzero;
    using sympp::is_positive;
    using sympp::is_real;
    using sympp::is_zero;
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};
    auto p = symbol("p", AssumptionMask{}.set_positive(true));
    auto n = symbol("n", AssumptionMask{}.set_negative(true));
    auto nz = symbol("nz", AssumptionMask{}.set_real(true).set_nonzero(true));
    auto im = symbol("im", AssumptionMask{}.set_imaginary(true));
    auto x = symbol("x");

    // exp is never zero, for any argument.
    REQUIRE(is_zero(exp(x)) == F);
    REQUIRE(is_zero(exp(im)) == F);
    REQUIRE(is_nonzero(exp(p)) == T);  // finite real ⇒ nonzero

    // sinh / tanh preserve sign.
    REQUIRE(is_positive(sinh(p)) == T);
    REQUIRE(is_negative(sinh(n)) == T);
    REQUIRE(is_positive(tanh(p)) == T);
    REQUIRE(is_negative(tanh(n)) == T);
    REQUIRE(is_nonzero(sinh(nz)) == T);
    REQUIRE(is_nonzero(tanh(nz)) == T);
    // Derived nonnegative/nonpositive follow.
    REQUIRE(is_nonnegative(tanh(p)) == T);
    REQUIRE(is_nonpositive(sinh(n)) == T);

    // No over-reach: a real argument of unknown sign, or an imaginary argument,
    // leaves the sign / zero Unknown; reality of the result still holds for real.
    REQUIRE(is_positive(sinh(symbol("r", AssumptionMask{}.set_real(true)))) == U);
    REQUIRE(is_zero(sinh(symbol("r2", AssumptionMask{}.set_real(true)))) == U);
    REQUIRE(is_positive(sinh(im)) == U);
    REQUIRE(is_real(sinh(p)) == T);
    REQUIRE(is_positive(cosh(symbol("rc", AssumptionMask{}.set_real(true))))
            == T);  // cosh unchanged
}

// FUNC-SIGN-2: atan is odd, strictly increasing on ℝ, and always real and finite
// (bounded in (−π/2, π/2)), so it preserves the argument's sign. Previously the
// sign was Unknown. Matches SymPy (atan(p>0).is_positive == True).
TEST_CASE("AssumptionMask: atan sign propagation (FUNC-SIGN-2)",
          "[2a][assumptions]") {
    using sympp::is_finite;
    using sympp::is_negative;
    using sympp::is_nonzero;
    using sympp::is_positive;
    using sympp::is_real;
    using sympp::is_zero;
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};
    auto p = symbol("p", AssumptionMask{}.set_positive(true));
    auto n = symbol("n", AssumptionMask{}.set_negative(true));
    auto nz = symbol("nz", AssumptionMask{}.set_real(true).set_nonzero(true));
    auto r = symbol("r", AssumptionMask{}.set_real(true));
    auto im = symbol("im", AssumptionMask{}.set_imaginary(true));

    REQUIRE(is_positive(atan(p)) == T);
    REQUIRE(is_negative(atan(n)) == T);
    REQUIRE(is_zero(atan(p)) == F);
    REQUIRE(is_nonzero(atan(nz)) == T);
    REQUIRE(is_nonnegative(atan(p)) == T);  // derived
    // atan(real) is always real and finite (the function is bounded).
    REQUIRE(is_real(atan(r)) == T);
    REQUIRE(is_finite(atan(r)) == T);
    // No over-reach: imaginary or unknown-sign-real arguments stay Unknown.
    REQUIRE(is_positive(atan(im)) == U);
    REQUIRE(is_positive(atan(r)) == U);
}

// MUL-CONSEC-EVEN-1: a product of integers containing two consecutive factors —
// n·(n−1), n·(n+1)·(n+2), 2n·(2n+1) — is even, since one of two consecutive
// integers is even. Previously Unknown. Matches SymPy ((n*(n-1)).is_even == True),
// which likewise infers this only for the factored product (the expanded n²−n
// stays Unknown in both systems).
TEST_CASE("AssumptionMask: product of consecutive integers is even (MUL-CONSEC-EVEN-1)",
          "[2a][assumptions]") {
    using sympp::is_even;
    using sympp::is_odd;
    const auto T = std::optional<bool>{true};
    const auto F = std::optional<bool>{false};
    const auto U = std::optional<bool>{};
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    auto m = symbol("m", AssumptionMask{}.set_integer(true));
    auto m1 = [&](const Expr& a, long c) { return add(a, integer(c)); };

    REQUIRE(is_even(mul(n, m1(n, -1))) == T);          // n·(n−1)
    REQUIRE(is_even(mul(n, m1(n, 1))) == T);           // n·(n+1)
    REQUIRE(is_odd(mul(n, m1(n, -1))) == F);
    REQUIRE(is_even(mul(mul(n, m1(n, 1)), m1(n, 2))) == T);  // n(n+1)(n+2)
    // Scaled non-constant part: 2n and 2n+1 are consecutive.
    auto twon = mul(integer(2), n);
    REQUIRE(is_even(mul(twon, m1(twon, 1))) == T);
    // A spectator integer factor does not break the inference.
    REQUIRE(is_even(mul(mul(n, m1(n, 1)), integer(7))) == T);

    // Not consecutive (offset 2) or different non-constant parts ⇒ Unknown.
    REQUIRE(is_even(mul(n, m1(n, -2))) == U);
    REQUIRE(is_even(mul(n, m)) == U);
    REQUIRE(is_even(mul(n, m1(m, 1))) == U);
    // Pre-existing parity arithmetic is unaffected.
    auto o = symbol("o", AssumptionMask{}.set_odd(true));
    REQUIRE(is_odd(mul(o, o)) == T);
    REQUIRE(is_even(mul(integer(2), n)) == T);
}
