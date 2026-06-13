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
#include <sympp/functions/trigonometric.hpp>

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

// Regression (SQRT-1): a perfect-square *rational* must reduce the same way a
// perfect-square integer does. sqrt(4)==2 already worked, but sqrt(1/4) used to
// stay as the unevaluated (1/4)**(1/2).
TEST_CASE("sqrt: perfect-square rationals fold to rationals",
          "[3d][sqrt][regression]") {
    // Reference: sympy.sqrt(Rational(1,4)) == Rational(1,2),
    //            sympy.sqrt(Rational(9,4)) == Rational(3,2).
    REQUIRE(sqrt(rational(1, 4)) == rational(1, 2));
    REQUIRE(sqrt(rational(9, 4)) == rational(3, 2));
    REQUIRE(sqrt(rational(9, 16)) == rational(3, 4));
    // Denominator becomes 1 → collapses to an Integer.
    REQUIRE(sqrt(rational(16, 1)) == integer(4));
}

// Regression (SQRT-2): square-factor extraction / denominator rationalisation.
// √8 → 2√2, √(1/2) → √2/2, √(2/3) → √6/3, √(8/9) → 2√2/3. A prime radicand
// (√7) and a square-free rational stay put in radical form (but rationalised).
TEST_CASE("sqrt: pulls out square factors and rationalises",
          "[3d][sqrt][regression][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(sqrt(integer(8))->str(), "2*sqrt(2)"));
    REQUIRE(oracle.equivalent(sqrt(integer(12))->str(), "2*sqrt(3)"));
    REQUIRE(oracle.equivalent(sqrt(integer(50))->str(), "5*sqrt(2)"));
    REQUIRE(oracle.equivalent(sqrt(rational(1, 2))->str(), "sqrt(2)/2"));
    REQUIRE(oracle.equivalent(sqrt(rational(8, 9))->str(), "2*sqrt(2)/3"));
    REQUIRE(oracle.equivalent(sqrt(rational(2, 3))->str(), "sqrt(6)/3"));
}

// Regression (NROOT-1): n-th-root factor extraction generalises SQRT-2 to any
// 1/n power. cbrt(16) → 2·cbrt(2), 96^(1/5) → 2·3^(1/5), 48^(1/4) → 2·3^(1/4),
// and the rational form (16/27)^(1/3) → 2·2^(1/3)/3. A prime radicand stays put.
TEST_CASE("nth-root: pulls out n-th-power factors and rationalises",
          "[3d][sqrt][nthroot][regression][oracle]") {
    auto& oracle = Oracle::instance();
    auto root = [](int b, int n) {
        return pow(integer(b), rational(1, n))->str();
    };
    REQUIRE(oracle.equivalent(root(16, 3), "2*2**Rational(1,3)"));
    REQUIRE(oracle.equivalent(root(54, 3), "3*2**Rational(1,3)"));
    REQUIRE(oracle.equivalent(root(24, 3), "2*3**Rational(1,3)"));
    REQUIRE(oracle.equivalent(root(96, 5), "2*3**Rational(1,5)"));
    REQUIRE(oracle.equivalent(root(48, 4), "2*3**Rational(1,4)"));
    REQUIRE(oracle.equivalent(
        pow(rational(16, 27), rational(1, 3))->str(), "2*2**Rational(1,3)/3"));
    REQUIRE(oracle.equivalent(
        pow(rational(2, 3), rational(1, 3))->str(), "18**Rational(1,3)/3"));
    // Prime radicand stays put (nothing pulls out).
    REQUIRE(oracle.equivalent(pow(integer(7), rational(1, 3))->str(),
                              "7**Rational(1,3)"));
}

// Regression (SQRT-3): principal square root of a negative number is
// imaginary — √(−a) = I·√a. SQRT-1/SQRT-2 deferred negative bases; this
// handles them for the ½ power, reusing the magnitude reduction.
TEST_CASE("sqrt: negative numbers fold to imaginary",
          "[3d][sqrt][regression][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(sqrt(integer(-1))->str(), "I"));
    REQUIRE(oracle.equivalent(sqrt(integer(-4))->str(), "2*I"));
    REQUIRE(oracle.equivalent(sqrt(integer(-8))->str(), "2*sqrt(2)*I"));
    REQUIRE(oracle.equivalent(sqrt(rational(-1, 4))->str(), "I/2"));
    REQUIRE(oracle.equivalent(sqrt(rational(-2, 3))->str(), "sqrt(6)*I/3"));
}

TEST_CASE("sqrt: prime radicand stays an irreducible Pow",
          "[3d][sqrt][regression]") {
    // Nothing to pull out of √7 — it must remain symbolic, not loop or mis-fold.
    auto e = sqrt(integer(7));
    REQUIRE(e->type_id() == TypeId::Pow);
    REQUIRE(e->str() == "7**(1/2)");
}

// Regression (SQRT-4): inverse square root rationalises the denominator.
// n^(-1/2) = sqrt(n)/n, with square-factor extraction (12^(-1/2) = sqrt(3)/6).
// Knock-on: atan(1/sqrt(3)) now folds to pi/6 (the argument is rationalised
// to the form the inverse-trig table recognises).
TEST_CASE("sqrt: inverse square root is rationalised",
          "[3d][sqrt][regression][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(pow(integer(3), rational(-1, 2))->str(),
                              "sqrt(3)/3"));
    REQUIRE(oracle.equivalent(pow(integer(2), rational(-1, 2))->str(),
                              "sqrt(2)/2"));
    REQUIRE(oracle.equivalent(pow(integer(12), rational(-1, 2))->str(),
                              "sqrt(3)/6"));
    REQUIRE(oracle.equivalent(pow(rational(2, 3), rational(-1, 2))->str(),
                              "sqrt(6)/2"));
    REQUIRE(oracle.equivalent(atan(pow(integer(3), rational(-1, 2)))->str(),
                              "pi/6"));
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

// Regression (ABS-1): Abs(c·rest) = |c|·Abs(rest) for a numeric coefficient.
// Previously only a leading −1 was pulled out; Abs(−2·x) etc. stayed put.
TEST_CASE("Abs: pulls out a numeric coefficient", "[3d][abs][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(oracle.equivalent(abs(mul(integer(-2), x))->str(), "2*Abs(x)"));
    REQUIRE(oracle.equivalent(abs(mul(integer(2), x))->str(), "2*Abs(x)"));
    REQUIRE(oracle.equivalent(abs(mul(rational(1, 2), x))->str(), "Abs(x)/2"));
    REQUIRE(oracle.equivalent(abs(mul(rational(-1, 3), x))->str(), "Abs(x)/3"));
    REQUIRE(oracle.equivalent(abs(mul(integer(-2), mul(x, y)))->str(),
                              "2*Abs(x*y)"));
}

TEST_CASE("Abs: pulls out positive/negative symbolic factors (ASSUME-5)",
          "[3d][abs][assumptions][regression]") {
    auto p = symbol("p", AssumptionMask{}.set_positive(true));
    auto q = symbol("q", AssumptionMask{}.set_positive(true));
    auto n = symbol("n", AssumptionMask{}.set_negative(true));
    auto x = symbol("x");
    auto y = symbol("y");
    // |p·x| = p·|x|; |p·q·x| = p·q·|x|; |n·x| = −n·|x| (|n| = −n for n<0).
    REQUIRE(abs(p * x) == p * abs(x));
    REQUIRE(abs(p * q * x) == p * q * abs(x));
    REQUIRE(abs(n * x) == mul(S::NegativeOne(), n) * abs(x));
    // A generic product is NOT split (modulus of x, y individually unknown).
    REQUIRE(abs(x * y)->str() == "Abs(x*y)");
    // Positive factor pulled, the rest kept under one Abs.
    REQUIRE(abs(p * x * y) == p * abs(x * y));
}

// Regression (ABS-2): Abs of a numeric complex number a+b·I → sqrt(a²+b²).
// Previously only real numbers reduced; a complex literal stayed symbolic.
TEST_CASE("Abs: complex modulus of a numeric a+b·I",
          "[3d][abs][oracle][regression]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(abs(integer(3) + integer(4) * S::I())->str(), "5"));
    REQUIRE(oracle.equivalent(abs(integer(1) + S::I())->str(), "sqrt(2)"));
    REQUIRE(oracle.equivalent(abs(integer(2) + integer(3) * S::I())->str(),
                              "sqrt(13)"));
    REQUIRE(oracle.equivalent(abs(integer(2) * S::I())->str(), "2"));
    REQUIRE(oracle.equivalent(abs(S::I())->str(), "1"));
    REQUIRE(oracle.equivalent(
        abs(integer(-3) - integer(4) * S::I())->str(), "5"));
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

// Regression (REIM-1): re/im/conjugate of a numeric complex a+b·I were left
// unevaluated (no Complex type). They now split the rational real/imaginary
// parts: re→a, im→b, conjugate→a−b·I.
TEST_CASE("re/im/conjugate: numeric complex a+b·I",
          "[3h][complex][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto z = integer(3) + integer(4) * S::I();
    REQUIRE(oracle.equivalent(re(z)->str(), "3"));
    REQUIRE(oracle.equivalent(im(z)->str(), "4"));
    REQUIRE(oracle.equivalent(conjugate(z)->str(), "3 - 4*I"));
    // Pure imaginary.
    REQUIRE(oracle.equivalent(re(integer(2) * S::I())->str(), "0"));
    REQUIRE(oracle.equivalent(im(integer(2) * S::I())->str(), "2"));
    REQUIRE(oracle.equivalent(conjugate(integer(2) * S::I())->str(), "-2*I"));
    // Rational parts.
    auto w = rational(1, 2) + rational(1, 3) * S::I();
    REQUIRE(oracle.equivalent(re(w)->str(), "1/2"));
    REQUIRE(oracle.equivalent(conjugate(w)->str(), "1/2 - I/3"));
}

// CONJ-DIST-1: conjugate is a ring homomorphism — distributes over products,
// sums and integer powers (conjugate(I·x) = −I·conjugate(x)). Matches SymPy.
TEST_CASE("conjugate: distributes over Mul/Add/Pow (CONJ-DIST-1)",
          "[3h][conjugate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(oracle.equivalent(conjugate(S::I() * x)->str(),
                              "-I*conjugate(x)"));
    REQUIRE(oracle.equivalent(conjugate(integer(2) * x)->str(),
                              "2*conjugate(x)"));
    REQUIRE(oracle.equivalent(conjugate(x * y)->str(),
                              "conjugate(x)*conjugate(y)"));
    REQUIRE(oracle.equivalent(conjugate(x + y)->str(),
                              "conjugate(x) + conjugate(y)"));
    REQUIRE(oracle.equivalent(conjugate(x + S::I())->str(),
                              "conjugate(x) - I"));
    REQUIRE(oracle.equivalent(conjugate(pow(x, integer(2)))->str(),
                              "conjugate(x)**2"));
}

// ABS-I-1: |I·x| = |x| (I has modulus 1) — the imaginary unit is pulled out of
// an Abs product alongside the numeric / known-sign factors.
TEST_CASE("Abs: pulls out the imaginary unit (ABS-I-1)",
          "[3d][abs][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(oracle.equivalent(abs(S::I() * x)->str(), "Abs(x)"));
    REQUIRE(oracle.equivalent(abs(integer(2) * S::I() * x)->str(), "2*Abs(x)"));
    REQUIRE(oracle.equivalent(abs(S::I() * x * y)->str(), "Abs(x*y)"));
}

// ABS-MOD-1: symbolic complex modulus |a + b·I| = sqrt(a² + b²) when the real
// and imaginary parts resolve (real symbols); a generic Abs(z) stays put.
TEST_CASE("Abs: symbolic complex modulus (ABS-MOD-1)",
          "[3d][abs][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    auto y = symbol("y", AssumptionMask{}.set_real(true));
    REQUIRE(oracle.equivalent(abs(x + S::I() * y)->str(), "sqrt(x**2 + y**2)"));
    REQUIRE(oracle.equivalent(abs(integer(2) + S::I() * y)->str(),
                              "sqrt(y**2 + 4)"));
    // A generic (unknown-reality) symbol keeps its Abs.
    auto z = symbol("z");
    REQUIRE(abs(z)->str() == "Abs(z)");
    // A real symbol keeps its Abs (no imaginary part).
    REQUIRE(abs(x)->str() == "Abs(x)");
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

// REIM-DIST-1: re/im are linear over sums and rotate by I —
//   re(I·w) = −im(w), im(I·w) = re(w), with a real coefficient pulled out.
// For real x, y: re(I·x)=0, im(I·x)=x, re(x+I·y)=x, im(x+I·y)=y; for generic z,
// re(I·z) = −im(z). Matches SymPy.
TEST_CASE("re/im: linearity and I-rotation (REIM-DIST-1)",
          "[3h][complex][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    auto y = symbol("y", AssumptionMask{}.set_real(true));
    auto z = symbol("z");  // generic
    REQUIRE(re(S::I() * x) == S::Zero());
    REQUIRE(im(S::I() * x) == x);
    REQUIRE(re(x + S::I() * y) == x);
    REQUIRE(im(x + S::I() * y) == y);
    REQUIRE(re(integer(2) * x) == integer(2) * x);
    // Generic z: the I-rotation surfaces re/im of z.
    REQUIRE(oracle.equivalent(re(S::I() * z)->str(), "-im(z)"));
    REQUIRE(oracle.equivalent(im(S::I() * z)->str(), "re(z)"));
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
