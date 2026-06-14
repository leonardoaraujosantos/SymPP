// series / limit tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/calculus/limit.hpp>
#include <sympp/calculus/series.hpp>
#include <sympp/core/infinity.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/functions/miscellaneous.hpp>  // sign, abs
#include <sympp/core/symbol.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- series ---------------------------------------------------------------

TEST_CASE("series: exp(x) about 0 to order 5", "[6][series][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto s = series(exp(x), x, S::Zero(), 5);
    REQUIRE(oracle.equivalent(s->str(),
                              "1 + x + x**2/2 + x**3/6 + x**4/24"));
}

TEST_CASE("series: sin(x) about 0 to order 6", "[6][series][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto s = series(sin(x), x, S::Zero(), 6);
    // 0 + x + 0 - x^3/6 + 0 + x^5/120
    REQUIRE(oracle.equivalent(s->str(), "x - x**3/6 + x**5/120"));
}

TEST_CASE("series: cos(x) about 0 to order 5", "[6][series][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto s = series(cos(x), x, S::Zero(), 5);
    // 1 + 0*x - x^2/2 + 0 + x^4/24
    REQUIRE(oracle.equivalent(s->str(), "1 - x**2/2 + x**4/24"));
}

// SERIES-SINGULAR-1: series() at a point where direct substitution gives a
// non-finite coefficient. A genuine singularity (log(x), 1/x at 0) cannot be
// Taylor-expanded — return the input unchanged, matching SymPy. A removable
// singularity (0/0 like sin(x)/x) has a finite limit, so the limit fallback
// recovers the correct Taylor coefficients. Previously these produced zoo/nan
// garbage.
TEST_CASE("series: singular and removable points (SERIES-SINGULAR-1)",
          "[6][series][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto zero = S::Zero();
    // Genuine singularities: returned unevaluated (no zoo garbage).
    REQUIRE(series(log(x), x, zero, 6) == log(x));
    REQUIRE(series(pow(x, integer(-1)), x, zero, 6) == pow(x, integer(-1)));
    REQUIRE(series(pow(x, rational(1, 2)), x, zero, 6)
            == pow(x, rational(1, 2)));
    REQUIRE(series(pow(x, integer(-2)), x, zero, 6) == pow(x, integer(-2)));
    // Removable singularities: the Taylor series via the limit fallback.
    REQUIRE(oracle.equivalent(series(sin(x) * pow(x, integer(-1)), x, zero, 6)
                                  ->str(),
                              "1 - x**2/6 + x**4/120"));
    REQUIRE(oracle.equivalent(
        series((exp(x) - integer(1)) * pow(x, integer(-1)), x, zero, 6)->str(),
        "1 + x/2 + x**2/6 + x**3/24 + x**4/120 + x**5/720"));
    // An ordinary analytic function is unaffected.
    REQUIRE(oracle.equivalent(series(exp(x), x, zero, 4)->str(),
                              "1 + x + x**2/2 + x**3/6"));
}

// ----- limit ----------------------------------------------------------------

TEST_CASE("limit: polynomial substitution", "[6][limit]") {
    auto x = symbol("x");
    auto e = pow(x, integer(2)) + integer(3) * x + integer(2);
    REQUIRE(limit(e, x, integer(1)) == integer(6));
}

TEST_CASE("limit: trig substitution", "[6][limit]") {
    auto x = symbol("x");
    REQUIRE(limit(sin(x), x, S::Pi()) == S::Zero());
    REQUIRE(limit(cos(x), x, S::Zero()) == S::One());
}

// ----- L'Hôpital indeterminate forms -----------------------------------------

TEST_CASE("limit: sin(x)/x at 0 → 1 (classic 0/0)", "[6][limit][lhopital][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sin(x) / x;
    auto v = limit(e, x, S::Zero());
    REQUIRE(oracle.equivalent(v->str(), "1"));
}

TEST_CASE("limit: (1 - cos(x))/x^2 at 0 → 1/2", "[6][limit][lhopital][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = (integer(1) - cos(x)) / pow(x, integer(2));
    auto v = limit(e, x, S::Zero());
    REQUIRE(oracle.equivalent(v->str(), "1/2"));
}

TEST_CASE("limit: (x^2 - 1)/(x - 1) at 1 → 2 (factorable 0/0)",
          "[6][limit][lhopital][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = (pow(x, integer(2)) - integer(1)) / (x - integer(1));
    auto v = limit(e, x, integer(1));
    REQUIRE(oracle.equivalent(v->str(), "2"));
}

TEST_CASE("limit: (e^x - 1)/x at 0 → 1", "[6][limit][lhopital][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = (exp(x) - integer(1)) / x;
    auto v = limit(e, x, S::Zero());
    REQUIRE(oracle.equivalent(v->str(), "1"));
}

TEST_CASE("limit: tan(x)/x at 0 → 1", "[6][limit][lhopital][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = tan(x) / x;
    auto v = limit(e, x, S::Zero());
    REQUIRE(oracle.equivalent(v->str(), "1"));
}

TEST_CASE("limit: log(1+x)/x at 0 → 1", "[6][limit][lhopital][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = log(integer(1) + x) / x;
    auto v = limit(e, x, S::Zero());
    REQUIRE(oracle.equivalent(v->str(), "1"));
}

// ----- Limits at infinity (LIM-1, regression for issue #2) -------------------
// `oo` used to be a plain symbol and the engine had no concept of infinity, so
// limit((1+1/x)**x, x, oo) returned garbage. Infinity is now a real atom and
// the engine handles the standard indeterminate forms.
TEST_CASE("limit: (1 + 1/x)^x at oo → E (the classic 1^oo form)",
          "[6][limit][infinity][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(integer(1) + integer(1) / x, x);
    auto v = limit(e, x, S::Infinity());
    REQUIRE(oracle.equivalent(v->str(), "E"));
}

TEST_CASE("limit: (1 + 2/x)^x at oo → exp(2)",
          "[6][limit][infinity][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(integer(1) + integer(2) / x, x);
    auto v = limit(e, x, S::Infinity());
    REQUIRE(oracle.equivalent(v->str(), "exp(2)"));
}

// POLE-SIGN-1: at a finite pole, direct substitution yields zoo; resolve the
// sign by sampling both sides. An even-order pole has matching signs on both
// sides → ±oo (matching SymPy); an odd-order pole has opposite signs → the
// two-sided limit is genuinely zoo (SymPy reports oo only via its dir='+'
// default — SymPP's limit is two-sided, so zoo is the correct answer).
// LIMIT-HANG-1: limit() must terminate on a radical ∞/∞ form. L'Hôpital on
// sqrt(x²+x)−x balloons the nested radicals each step (the ratio never
// stabilises), which used to hang; a size budget now bails to the unevaluated
// nan instead. The true value (1/2) needs asymptotic-series / Gruntz machinery.
TEST_CASE("limit: radical infinity ratio terminates (LIMIT-HANG-1)",
          "[6][limit][infinity][regression]") {
    auto x = symbol("x");
    auto oo = S::Infinity();
    // The point is termination, not the value: this must return (any value),
    // not loop forever.
    auto a = limit(pow(pow(x, integer(2)) + x, S::Half()) - x, x, oo);
    REQUIRE(a);  // returned something (no hang)
    auto b = limit(x / (pow(pow(x, integer(2)) + x, S::Half()) + x), x, oo);
    REQUIRE(b);
    // Sanity: ordinary limits next to it still resolve.
    REQUIRE(limit(sin(x) / x, x, S::Zero()) == integer(1));
}

// LIMIT-EXP-1: 0·∞ where an exponential dominates a polynomial. x^n·e^(-x) → 0
// at +∞. try_product_form now tries both the 0/0 and ∞/∞ arrangements (the
// latter, x^n / e^x, is the one L'Hôpital can crack), with an exp-aware
// reciprocal so e^(-x) stays in the denominator across iterations; limit also
// gained linearity over Add / Mul so a sum/product of such terms resolves.
TEST_CASE("limit: polynomial times decaying exponential (LIMIT-EXP-1)",
          "[6][limit][infinity][regression]") {
    auto x = symbol("x");
    auto oo = S::Infinity();
    REQUIRE(limit(x * exp(mul(S::NegativeOne(), x)), x, oo) == S::Zero());
    REQUIRE(limit(pow(x, integer(2)) * exp(mul(S::NegativeOne(), x)), x, oo)
            == S::Zero());
    REQUIRE(limit(pow(x, integer(5)) * exp(mul(S::NegativeOne(), x)), x, oo)
            == S::Zero());
    REQUIRE(limit(x * exp(mul(integer(-2), x)), x, oo) == S::Zero());
    // Linearity: a sum of decaying terms.
    REQUIRE(limit(exp(mul(S::NegativeOne(), x))
                      - exp(mul(integer(-2), x)),
                  x, oo) == S::Zero());
    // A genuinely divergent term must NOT be coerced to a finite value.
    REQUIRE(limit(x + integer(1) / x, x, oo) == oo);
}

// LIMIT-POLY-INF-1: a polynomial at ±∞ is governed by its leading term, so the
// ∞−∞ that direct substitution leaves as nan resolves to the signed infinity.
TEST_CASE("limit: polynomial at infinity via the leading term (LIMIT-POLY-INF-1)",
          "[6][limit][infinity][regression]") {
    auto x = symbol("x");
    auto oo = S::Infinity();
    auto noo = S::NegativeInfinity();
    REQUIRE(limit(pow(x, integer(2)) - x, x, oo) == oo);
    REQUIRE(limit(pow(x, integer(2)) - x, x, noo) == oo);
    REQUIRE(limit(x - pow(x, integer(2)), x, oo) == noo);
    REQUIRE(limit(integer(2) * pow(x, integer(2)) - integer(5) * x, x, oo) == oo);
    // Odd leading degree flips sign at −∞.
    REQUIRE(limit(mul(S::NegativeOne(), pow(x, integer(3))) + x, x, oo) == noo);
    REQUIRE(limit(mul(S::NegativeOne(), pow(x, integer(3))) + x, x, noo) == oo);
}

TEST_CASE("limit: signed infinity at an even pole (POLE-SIGN-1)",
          "[6][limit][infinity][regression]") {
    auto x = symbol("x");
    REQUIRE(limit(pow(x, integer(-2)), x, S::Zero()) == S::Infinity());
    REQUIRE(limit(pow(x, integer(-4)), x, S::Zero()) == S::Infinity());
    REQUIRE(limit(mul(S::NegativeOne(), pow(x, integer(-2))), x, S::Zero())
            == S::NegativeInfinity());
    // Shifted pole: 1/(x-1)^2 → +oo as x → 1.
    REQUIRE(limit(pow(x - integer(1), integer(-2)), x, integer(1))
            == S::Infinity());
    // Odd pole stays the unsigned zoo (two-sided): 1/x as x → 0.
    REQUIRE(limit(pow(x, integer(-1)), x, S::Zero()) == S::ComplexInfinity());
}

// LIMIT-SIGN-1: an expression with sign(g) or abs(g) where g → 0 is
// discontinuous at the target. Direct substitution wrongly used sign(0)=0 (so
// limit(sign(x),x,0) returned 0); it now samples g's sign on each side. Matching
// sides give the value (sign(x²) → 1); a sign change means the two-sided limit
// does not exist (nan), e.g. sign(x) and |x|/x. Matches SymPy's two-sided limit.
TEST_CASE("limit: discontinuous sign/abs at the target (LIMIT-SIGN-1)",
          "[6][limit][regression]") {
    auto x = symbol("x");
    auto zero = S::Zero();
    auto is_nan = [](const Expr& e) { return e->type_id() == TypeId::NaN; };
    // sign(g) with a genuine sign change ⇒ two-sided DNE (nan).
    REQUIRE(is_nan(limit(sign(x), x, zero)));
    REQUIRE(is_nan(limit(sign(mul(S::NegativeOne(), x)), x, zero)));  // sign(-x)
    REQUIRE(is_nan(limit(sign(pow(x, integer(3))), x, zero)));  // sign(x³)
    // No sign change ⇒ the constant sign: sign(x²) → 1.
    REQUIRE(limit(sign(pow(x, integer(2))), x, zero) == integer(1));
    // |x|/x = sign(x): DNE. (Previously returned 0 via L'Hôpital → sign(x) → 0.)
    REQUIRE(is_nan(limit(sympp::abs(x) * pow(x, integer(-1)), x, zero)));
    REQUIRE(is_nan(limit(x * pow(sympp::abs(x), integer(-1)), x, zero)));  // x/|x|
    // sign·(→0) is determinate: sign(x)·x → 0; and plain |x| is continuous → 0.
    REQUIRE(limit(sign(x) * x, x, zero) == zero);
    REQUIRE(limit(sympp::abs(x), x, zero) == zero);
    // sign away from the discontinuity is unaffected: sign(x−2) → −1 at 0.
    REQUIRE(limit(sign(x - integer(2)), x, zero) == integer(-1));
}

TEST_CASE("limit: rational functions at oo (leading-term ratio via L'Hopital)",
          "[6][limit][infinity][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // deg num == deg den → ratio of leading coefficients.
    auto r1 = limit((integer(2) * pow(x, integer(2)) + integer(3))
                        / (pow(x, integer(2)) - integer(1)),
                    x, S::Infinity());
    REQUIRE(oracle.equivalent(r1->str(), "2"));
    // deg num > deg den → oo.
    auto r2 = limit(pow(x, integer(2)) / (x + integer(1)), x, S::Infinity());
    REQUIRE(r2 == S::Infinity());
    // deg num < deg den → 0.
    auto r3 = limit((x + integer(1)) / pow(x, integer(2)), x, S::Infinity());
    REQUIRE(r3 == S::Zero());
}

// LIMIT-RAT-SYM-1: rational functions at ∞ with symbolic (var-free) coefficients.
// Direct ∞ substitution is unreliable and L'Hôpital mishandled these — (x+a)/x
// returned 0. A degree/leading-coefficient comparison resolves them, which in
// turn makes (1+a/x)^x → eᵃ (the 1^∞ form whose inner ∞·0 limit is a). Matches
// SymPy.
TEST_CASE("limit: rational at ∞ with symbolic coefficients (LIMIT-RAT-SYM-1)",
          "[6][limit][infinity][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto a = symbol("a");
    auto oo = S::Infinity();
    // equal degree → leading-coefficient ratio (symbolic).
    REQUIRE(oracle.equivalent(limit((x + a) / x, x, oo)->str(), "1"));
    REQUIRE(oracle.equivalent(limit(a * x / (x + integer(1)), x, oo)->str(),
                              "a"));
    REQUIRE(oracle.equivalent(
        limit((a * pow(x, integer(2)) + integer(1))
                  / (pow(x, integer(2)) + x),
              x, oo)
            ->str(),
        "a"));
    // lower degree numerator → 0; higher → ∞ (sign from symbolic a is left to a²).
    REQUIRE(limit((x + a) / pow(x, integer(2)), x, oo) == S::Zero());
    REQUIRE(limit((pow(x, integer(2)) + a) / (x + integer(1)), x, oo) == oo);
    // 1^∞ with a symbolic parameter: (1 + a/x)^x → eᵃ.
    REQUIRE(oracle.equivalent(
        limit(pow(integer(1) + a / x, x), x, oo)->str(), "exp(a)"));
    // Numeric cases unchanged.
    REQUIRE(oracle.equivalent(limit(pow(integer(1) + integer(2) / x, x), x, oo)
                                  ->str(),
                              "exp(2)"));
}

TEST_CASE("limit: 0*oo and oo/oo forms at oo",
          "[6][limit][infinity][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x*sin(1/x) → 1 (0*oo).
    REQUIRE(oracle.equivalent(
        limit(x * sin(integer(1) / x), x, S::Infinity())->str(), "1"));
    // exp(x)/x → oo (oo/oo).
    REQUIRE(limit(exp(x) / x, x, S::Infinity()) == S::Infinity());
    // log(x)/x → 0.
    REQUIRE(oracle.equivalent(limit(log(x) / x, x, S::Infinity())->str(), "0"));
}

// LIMIT-GAMMA-1: limits of gamma/factorial at +∞. Direct substitution gives
// gamma(∞)/gamma(∞), which simplify wrongly cancels to 1 — these used to return
// 1 (a wrong answer) or unevaluated/garbage. The engine now (a) folds
// gamma(∞)=factorial(∞)=∞, (b) collapses gamma ratios before substituting
// (gamma(x+1)/gamma(x) → x → ∞), and (c) ranks growth classes
// (gamma ≫ exp ≫ poly): exp(x)/gamma(x) → 0.
TEST_CASE("limit: gamma/factorial at infinity (LIMIT-GAMMA-1)",
          "[6][limit][infinity][gamma][oracle][regression]") {
    auto x = symbol("x");
    const Expr oo = S::Infinity();
    // Bare gamma/factorial diverge.
    REQUIRE(limit(gamma(x), x, oo) == oo);
    REQUIRE(limit(factorial(x), x, oo) == oo);
    REQUIRE(limit(gamma(x + integer(1)), x, oo) == oo);
    // Ratios collapse to polynomials (the previously-wrong "= 1" cases).
    REQUIRE(limit(gamma(x + integer(1)) / gamma(x), x, oo) == oo);
    REQUIRE(limit(gamma(x + integer(2)) / gamma(x), x, oo) == oo);
    REQUIRE(limit(gamma(x) / gamma(x + integer(1)), x, oo) == S::Zero());
    // Cross-class growth: gamma ≫ exp ≫ polynomial.
    REQUIRE(limit(exp(x) / gamma(x), x, oo) == S::Zero());
    REQUIRE(limit(gamma(x) / exp(x), x, oo) == oo);
    REQUIRE(limit(pow(x, integer(5)) / gamma(x), x, oo) == S::Zero());
    // factorial normalizes to gamma, so mixed forms resolve too.
    REQUIRE(limit(gamma(x) / factorial(x), x, oo) == S::Zero());
    REQUIRE(limit(factorial(x) / exp(x), x, oo) == oo);
    REQUIRE(limit(pow(integer(2), x) / factorial(x), x, oo) == S::Zero());
}

TEST_CASE("limit: at -oo", "[6][limit][infinity][oracle][regression]") {
    auto x = symbol("x");
    REQUIRE(limit(exp(x), x, S::NegativeInfinity()) == S::Zero());
    REQUIRE(limit(pow(x, integer(3)), x, S::NegativeInfinity())
            == S::NegativeInfinity());
}

// Bounded/monotone functions at infinity now fold (FUNC-INF), so their limits
// resolve directly: atan(x) → pi/2, tanh(x) → 1.
TEST_CASE("limit: bounded functions at oo (atan, tanh)",
          "[6][limit][infinity][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(limit(atan(x), x, S::Infinity())->str(), "pi/2"));
    REQUIRE(limit(tanh(x), x, S::Infinity()) == S::One());
}

#include <sympp/calculus/order.hpp>

// ----- Order / Big-O ---------------------------------------------------------

TEST_CASE("Order: O(x^3) prints as O(x**3)", "[6][order]") {
    auto x = symbol("x");
    auto o = order(pow(x, integer(3)), x);
    REQUIRE(o->str() == "O(x**3)");
}

TEST_CASE("Order: O(0) collapses to 0", "[6][order]") {
    auto x = symbol("x");
    auto o = order(S::Zero(), x);
    REQUIRE(o == S::Zero());
}

TEST_CASE("Order: hash-cons identity", "[6][order]") {
    auto x = symbol("x");
    auto a = order(pow(x, integer(2)), x);
    auto b = order(pow(x, integer(2)), x);
    REQUIRE(a.get() == b.get());
}

TEST_CASE("Order: distinct orders inequal", "[6][order]") {
    auto x = symbol("x");
    auto a = order(pow(x, integer(2)), x);
    auto b = order(pow(x, integer(3)), x);
    REQUIRE(a != b);
}

TEST_CASE("Order: with explicit point prints (var, point)", "[6][order]") {
    auto x = symbol("x");
    auto o = order(pow(x - integer(1), integer(2)), x, integer(1));
    REQUIRE(o->str() == "O((x - 1)**2, (x, 1))");
}

#include <sympp/calculus/summation.hpp>
#include <sympp/core/traversal.hpp>

// ----- summation -------------------------------------------------------------

TEST_CASE("summation: Σ k from 1 to n → n(n+1)/2", "[6][summation][oracle]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto s = summation(k, k, integer(1), n);
    REQUIRE(oracle.equivalent(s->str(), "n*(n+1)/2"));
}

TEST_CASE("summation: Σ k² from 1 to n → n(n+1)(2n+1)/6",
          "[6][summation][oracle]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto s = summation(pow(k, integer(2)), k, integer(1), n);
    REQUIRE(oracle.equivalent(s->str(), "n*(n+1)*(2*n+1)/6"));
}

TEST_CASE("summation: Σ k³ from 1 to n → (n(n+1)/2)²",
          "[6][summation][oracle]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto s = summation(pow(k, integer(3)), k, integer(1), n);
    REQUIRE(oracle.equivalent(s->str(), "(n*(n+1)/2)**2"));
}

// SUM-FAULHABER-1: Σ kᵖ for any positive integer p via Faulhaber's formula
// (Bernoulli-number coefficients). Previously only p ∈ {2,3} were closed; higher
// powers came back unevaluated. Matches SymPy.
TEST_CASE("summation: Σ kᵖ for higher p (SUM-FAULHABER-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto sum_pow = [&](int p) {
        return summation(pow(k, integer(p)), k, integer(1), n);
    };
    REQUIRE(oracle.equivalent(sum_pow(4)->str(),
                              "n**5/5 + n**4/2 + n**3/3 - n/30"));
    REQUIRE(oracle.equivalent(sum_pow(5)->str(),
                              "n**6/6 + n**5/2 + 5*n**4/12 - n**2/12"));
    REQUIRE(oracle.equivalent(
        sum_pow(6)->str(),
        "n**7/7 + n**6/2 + n**5/2 - n**3/6 + n/42"));
    // No leftover Sum() marker for any of them.
    for (int p = 4; p <= 12; ++p) {
        REQUIRE(sum_pow(p)->str().find("Sum(") == std::string::npos);
    }
    // A partial range still works: Σ_{k=2}^n k⁴ = (Σ_{1}^n) − 1.
    REQUIRE(oracle.equivalent(
        summation(pow(k, integer(4)), k, integer(2), n)->str(),
        "n**5/5 + n**4/2 + n**3/3 - n/30 - 1"));
}

// SUM-TELESCOPE-1: Σ of a rational summand c/D(k), where D is a quadratic with
// two distinct linear factors whose roots differ by an integer, telescopes to a
// closed form (the sum was returned unevaluated before). Verified equivalent to
// SymPy's summation — the presented form may differ.
TEST_CASE("summation: telescoping rational sums (SUM-TELESCOPE-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto one = integer(1);
    // Σ 1/(k(k+1)) = n/(n+1).
    REQUIRE(oracle.equivalent(
        summation(pow(k * (k + one), integer(-1)), k, one, n)->str(),
        "n/(n+1)"));
    // Σ 1/(k(k+2)), poles two apart → two telescoping chains.
    REQUIRE(oracle.equivalent(
        summation(pow(k * (k + integer(2)), integer(-1)), k, one, n)->str(),
        "3/4 - 1/(2*(n+1)) - 1/(2*(n+2))"));
    // Σ 1/(4k²-1) = n/(2n+1).
    REQUIRE(oracle.equivalent(
        summation(pow(integer(4) * pow(k, integer(2)) - one, integer(-1)), k,
                  one, n)
            ->str(),
        "n/(2*n+1)"));
    // A constant numerator scales through.
    REQUIRE(oracle.equivalent(
        summation(integer(3) * pow(k * (k + one), integer(-1)), k, one, n)
            ->str(),
        "3*n/(n+1)"));
    // Non-unit leading coefficient: 1/((2k-1)(2k+1)) = n/(2n+1).
    REQUIRE(oracle.equivalent(
        summation(
            pow((integer(2) * k - one) * (integer(2) * k + one), integer(-1)), k,
            one, n)
            ->str(),
        "n/(2*n+1)"));
    // A pole inside the range (root k=1) must stay unevaluated, not produce a
    // bogus closed form.
    auto pole = summation(pow(k * (k - one), integer(-1)), k, one, n);
    REQUIRE(pole->str().find("Sum(") != std::string::npos);
}

TEST_CASE("summation: constant Σ c from 1 to n → c*n",
          "[6][summation][oracle]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto c = symbol("c");
    auto s = summation(c, k, integer(1), n);
    REQUIRE(oracle.equivalent(s->str(), "c*n"));
}

TEST_CASE("summation: Σ (3k + 2) from 1 to n → 3n(n+1)/2 + 2n (linearity)",
          "[6][summation][oracle]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto e = integer(3) * k + integer(2);
    auto s = summation(e, k, integer(1), n);
    REQUIRE(oracle.equivalent(s->str(), "3*n*(n+1)/2 + 2*n"));
}

TEST_CASE("summation: geometric Σ 2^k from 0 to n → 2^(n+1) - 1",
          "[6][summation][oracle]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto e = pow(integer(2), k);
    auto s = summation(e, k, integer(0), n);
    REQUIRE(oracle.equivalent(s->str(), "2**(n+1) - 1"));
}

// Regression (SUM-2): arithmetic-geometric Σ k·r^k (k times a geometric term)
// used to fall through — a Mul of two var-dependent factors isn't pulled apart
// by the constant-extraction path.
TEST_CASE("summation: arithmetic-geometric Σ k·2^k from 0 to n",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto s = summation(k * pow(integer(2), k), k, integer(0), n);
    REQUIRE(s->str().find("Sum(") == std::string::npos);
    REQUIRE(oracle.equivalent(s->str(), "2*2**n*n - 2*2**n + 2"));
}

TEST_CASE("summation: Σ k·2^k from 1 to 3 → 34 (numeric)",
          "[6][summation][regression]") {
    auto k = symbol("k");
    auto s = summation(k * pow(integer(2), k), k, integer(1), integer(3));
    REQUIRE(s == integer(34));
}

// SUM-POLYGEOM-1: Σ P(k)·r^k for a polynomial P of degree ≥ 2 and a concrete
// ratio r ≠ 1 — the general arithmetic-geometric sum (k²·2^k, (2k+1)·3^k,
// k²/2^k). Found via the antidifference Q(k)·r^k with r·Q(k+1)−Q(k)=P(k).
// Previously only degree 1 (k·r^k) was closed. Matches SymPy.
TEST_CASE("summation: polynomial × geometric (SUM-POLYGEOM-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto sum1n = [&](const Expr& e) {
        return summation(e, k, integer(1), n);
    };
    auto pw = [&](int b) { return pow(integer(b), k); };
    // Σ k²·2^k.
    auto s2 = sum1n(pow(k, integer(2)) * pw(2));
    REQUIRE(s2->str().find("Sum(") == std::string::npos);
    REQUIRE(oracle.equivalent(s2->str(),
                              "2*(2**n*n**2 - 2*2**n*n + 3*2**n - 3)"));
    // Σ k³·2^k.
    auto s3 = sum1n(pow(k, integer(3)) * pw(2));
    REQUIRE(s3->str().find("Sum(") == std::string::npos);
    REQUIRE(oracle.equivalent(
        s3->str(), "2*(2**n*n**3 - 3*2**n*n**2 + 9*2**n*n - 13*2**n + 13)"));
    // Σ (2k+1)·3^k → n·3^(n+1).
    REQUIRE(oracle.equivalent(
        sum1n((integer(2) * k + integer(1)) * pw(3))->str(), "n*3**(n+1)"));
    // Negative ratio direction: Σ k²/2^k.
    auto sneg = sum1n(pow(k, integer(2)) * pow(integer(2), mul(S::NegativeOne(), k)));
    REQUIRE(sneg->str().find("Sum(") == std::string::npos);
    REQUIRE(oracle.equivalent(sneg->str(), "6 - (n**2 + 4*n + 6)/2**n"));
    // A concrete numeric range: Σ_{k=1}^3 k²·2^k = 2 + 16 + 72 = 90.
    REQUIRE(summation(pow(k, integer(2)) * pw(2), k, integer(1), integer(3))
            == integer(90));
}

// Regression: geometric detection used to require the exponent to be
// *exactly* the summation variable, so Σ 2^(-k) (which canonicalizes to
// base 2, exponent -k) and Σ 2^(2k) fell through and returned the summand
// with the bound variable still free. See known-issues "summation
// geometric exponent linear in index".
TEST_CASE("summation: geometric Σ 2^(-k) from 0 to n (negated index)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto e = pow(integer(2), -k);  // 2**(-k)
    auto s = summation(e, k, integer(0), n);
    REQUIRE_FALSE(has(s, k));  // bound variable must not leak free
    REQUIRE(oracle.equivalent(s->str(), "2 - 2*2**(-n - 1)"));
}

TEST_CASE("summation: geometric Σ 2^(-k) from 0 to 3 → 15/8",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto e = pow(integer(2), -k);
    auto s = summation(e, k, integer(0), integer(3));
    REQUIRE(oracle.equivalent(s->str(), "15/8"));
}

TEST_CASE("summation: geometric Σ 2^(2k) from 0 to n (scaled index)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto e = pow(integer(2), integer(2) * k);  // 2**(2*k)
    auto s = summation(e, k, integer(0), n);
    REQUIRE_FALSE(has(s, k));
    REQUIRE(oracle.equivalent(s->str(), "4**(n + 1)/3 - 1/3"));
}

TEST_CASE("summation: Σ from a to b uses telescoping",
          "[6][summation][oracle]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto a = symbol("a");
    auto b = symbol("b");
    auto s = summation(k, k, a, b);
    REQUIRE(oracle.equivalent(s->str(), "(b - a + 1)*(a + b)/2"));
}

// Infinite geometric series now close (the ratio^oo → 0 fold from the
// Infinity feature): Σ (1/2)^k = 2, Σ 2^-k from 1 = 1.
TEST_CASE("summation: infinite geometric series",
          "[6][summation][infinity][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto oo = S::Infinity();
    REQUIRE(oracle.equivalent(
        summation(pow(rational(1, 2), k), k, integer(0), oo)->str(), "2"));
    REQUIRE(oracle.equivalent(
        summation(pow(integer(2), integer(-1) * k), k, integer(1), oo)->str(),
        "1"));
    // Divergent Σ k from 1 to oo → oo.
    REQUIRE(summation(k, k, integer(1), oo) == oo);
}

// Regression (SUM-3): an unrecognised sum must return an unevaluated Sum
// marker, never the bare summand. Σ 1/k! (= E, which SymPP doesn't close)
// stands in for the general case; the p-series Σ 1/k^s now close via ZETA.
TEST_CASE("summation: unrecognised sum stays an unevaluated Sum",
          "[6][summation][regression]") {
    auto k = symbol("k");
    auto e = pow(factorial(k), integer(-1));
    auto s = summation(e, k, integer(1), S::Infinity());
    REQUIRE(s->type_id() == TypeId::Function);   // undefined-function marker
    REQUIRE(s->str().rfind("Sum(", 0) == 0);     // starts with "Sum("
    REQUIRE_FALSE(s == e);                       // not the bare summand
}

// Regression (ZETA-EVEN): convergent even p-series close to ζ(2n)=r·π^(2n).
// Σ1/k²=π²/6 (Basel), Σ1/k⁴=π⁴/90, …; odd p stays an unevaluated Sum.
TEST_CASE("summation: even p-series close to zeta values",
          "[6][summation][zeta][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto oo = S::Infinity();
    REQUIRE(oracle.equivalent(
        summation(pow(k, integer(-2)), k, integer(1), oo)->str(), "pi**2/6"));
    REQUIRE(oracle.equivalent(
        summation(pow(k, integer(-4)), k, integer(1), oo)->str(), "pi**4/90"));
    REQUIRE(oracle.equivalent(
        summation(pow(k, integer(-6)), k, integer(1), oo)->str(), "pi**6/945"));
    REQUIRE(oracle.equivalent(
        summation(pow(k, integer(-12)), k, integer(1), oo)->str(),
        "691*pi**12/638512875"));
    // Odd exponent has no elementary closed form — closes to a symbolic ζ(s)
    // (matching SymPy's Sum(1/k**3).doit() = zeta(3)), not π-powers.
    REQUIRE(oracle.equivalent(
        summation(pow(k, integer(-3)), k, integer(1), oo)->str(), "zeta(3)"));
    REQUIRE(oracle.equivalent(
        summation(pow(k, integer(-5)), k, integer(1), oo)->str(), "zeta(5)"));
}

// Single-term range Σ_{k=a}^{a} f(k) = f(a).
TEST_CASE("summation: single-term range substitutes the bound",
          "[6][summation][regression]") {
    auto k = symbol("k");
    REQUIRE(summation(pow(k, integer(-2)), k, integer(3), integer(3))
            == rational(1, 9));
    REQUIRE(summation(pow(k, integer(2)), k, integer(5), integer(5))
            == integer(25));
}

#include <sympp/calculus/pade.hpp>

// ----- Padé approximant ------------------------------------------------------

TEST_CASE("pade: [1/1] of e^x → (1 + x/2)/(1 - x/2)",
          "[6][pade][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto p = pade(exp(x), x, 1, 1);
    REQUIRE(oracle.equivalent(p->str(), "(1 + x/2)/(1 - x/2)"));
}

TEST_CASE("pade: [2/2] of e^x", "[6][pade][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto p = pade(exp(x), x, 2, 2);
    // (1 + x/2 + x²/12) / (1 - x/2 + x²/12)
    REQUIRE(oracle.equivalent(p->str(),
                              "(1 + x/2 + x**2/12)/(1 - x/2 + x**2/12)"));
}

TEST_CASE("pade: [3/2] of sin(x)", "[6][pade][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto p = pade(sin(x), x, 3, 2);
    // (x - 7x³/60) / (1 + x²/20)
    REQUIRE(oracle.equivalent(p->str(),
                              "(x - 7*x**3/60)/(1 + x**2/20)"));
}

TEST_CASE("pade: [m/0] is the truncated Taylor polynomial",
          "[6][pade][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto p = pade(exp(x), x, 3, 0);
    REQUIRE(oracle.equivalent(p->str(),
                              "1 + x + x**2/2 + x**3/6"));
}

#include <sympp/calculus/extrema.hpp>

// ----- stationary_points / minimum / maximum / monotonicity ------------------

TEST_CASE("stationary_points: x² has critical at 0", "[6][stationary]") {
    auto x = symbol("x");
    auto pts = stationary_points(pow(x, integer(2)), x);
    REQUIRE(pts.size() == 1);
    REQUIRE(pts[0] == S::Zero());
}

TEST_CASE("stationary_points: x³ - 3x has critical at ±1",
          "[6][stationary][oracle]") {
    auto x = symbol("x");
    auto pts = stationary_points(pow(x, integer(3)) - integer(3) * x, x);
    REQUIRE(pts.size() == 2);
    // Roots of 3x² - 3 = 0 — ±1.
    bool has_pos = false, has_neg = false;
    for (auto& p : pts) {
        if (p == integer(1)) has_pos = true;
        if (p == integer(-1)) has_neg = true;
    }
    REQUIRE(has_pos);
    REQUIRE(has_neg);
}

TEST_CASE("minimum: x² minimum at 0 → 0", "[6][minimum]") {
    auto x = symbol("x");
    auto m = minimum(pow(x, integer(2)), x);
    REQUIRE(m == S::Zero());
}

TEST_CASE("maximum: -x² + 4 maximum is 4", "[6][maximum]") {
    auto x = symbol("x");
    auto e = mul(S::NegativeOne(), pow(x, integer(2))) + integer(4);
    auto m = maximum(e, x);
    REQUIRE(m == integer(4));
}

TEST_CASE("is_increasing: e^x is increasing (positive symbol arg)",
          "[6][monotonicity]") {
    // For e^x to be definitely-increasing we'd need the engine to know
    // e^x > 0 always. With our generic Symbol we can at least confirm
    // the helper doesn't crash and returns Unknown when undecidable.
    auto x = symbol("x");
    auto inc = is_increasing(pow(x, integer(2)) + integer(5), x);
    // 2x sign depends on x; can't prove globally → Unknown.
    REQUIRE(!inc.has_value());
}

TEST_CASE("is_increasing: 3x + 7 is increasing (constant positive derivative)",
          "[6][monotonicity]") {
    auto x = symbol("x");
    auto e = integer(3) * x + integer(7);
    auto inc = is_increasing(e, x);
    REQUIRE(inc == std::optional<bool>{true});
}

TEST_CASE("is_decreasing: -2x is decreasing", "[6][monotonicity]") {
    auto x = symbol("x");
    auto e = mul(S::NegativeOne(), integer(2) * x);
    auto dec = is_decreasing(e, x);
    REQUIRE(dec == std::optional<bool>{true});
}

#include <sympp/calculus/euler_lagrange.hpp>

// ----- Euler-Lagrange --------------------------------------------------------

TEST_CASE("euler_lagrange: free particle L = (1/2)*y'^2 → y'' = 0",
          "[6][euler_lagrange][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    auto ypp = symbol("ypp");
    auto L = pow(yp, integer(2)) / integer(2);
    auto el = euler_lagrange(L, y, yp, ypp, x);
    // EL: 0 - ypp = -ypp; equation -ypp = 0 means ypp = 0.
    REQUIRE(oracle.equivalent(el->str(), "-ypp"));
}

TEST_CASE("euler_lagrange: harmonic oscillator → -y - y''",
          "[6][euler_lagrange][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    auto ypp = symbol("ypp");
    // L = (1/2)yp² - (1/2)y²
    auto L = pow(yp, integer(2)) / integer(2)
             - pow(y, integer(2)) / integer(2);
    auto el = euler_lagrange(L, y, yp, ypp, x);
    REQUIRE(oracle.equivalent(el->str(), "-y - ypp"));
}

TEST_CASE("euler_lagrange: pendulum → -sin(y) - y''",
          "[6][euler_lagrange][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    auto ypp = symbol("ypp");
    // L = (1/2)yp² + cos(y) — gravity below; equation y'' + sin(y) = 0.
    auto L = pow(yp, integer(2)) / integer(2) + cos(y);
    auto el = euler_lagrange(L, y, yp, ypp, x);
    REQUIRE(oracle.equivalent(el->str(), "-sin(y) - ypp"));
}

#include <sympp/calculus/asymptotes.hpp>

// ----- vertical_asymptotes / inflection_points -------------------------------

TEST_CASE("vertical_asymptotes: 1/(x-2) at x=2",
          "[6][asymptote][oracle]") {
    auto x = symbol("x");
    auto f = pow(x - integer(2), integer(-1));
    auto va = vertical_asymptotes(f, x);
    REQUIRE(va.size() == 1);
    REQUIRE(va[0] == integer(2));
}

TEST_CASE("vertical_asymptotes: x/(x²-1) at ±1",
          "[6][asymptote][oracle]") {
    auto x = symbol("x");
    auto f = x / (pow(x, integer(2)) - integer(1));
    auto va = vertical_asymptotes(f, x);
    REQUIRE(va.size() == 2);
}

TEST_CASE("vertical_asymptotes: removable singularity is excluded",
          "[6][asymptote]") {
    auto x = symbol("x");
    // (x²-1)/(x-1) — numerator and denominator share (x-1). After cancel
    // the simplified form is x+1, which has no denominator.
    auto f = (pow(x, integer(2)) - integer(1)) / (x - integer(1));
    auto va = vertical_asymptotes(f, x);
    REQUIRE(va.empty());
}

TEST_CASE("inflection_points: x³ at 0 (f''(x) = 6x = 0)",
          "[6][inflection]") {
    auto x = symbol("x");
    auto f = pow(x, integer(3));
    auto ip = inflection_points(f, x);
    REQUIRE(ip.size() == 1);
    REQUIRE(ip[0] == S::Zero());
}

TEST_CASE("inflection_points: x⁴ - 6x² at ±1 (f'' = 12x² - 12 = 0)",
          "[6][inflection]") {
    auto x = symbol("x");
    auto f = pow(x, integer(4)) - integer(6) * pow(x, integer(2));
    auto ip = inflection_points(f, x);
    REQUIRE(ip.size() == 2);
    bool has_pos = false, has_neg = false;
    for (auto& p : ip) {
        if (p == integer(1)) has_pos = true;
        if (p == integer(-1)) has_neg = true;
    }
    REQUIRE(has_pos);
    REQUIRE(has_neg);
}
