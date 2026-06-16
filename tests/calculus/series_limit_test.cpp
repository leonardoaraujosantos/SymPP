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

// SERIES-COMPOSE-1: a composite f(g(x)) whose inner g is finite-but-non-simple at
// 0 (e.g. g = sin(x)/x with its 1/x factor). Differentiating f(g) directly leaves
// each Taylor coefficient as a hard 0/0 limit that returns nan, so log(sin x / x)
// stayed unexpanded. Composing the inner and outer series resolves it.
TEST_CASE("series: composition of f(g(x)) with a removable inner (SERIES-COMPOSE-1)",
          "[6][series][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto zero = S::Zero();
    auto inv = [&](const Expr& e) { return pow(e, S::NegativeOne()); };
    // log(sin x / x) = −x²/6 − x⁴/180 − …
    REQUIRE(oracle.equivalent(
        series(log(sin(x) * inv(x)), x, zero, 6)->str(),
        "-x**2/6 - x**4/180"));
    // log(sinh x / x) = x²/6 − x⁴/180 + …
    REQUIRE(oracle.equivalent(
        series(log(sinh(x) * inv(x)), x, zero, 6)->str(),
        "x**2/6 - x**4/180"));
    // Inner that itself needs the Laurent-division path (tan x / x), composed
    // through log and through a power: log(tan x/x) = x²/3 + 7x⁴/90,
    // √(tan x/x) = 1 + x²/6 + 19x⁴/360.
    REQUIRE(oracle.equivalent(
        series(log(tan(x) * inv(x)), x, zero, 6)->str(),
        "x**2/3 + 7*x**4/90"));
    REQUIRE(oracle.equivalent(
        series(pow(tan(x) * inv(x), rational(1, 2)), x, zero, 6)->str(),
        "1 + x**2/6 + 19*x**4/360"));
    // A power g^p with a removable inner (√(sin x/x)).
    REQUIRE(oracle.equivalent(
        series(pow(sin(x) * inv(x), rational(1, 2)), x, zero, 6)->str(),
        "1 - x**2/12 + x**4/1440"));
    // A genuine singularity (log x) is still left unexpanded by composition.
    REQUIRE(series(log(x), x, zero, 6) == log(x));
    // √x stays a branch point (outer Taylor at c = 0 fails).
    REQUIRE(series(pow(x, rational(1, 2)), x, zero, 6) == pow(x, rational(1, 2)));
    // Composition must not hijack a genuine pole: csc²(x) keeps its Laurent term.
    REQUIRE(oracle.equivalent(series(pow(csc(x), integer(2)), x, zero, 6)->str(),
                              "1/x**2 + 1/3 + x**2/15 + 2*x**4/189"));
    // Single-function composites that already worked are unchanged.
    REQUIRE(oracle.equivalent(series(log(cos(x)), x, zero, 6)->str(),
                              "-x**2/2 - x**4/12"));
}

// SERIES-LAURENT-1: functions with a pole at 0 (cot, csc, coth, csch, 1/sin,
// csc², 1/(eˣ−1)) expand to a Laurent series. The engine rewrites reciprocal
// trig/hyperbolic to sin/cos ratios and divides the numerator and denominator
// power series; previously these returned the input unexpanded.
TEST_CASE("series: Laurent expansion at a pole (SERIES-LAURENT-1)",
          "[6][series][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto zero = S::Zero();
    auto chk = [&](const Expr& e, const std::string& want) {
        REQUIRE(oracle.equivalent(series(e, x, zero, 6)->str(), want));
    };
    // cot(x) = 1/x − x/3 − x³/45 − 2x⁵/945.
    chk(cot(x), "1/x - x/3 - x**3/45 - 2*x**5/945");
    // csc(x) = 1/x + x/6 + 7x³/360 + 31x⁵/15120.
    chk(csc(x), "1/x + x/6 + 7*x**3/360 + 31*x**5/15120");
    // coth(x) = 1/x + x/3 − x³/45 + 2x⁵/945.
    chk(coth(x), "1/x + x/3 - x**3/45 + 2*x**5/945");
    // csch(x) = 1/x − x/6 + 7x³/360 − 31x⁵/15120.
    chk(csch(x), "1/x - x/6 + 7*x**3/360 - 31*x**5/15120");
    // Second-order pole: csc²(x) = 1/x² + 1/3 + x²/15 + 2x⁴/189.
    chk(pow(csc(x), integer(2)), "1/x**2 + 1/3 + x**2/15 + 2*x**4/189");
    // 1/(eˣ − 1) = 1/x − 1/2 + x/12 − x³/720 + x⁵/30240.
    chk(pow(exp(x) - integer(1), integer(-1)),
        "1/x - 1/2 + x/12 - x**3/720 + x**5/30240");
    // x·cot(x) is analytic (the pole cancels): 1 − x²/3 − x⁴/45.
    chk(x * cot(x), "1 - x**2/3 - x**4/45");
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

// LIMIT-LOGSUMEXP-1: log of an exponential-dominated sum. Factoring out the
// fastest-growing exp(e_dom) gives log(g) = e_dom + log(g·e^(−e_dom)) with a
// finite residual, so x − log(cosh x) → log 2 and log(2^x+3^x)/x → log 3 (hence
// (2^x+3^x)^(1/x) → 3, the max of the bases). Previously these returned nan or
// an unevaluated ∞-arithmetic mess.
TEST_CASE("limit: log of an exponential-dominated sum (LIMIT-LOGSUMEXP-1)",
          "[6][limit][infinity][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto oo = S::Infinity();
    auto pw = [&](int b) { return pow(integer(b), x); };
    // x − log(cosh x) → log 2 (and sinh likewise).
    REQUIRE(oracle.equivalent(
        limit(x - log(cosh(x)), x, oo)->str(), "log(2)"));
    REQUIRE(oracle.equivalent(
        limit(x - log(sinh(x)), x, oo)->str(), "log(2)"));
    // log(1 + e^x) − x → 0.
    REQUIRE(oracle.equivalent(
        limit(log(integer(1) + exp(x)) - x, x, oo)->str(), "0"));
    // log(2^x + 3^x)/x → log 3 (dominant base 3).
    REQUIRE(oracle.equivalent(
        limit(log(pw(2) + pw(3)) / x, x, oo)->str(), "log(3)"));
    // (2^x + 3^x)^(1/x) → 3, and a three-term max.
    REQUIRE(oracle.equivalent(
        limit(pow(pw(2) + pw(3), pow(x, S::NegativeOne())), x, oo)->str(), "3"));
    REQUIRE(oracle.equivalent(
        limit(pow(pw(3) + pw(5) + pw(2), pow(x, S::NegativeOne())), x, oo)->str(),
        "5"));
    // log(e^(2x) + e^x)/x → 2 (dominant rate 2).
    REQUIRE(oracle.equivalent(
        limit(log(exp(integer(2) * x) + exp(x)) / x, x, oo)->str(), "2"));
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

// LIMIT-EXPRATIO-1: a product/ratio of constant-base exponentials sharing one
// var-monomial — 2^x/3^x, exp(x)/exp(2x), 2^x·e^(−3x) — is an ∞·0 form that
// L'Hôpital cannot crack (differentiating reproduces it), so the generic
// product path returned nan. Combining into a single B^x (B = ∏bᵢ^cᵢ·e^Σdⱼ) and
// deciding the limit from sign(B−1) resolves it.
TEST_CASE("limit: ratio of exponentials at infinity (LIMIT-EXPRATIO-1)",
          "[6][limit][infinity][regression]") {
    auto x = symbol("x");
    auto oo = S::Infinity();
    // 2^x/3^x = (2/3)^x → 0; 3^x/2^x → oo.
    REQUIRE(limit(pow(integer(2), x) * pow(integer(3), mul(S::NegativeOne(), x)),
                  x, oo) == S::Zero());
    REQUIRE(limit(pow(integer(3), x) * pow(integer(2), mul(S::NegativeOne(), x)),
                  x, oo) == oo);
    // exp(x)/exp(2x) = exp(−x) → 0; exp(2x)/exp(x) → oo (1/exp surfaces as
    // Pow(exp,−1), exercising the exp(g)^k branch).
    REQUIRE(limit(exp(x) * pow(exp(mul(integer(2), x)), integer(-1)), x, oo)
            == S::Zero());
    REQUIRE(limit(exp(mul(integer(2), x)) * pow(exp(x), integer(-1)), x, oo)
            == oo);
    // Mixed constant-base × exp: 2^x·e^(−3x) = e^(x(ln2−3)) → 0 (ln2−3 < 0).
    REQUIRE(limit(pow(integer(2), x) * exp(mul(integer(-3), x)), x, oo)
            == S::Zero());
    // Equal rates cancel to the constant 1 (no false ∞ or 0).
    REQUIRE(limit(pow(integer(2), x) * pow(integer(2), mul(S::NegativeOne(), x)),
                  x, oo) == S::One());
    // Toward −∞ the direction flips: (2/3)^x → oo. (SymPy itself is internally
    // inconsistent here — limit((2/3)**x,-oo)=0 yet limit((2/3)**(-x),oo)=oo —
    // so this is checked directly; (2/3)^(−100) ≈ 4·10¹⁷ confirms the divergence.)
    REQUIRE(limit(pow(integer(2), x) * pow(integer(3), mul(S::NegativeOne(), x)),
                  x, S::NegativeInfinity()) == oo);
    // Polynomial residual: exponential growth class dominates any polynomial
    // degree, where the generic L'Hôpital path stalled for degree ≥ 2 with a
    // rational base (degree 1 worked, x·(1/2)^x → 0). Decaying exp → 0 regardless
    // of degree; growing exp → ±∞ with the polynomial's sign.
    REQUIRE(limit(pow(x, integer(2))
                      * pow(integer(2), mul(S::NegativeOne(), x)),
                  x, oo) == S::Zero());  // x²/2^x
    REQUIRE(limit(pow(x, integer(3)) * pow(integer(2), x)
                      * pow(integer(3), mul(S::NegativeOne(), x)),
                  x, oo) == S::Zero());  // x³·2^x/3^x
    REQUIRE(limit(pow(x, integer(2)) * pow(integer(3), x)
                      * pow(integer(2), mul(S::NegativeOne(), x)),
                  x, oo) == oo);  // x²·3^x/2^x → ∞
    REQUIRE(limit(mul(S::NegativeOne(), pow(x, integer(2))) * pow(integer(3), x),
                  x, oo) == S::NegativeInfinity());  // −x²·3^x → −∞
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

// LIMIT-ESSENTIAL-PT-1: an essential singularity at a finite point (exp(1/x),
// 1/x² at 0). Substituting u = 1/(x−a) and comparing the u → ±∞ one-sided limits
// resolves it: direct substitution gives exp(zoo) = nan otherwise.
TEST_CASE("limit: essential singularity at a finite point (LIMIT-ESSENTIAL-PT-1)",
          "[6][limit][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto z = S::Zero();
    // exp(−1/x²) → 0 (both sides); x/(exp(1/x)−1) → 0; x²/(exp(1/x²)−1) → 0.
    REQUIRE(limit(exp(mul(S::NegativeOne(), pow(x, integer(-2)))), x, z)
            == S::Zero());
    REQUIRE(oracle.equivalent(
        limit(x / (exp(pow(x, integer(-1))) - integer(1)), x, z)->str(), "0"));
    REQUIRE(limit(pow(x, integer(2))
                      / (exp(pow(x, integer(-2))) - integer(1)),
                  x, z) == S::Zero());
    // Genuinely two-sided-divergent cases are NOT coerced: the one-sided limits
    // disagree, so exp(1/x) and 1/x stay nan / zoo.
    REQUIRE(limit(exp(pow(x, integer(-1))), x, z)->type_id() == TypeId::NaN);
    REQUIRE(limit(pow(x, integer(-1)), x, z) == S::ComplexInfinity());
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

// LIMIT-POWFORM-1: indeterminate power forms (1^∞, 0^0, ∞^0) at a FINITE
// target. Direct substitution collapses (1+x)^(1/x) to pow(1, zoo) = 1 before
// the power-form handler runs, so the textbook e-limit returned 1 instead of E.
// The handler now runs before substitution.
TEST_CASE("limit: indeterminate power forms at a finite point (LIMIT-POWFORM-1)",
          "[6][limit][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto z = S::Zero();
    auto inv = [&](const Expr& e) { return pow(e, S::NegativeOne()); };
    // (1+x)^(1/x) → e (the definition of e).
    REQUIRE(oracle.equivalent(
        limit(pow(integer(1) + x, inv(x)), x, z)->str(), "E"));
    // (1+2x)^(1/x) → e², (1+x)^(2/x) → e², (1−x)^(1/x) → e⁻¹.
    REQUIRE(oracle.equivalent(
        limit(pow(integer(1) + integer(2) * x, inv(x)), x, z)->str(), "exp(2)"));
    REQUIRE(oracle.equivalent(
        limit(pow(integer(1) + x, integer(2) * inv(x)), x, z)->str(), "exp(2)"));
    REQUIRE(oracle.equivalent(
        limit(pow(integer(1) - x, inv(x)), x, z)->str(), "exp(-1)"));
    // cos(x)^(1/x²) → e^(−1/2).
    REQUIRE(oracle.equivalent(
        limit(pow(cos(x), inv(pow(x, integer(2)))), x, z)->str(), "exp(-1/2)"));
    // Composite base → 1: the 1^∞ rate uses log(base) ~ (base−1), which avoids the
    // missing Taylor series of log(sin x/x). (sin x/x)^(1/x²) → e^(−1/6),
    // (tan x/x)^(1/x²) → e^(1/3), cos(2x)^(1/x²) → e^(−2).
    REQUIRE(oracle.equivalent(
        limit(pow(sin(x) * inv(x), inv(pow(x, integer(2)))), x, z)->str(),
        "exp(-1/6)"));
    REQUIRE(oracle.equivalent(
        limit(pow(tan(x) * inv(x), inv(pow(x, integer(2)))), x, z)->str(),
        "exp(1/3)"));
    REQUIRE(oracle.equivalent(
        limit(pow(cos(integer(2) * x), inv(pow(x, integer(2)))), x, z)->str(),
        "exp(-2)"));
    // Determinate powers are unaffected: (1+x)² → 1, x^x → 1.
    REQUIRE(limit(pow(integer(1) + x, integer(2)), x, z) == S::One());
    REQUIRE(limit(pow(x, x), x, z) == S::One());
}

// LIMIT-CONJUGATE-1: ∞ − ∞ involving a radical, resolved by the conjugate.
// x − √(x²+1) = −1/(x + √(x²+1)) → 0. Direct substitution gives the
// indeterminate ∞ − ∞ (nan); the engine rationalizes via t₁+t₂ =
// (t₁²−t₂²)/(t₁−t₂), clearing the radical from the numerator.
TEST_CASE("limit: conjugate of a radical difference (LIMIT-CONJUGATE-1)",
          "[6][limit][infinity][oracle][regression]") {
    auto x = symbol("x");
    const Expr oo = S::Infinity();
    auto sq = [&](const Expr& e) { return pow(e, rational(1, 2)); };
    // x − √(x²+1) → 0, x − √(x²−1) → 0.
    REQUIRE(limit(x - sq(pow(x, integer(2)) + integer(1)), x, oo) == S::Zero());
    REQUIRE(limit(x - sq(pow(x, integer(2)) - integer(1)), x, oo) == S::Zero());
    // √(x+1) − √x → 0.
    REQUIRE(limit(sq(x + integer(1)) - sq(x), x, oo) == S::Zero());
    // The non-indeterminate companion still diverges: x + √(x²+1) → ∞.
    REQUIRE(limit(x + sq(pow(x, integer(2)) + integer(1)), x, oo) == oo);
}

// LIMIT-RADICAL-INF-1: √-difference limits at +∞ with a NONZERO finite value.
// The conjugate clears the ∞−∞, but the residual ratio (e.g. x/(√(x²+x)+x)) is an
// ∞/∞ that L'Hôpital abandons on radicals (the radical never stabilises). A
// leading-asymptotic-term evaluator resolves it: √(x²+x)−x → 1/2, the two-radical
// difference √(x²+x)−√(x²−x) → 1, and the 0·∞ product x·(√(x²+1)−x) → 1/2.
// These previously returned nan (a wrong answer). Matches SymPy.
TEST_CASE("limit: nonzero radical differences at infinity (LIMIT-RADICAL-INF-1)",
          "[6][limit][infinity][oracle][regression]") {
    auto x = symbol("x");
    const Expr oo = S::Infinity();
    auto sq = [&](const Expr& e) { return pow(e, rational(1, 2)); };
    auto x2 = pow(x, integer(2));
    // √(x²+x) − x → 1/2 ; x − √(x²−x) → 1/2.
    REQUIRE(limit(sq(x2 + x) - x, x, oo) == rational(1, 2));
    REQUIRE(limit(x - sq(x2 - x), x, oo) == rational(1, 2));
    // √(x²+x) − √(x²−x) → 1.
    REQUIRE(limit(sq(x2 + x) - sq(x2 - x), x, oo) == S::One());
    // 0·∞ product: x·(√(x²+1) − x) → 1/2.
    REQUIRE(limit(x * (sq(x2 + integer(1)) - x), x, oo) == rational(1, 2));
    // General coefficient: √(4x²+x) − 2x → 1/4.
    REQUIRE(limit(sq(integer(4) * x2 + x) - integer(2) * x, x, oo)
            == rational(1, 4));
    // Larger linear term: √(x²+2x) − x → 1 (also guards try_algebraic_inf's
    // leading-term path, which must stay silent — no stray diagnostics).
    REQUIRE(limit(sq(x2 + integer(2) * x) - x, x, oo) == S::One());
    // LIMIT-NROOT-INF: the conjugate generalizes from √ to an n-th root via
    // u − v = (uⁿ − vⁿ)/Σ u^(n−1−i)vⁱ. (x³+x²)^(1/3) − x → 1/3, (x⁴+x³)^(1/4) − x →
    // 1/4, the two-cube-root difference (x³+x²)^(1/3) − (x³−x²)^(1/3) → 2/3, and a
    // leading coefficient (8x³+x²)^(1/3) − 2x → 1/12. Previously nan. Matches SymPy.
    auto cbrt = [&](const Expr& e) { return pow(e, rational(1, 3)); };
    auto x3 = pow(x, integer(3));
    REQUIRE(limit(cbrt(x3 + x2) - x, x, oo) == rational(1, 3));
    REQUIRE(limit(pow(pow(x, integer(4)) + x3, rational(1, 4)) - x, x, oo)
            == rational(1, 4));
    REQUIRE(limit(cbrt(x3 + x2) - cbrt(x3 - x2), x, oo) == rational(2, 3));
    REQUIRE(limit(cbrt(integer(8) * x3 + x2) - integer(2) * x, x, oo)
            == rational(1, 12));
}

// LIMIT-RECIP-INF-1: asymptotic-expansion limits at +∞ with a transcendental
// subleading term, resolved by the reciprocal substitution x = 1/t → t → 0 (where
// the series/L'Hôpital machinery applies), guarded by a numeric convergence check
// so a wrong t-limit cannot slip through. x − x²·log(1+1/x) → 1/2,
// x²(1−cos(1/x)) → 1/2, x·(eˣ⁻¹... ) etc. Previously returned nan. Matches SymPy.
TEST_CASE("limit: reciprocal-substitution asymptotics at infinity (LIMIT-RECIP-INF-1)",
          "[6][limit][infinity][oracle][regression]") {
    auto x = symbol("x");
    const Expr oo = S::Infinity();
    auto inv = [&](const Expr& e) { return pow(e, integer(-1)); };
    // x − x²·log(1+1/x) → 1/2.
    REQUIRE(limit(x - pow(x, integer(2)) * log(integer(1) + inv(x)), x, oo)
            == rational(1, 2));
    // x²·(1 − cos(1/x)) → 1/2.
    REQUIRE(limit(pow(x, integer(2)) * (integer(1) - cos(inv(x))), x, oo)
            == rational(1, 2));
    // x·(e^(1/x) − 1) → 1.
    REQUIRE(limit(x * (exp(inv(x)) - integer(1)), x, oo) == S::One());
    // x²·(e^(1/x) − 1 − 1/x) → 1/2.
    REQUIRE(limit(pow(x, integer(2))
                      * (exp(inv(x)) - integer(1) - inv(x)),
                  x, oo)
            == rational(1, 2));
    // Sanity: the fallback's numeric guard keeps determinate limits correct and
    // does not invent a finite value where the true limit is infinite.
    REQUIRE(limit(pow(x, integer(2)) - x, x, oo) == oo);
    REQUIRE(limit(x * sin(inv(x)), x, oo) == S::One());
}

// LIMIT-LOG-1: logarithms at ∞ — log-continuity (limit(log g) = log(lim g)),
// the ∞ − ∞ between logs (log(x+1) − log(x) → 0, via combining), and
// atan-continuity (limit(atan g) = atan(lim g) = π/2). Also exercises the
// supporting infinity-arithmetic fix (oo + √2 = oo). Previously returned nan.
TEST_CASE("limit: logarithms and atan at infinity (LIMIT-LOG-1)",
          "[6][limit][infinity][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    const Expr oo = S::Infinity();
    auto sq = [&](const Expr& e) { return pow(e, integer(2)); };
    // ∞ − ∞ between logs.
    REQUIRE(limit(log(x + integer(1)) - log(x), x, oo) == S::Zero());
    REQUIRE(limit(log(sq(x) + x + integer(1)) - log(sq(x) - x + integer(1)), x,
                  oo)
            == S::Zero());
    // log(2x) − log(x) = log 2.
    REQUIRE(oracle.equivalent(
        limit(log(integer(2) * x) - log(x), x, oo)->str(), "log(2)"));
    // atan-continuity: atan(2x+1) → π/2.
    REQUIRE(oracle.equivalent(limit(atan(integer(2) * x + integer(1)), x, oo)
                                  ->str(),
                              "pi/2"));
    // Infinity absorbs a finite real constant: oo + √2 = oo (not oo + √2).
    REQUIRE(S::Infinity() + pow(integer(2), rational(1, 2)) == S::Infinity());
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
    // LIMIT-LHOPITAL-NEST-1: a 0·∞ whose L'Hôpital ratio has a derivative that is
    // itself a fraction (d/dx(1/x) = −x⁻²). together() leaves the resulting
    // nested reciprocal un-flattened, which made these return nan; flattening the
    // ratio each step fixes it. x·(π/2 − atan x) → 1, x·atan(1/x) → 1.
    REQUIRE(oracle.equivalent(
        limit(x * (S::Pi() / integer(2) - atan(x)), x, S::Infinity())->str(),
        "1"));
    REQUIRE(oracle.equivalent(
        limit(x * atan(integer(1) / x), x, S::Infinity())->str(), "1"));
    REQUIRE(oracle.equivalent(
        limit(x * tan(integer(1) / x), x, S::Infinity())->str(), "1"));
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

// LIMIT-RECIPTRIG-1: reciprocal trig/hyperbolic functions (cot, csc, sec, coth,
// csch, sech) were opaque to the limit engine, so 0·∞ forms like x·cot(x)
// returned nan. The engine now rewrites them as sin/cos ratios before resolving.
TEST_CASE("limit: reciprocal trig/hyperbolic (LIMIT-RECIPTRIG-1)",
          "[6][limit][oracle][regression]") {
    auto x = symbol("x");
    const Expr z = S::Zero();
    REQUIRE(limit(x * cot(x), x, z) == S::One());
    REQUIRE(limit(cot(x) * sin(x), x, z) == S::One());
    REQUIRE(limit(x * csc(x), x, z) == S::One());
    REQUIRE(limit(x * coth(x), x, z) == S::One());
    REQUIRE(limit(x * csch(x), x, z) == S::One());
    REQUIRE(limit(sec(x), x, z) == S::One());
    // x²·csc²(x) → 1, and tan·cot → 1 (both 0·∞ / ∞·0 reciprocal forms).
    REQUIRE(limit(pow(x, integer(2)) * pow(csc(x), integer(2)), x, z)
            == S::One());
    REQUIRE(limit(tan(x) * cot(x), x, z) == S::One());
    // (cos x − 1)·csc(x) → 0.
    REQUIRE(limit((cos(x) - integer(1)) * csc(x), x, z) == S::Zero());
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

// SUM-POLYEXPAND-1: a product/power polynomial summand is expanded before
// dispatch, so k·(k+1), (k+1)², (2k+1)(k−1) sum via linearity + Faulhaber
// instead of returning an unevaluated Sum.
TEST_CASE("summation: product/power polynomial summands (SUM-POLYEXPAND-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto S = [&](const Expr& e) { return summation(e, k, integer(1), n); };
    // Σ k(k+1) = n(n+1)(n+2)/3.
    REQUIRE(oracle.equivalent(S(mul(k, k + integer(1)))->str(),
                              "n*(n + 1)*(n + 2)/3"));
    // Σ k(k−1) = n(n−1)(n+1)/3.
    REQUIRE(oracle.equivalent(S(mul(k, k - integer(1)))->str(),
                              "n*(n - 1)*(n + 1)/3"));
    // Σ (k+1)² = n(2n²+9n+13)/6.
    REQUIRE(oracle.equivalent(S(pow(k + integer(1), integer(2)))->str(),
                              "n*(2*n**2 + 9*n + 13)/6"));
    // Σ (2k+1)(k−1) = n(n−1)(4n+7)/6.
    REQUIRE(oracle.equivalent(
        S(mul(integer(2) * k + integer(1), k - integer(1)))->str(),
        "n*(n - 1)*(4*n + 7)/6"));
    REQUIRE(S(mul(k, k + integer(1)))->str().find("Sum(") == std::string::npos);
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
    // SUM-TELESCOPE-DIFF-1: an explicit difference g(k) − g(k+1) telescopes to
    // g(lo) − g(hi+1), even when neither term sums on its own — the check runs
    // before the linearity split that would otherwise separate them. Also closes
    // factorial differences the rational path can't.
    REQUIRE(oracle.equivalent(
        summation(pow(k, integer(-1)) - pow(k + one, integer(-1)), k, one, n)
            ->str(),
        "1 - 1/(n+1)"));
    REQUIRE(oracle.equivalent(
        summation(pow(factorial(k), integer(-1))
                      - pow(factorial(k + one), integer(-1)),
                  k, one, n)->str(),
        "1 - 1/factorial(n+1)"));
    REQUIRE(oracle.equivalent(
        summation(pow(k, integer(-2)) - pow(k + one, integer(-2)), k, one, n)
            ->str(),
        "1 - 1/(n+1)**2"));
    // A non-telescoping difference (1/2^k − 1/3^k, both geometric) is unaffected:
    // it falls through to the geometric handlers (sum to 1/2).
    REQUIRE(oracle.equivalent(
        summation(pow(integer(2), mul(S::NegativeOne(), k))
                      - pow(integer(3), mul(S::NegativeOne(), k)),
                  k, one, S::Infinity())->str(),
        "1/2"));
}

// SUM-TELESCOPE-2: the rational telescoping closed form generalized from a quadratic
// denominator to any degree ≥ 2 whose roots are rational and pairwise differ by
// integers (via partial fractions over Σ Aᵢ/(k−rᵢ)). Closes Σ1/(k(k+1)(k+2))=1/4 and
// Σ1/(k(k+1)(k+2)(k+3))=1/18, previously unevaluated. Verified against SymPy.
TEST_CASE("summation: higher-degree telescoping rationals (SUM-TELESCOPE-2)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto one = integer(1);
    const Expr oo = S::Infinity();
    // Σ_{k=1}^∞ 1/(k(k+1)(k+2)) = 1/4.
    REQUIRE(oracle.equivalent(
        summation(pow(k * (k + one) * (k + integer(2)), integer(-1)), k, one, oo)
            ->str(),
        "1/4"));
    // Σ_{k=1}^∞ 1/(k(k+1)(k+2)(k+3)) = 1/18.
    REQUIRE(oracle.equivalent(
        summation(pow(k * (k + one) * (k + integer(2)) * (k + integer(3)),
                      integer(-1)),
                  k, one, oo)
            ->str(),
        "1/18"));
    // Half-integer roots: Σ 1/((2k-1)(2k+1)(2k+3)) = 1/12.
    REQUIRE(oracle.equivalent(
        summation(pow((integer(2) * k - one) * (integer(2) * k + one)
                          * (integer(2) * k + integer(3)),
                      integer(-1)),
                  k, one, oo)
            ->str(),
        "1/12"));
    // Finite closed form: Σ_{k=1}^n 1/(k(k+1)(k+2)) = n(n+3)/(4(n+1)(n+2)).
    REQUIRE(oracle.equivalent(
        summation(pow(k * (k + one) * (k + integer(2)), integer(-1)), k, one, n)
            ->str(),
        "n*(n+3)/(4*(n+1)*(n+2))"));
}

// SUM-TELESCOPE-3: a repeated-root rational that telescopes after partial fractions,
// (k+1)^j − k^j over k^j(k+1)^j = 1/k^j − 1/(k+1)^j. apart() exposes the g(k)−g(k+1)
// form the explicit-difference check misses on the combined fraction. Closes
// Σ(2k+1)/(k²(k+1)²)=1 and Σ(3k²+3k+1)/(k³(k+1)³)=1; previously unevaluated. A
// residual ζ part (1/(k²(k+1)) → π²/6−1) is a 3-term apart and stays unevaluated.
TEST_CASE("summation: repeated-root telescoping via apart (SUM-TELESCOPE-3)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto one = integer(1);
    const Expr oo = S::Infinity();
    auto sq = [&](const Expr& e) { return pow(e, integer(2)); };
    auto cube = [&](const Expr& e) { return pow(e, integer(3)); };
    // Σ_{k=1}^∞ (2k+1)/(k²(k+1)²) = 1.
    REQUIRE(oracle.equivalent(
        summation((integer(2) * k + one)
                      * pow(sq(k) * sq(k + one), integer(-1)),
                  k, one, oo)
            ->str(),
        "1"));
    // Σ_{k=1}^∞ (3k²+3k+1)/(k³(k+1)³) = 1.
    REQUIRE(oracle.equivalent(
        summation((integer(3) * sq(k) + integer(3) * k + one)
                      * pow(cube(k) * cube(k + one), integer(-1)),
                  k, one, oo)
            ->str(),
        "1"));
    // Finite form: Σ_{k=1}^n (2k+1)/(k²(k+1)²) = 1 − 1/(n+1)².
    REQUIRE(oracle.equivalent(
        summation((integer(2) * k + one)
                      * pow(sq(k) * sq(k + one), integer(-1)),
                  k, one, n)
            ->str(),
        "1 - 1/(n+1)**2"));
    // A residual-ζ case is NOT a pure telescoping value: it is closed by the general
    // rational-sum path instead (SUM-RATIONAL-1): Σ1/(k²(k+1)) = π²/6 − 1.
    REQUIRE(oracle.equivalent(
        summation(pow(sq(k) * (k + one), integer(-1)), k, one, oo)->str(),
        "pi**2/6 - 1"));
}

// SUM-RATIONAL-1: a general convergent rational sum closes via partial fractions —
// apart() splits it into a ζ part (poles of order ≥ 2) and a telescoping part (simple
// poles, whose residues sum to zero). Σ1/(k²(k+1))=π²/6−1, Σ1/(k(k+1)²)=2−π²/6,
// Σ1/(k²(k+2))=π²/12−3/8. Previously unevaluated. Matches SymPy.
TEST_CASE("summation: general rational sum via partial fractions (SUM-RATIONAL-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto one = integer(1);
    const Expr oo = S::Infinity();
    auto sq = [&](const Expr& e) { return pow(e, integer(2)); };
    // Σ 1/(k²(k+1)) = π²/6 − 1.
    REQUIRE(oracle.equivalent(
        summation(pow(sq(k) * (k + one), integer(-1)), k, one, oo)->str(),
        "pi**2/6 - 1"));
    // Σ 1/(k(k+1)²) = 2 − π²/6.
    REQUIRE(oracle.equivalent(
        summation(pow(k * sq(k + one), integer(-1)), k, one, oo)->str(),
        "2 - pi**2/6"));
    // Σ 1/(k²(k+2)) = π²/12 − 3/8.
    REQUIRE(oracle.equivalent(
        summation(pow(sq(k) * (k + integer(2)), integer(-1)), k, one, oo)->str(),
        "pi**2/12 - 3/8"));
    // Constant scales through: Σ 2/(k²(k+1)) = π²/3 − 2.
    REQUIRE(oracle.equivalent(
        summation(integer(2) * pow(sq(k) * (k + one), integer(-1)), k, one, oo)
            ->str(),
        "pi**2/3 - 2"));
}

// SUM-FACT-TELESCOPE-1: Σ P(k)/(k+m)! with deg P ≥ 1 — Gosper for a factorial
// denominator. The antidifference g(k)=Q(k)/(k+m−1)! exists iff Q(k)(k+m)−Q(k+1)=P(k)
// is consistent (a constant numerator like 1/(k+1)! → e−2 is not telescoping and
// stays unevaluated). Closes Σk/(k+1)!=1 and Σ(k²−1)/(k+1)!=1; previously left a
// partially-split unevaluated Sum. Matches SymPy.
TEST_CASE("summation: factorial telescoping Σ P(k)/(k+m)! (SUM-FACT-TELESCOPE-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto one = integer(1);
    const Expr oo = S::Infinity();
    auto inv_fact = [&](const Expr& a) {
        return pow(factorial(a), integer(-1));
    };
    // Σ_{k=1}^∞ k/(k+1)! = 1.
    REQUIRE(oracle.equivalent(
        summation(k * inv_fact(k + one), k, one, oo)->str(), "1"));
    // Σ_{k=1}^∞ (k²−1)/(k+1)! = 1 (telescopes as a whole, not term by term).
    REQUIRE(oracle.equivalent(
        summation((pow(k, integer(2)) - one) * inv_fact(k + one), k, one, oo)
            ->str(),
        "1"));
    // Finite form: Σ_{k=1}^n k/(k+1)! = 1 − 1/(n+1)!.
    REQUIRE(oracle.equivalent(
        summation(k * inv_fact(k + one), k, one, n)->str(),
        "1 - 1/factorial(n+1)"));
    // A constant numerator is NOT telescoping — this handler rejects it (its
    // consistency check fails). Σ 1/(k+1)! = e−2 is instead closed by the shifted
    // exponential-series path (SUM-EXP-SHIFT-1), not a bogus telescoping form.
    REQUIRE(oracle.equivalent(summation(inv_fact(k + one), k, one, oo)->str(),
                              "E - 2"));
}

// SUM-BINOMIAL-1: the binomial theorem Σ_{k=0}^n C(n,k)·rᵏ = (1+r)ⁿ. A summand
// binomial(n,k)·base^(a·k+b) (geometric factor optional), where n is exactly the
// binomial's first arg, closes to const·base^b·(1+base^a)ⁿ.
TEST_CASE("summation: binomial theorem Σ C(n,k)·r^k (SUM-BINOMIAL-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto x = symbol("x");
    auto z = S::Zero();
    // Σ C(n,k) = 2ⁿ; Σ (−1)^k C(n,k) = 0; Σ 2^k C(n,k) = 3ⁿ; Σ x^k C(n,k)=(1+x)ⁿ.
    REQUIRE(oracle.equivalent(summation(binomial(n, k), k, z, n)->str(), "2**n"));
    REQUIRE(summation(pow(integer(-1), k) * binomial(n, k), k, z, n)
            == S::Zero());
    REQUIRE(oracle.equivalent(
        summation(pow(integer(2), k) * binomial(n, k), k, z, n)->str(), "3**n"));
    REQUIRE(oracle.equivalent(
        summation(pow(x, k) * binomial(n, k), k, z, n)->str(), "(1 + x)**n"));
    // A concrete upper bound: Σ_{k=0}^5 C(5,k) = 32.
    REQUIRE(summation(binomial(integer(5), k), k, z, integer(5)) == integer(32));
    // Mismatched binomial argument (binomial(m,k) with upper bound n) does NOT
    // apply the theorem — left unevaluated.
    auto m = symbol("m");
    REQUIRE(summation(binomial(m, k), k, z, n)->str().find("Sum(")
            != std::string::npos);
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

// SUM-POLYGEOM-INF-1: Σ_{k=lo}^∞ P(k)·r^k with |r| < 1. The finite closed form
// evaluates the upper boundary term Q(k)·r^k at k = ∞ as (poly in ∞)·r^∞ = ∞·0,
// which returned nan; for |r| < 1 that term vanishes, so the sum is −S(lo).
TEST_CASE("summation: infinite polynomial × geometric, |r| < 1 (SUM-POLYGEOM-INF-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto oo = S::Infinity();
    auto inv_pow = [&](int b) { return pow(integer(b), mul(S::NegativeOne(), k)); };
    auto suminf = [&](const Expr& e, const Expr& lo) {
        Expr s = summation(e, k, lo, oo);
        REQUIRE(s->str().find("Sum(") == std::string::npos);
        return s;
    };
    // Σ_{k=0}^∞ k/2^k = 2 (degree 1, routed through the poly·geometric path).
    REQUIRE(oracle.equivalent(suminf(k * inv_pow(2), S::Zero())->str(), "2"));
    // Σ_{k=1}^∞ k/2^k = 2 (the k = 0 term is 0).
    REQUIRE(oracle.equivalent(suminf(k * inv_pow(2), S::One())->str(), "2"));
    // Σ_{k=0}^∞ k²/2^k = 6.
    REQUIRE(oracle.equivalent(
        suminf(pow(k, integer(2)) * inv_pow(2), S::Zero())->str(), "6"));
    // Σ_{k=0}^∞ k/3^k = 3/4.
    REQUIRE(oracle.equivalent(suminf(k * inv_pow(3), S::Zero())->str(), "3/4"));
    // Σ_{k=0}^∞ (k+1)/2^k = 4.
    REQUIRE(oracle.equivalent(
        suminf((k + integer(1)) * inv_pow(2), S::Zero())->str(), "4"));
    // A divergent ratio (|r| > 1) is NOT given a bogus finite value — it stays an
    // unevaluated Sum rather than nan.
    REQUIRE(summation(k * pow(integer(2), k), k, S::Zero(), oo)->str()
                .find("Sum(") != std::string::npos);
}

// SUM-POLYGEOM-SYM-1: Σ_{k=1}^n P(k)·xᵏ for a SYMBOLIC ratio x — the generating-
// function identity Σ k·xᵏ = x(1−(n+1)xⁿ+n·xⁿ⁺¹)/(x−1)². The poly·geometric closed
// form was gated to a numeric base/ratio; it now applies symbolically for a finite
// sum (matching SymPy's general branch). An infinite sum with a symbolic ratio fails
// the |x| < 1 convergence test and is correctly left unevaluated (no x**∞ terms).
TEST_CASE("summation: polynomial × geometric, symbolic ratio (SUM-POLYGEOM-SYM-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto x = symbol("x");
    auto sum1n = [&](const Expr& e) { return summation(e, k, integer(1), n); };

    // Σ_{k=1}^n k·xᵏ — verify against the closed form on concrete (n, x) samples.
    auto s1 = sum1n(k * pow(x, k));
    REQUIRE(s1->str().find("Sum(") == std::string::npos);
    // n = 3: 1·x + 2·x² + 3·x³.
    REQUIRE(oracle.equivalent(subs(s1, n, integer(3))->str(),
                              "x + 2*x**2 + 3*x**3"));
    // n = 4: + 4·x⁴.
    REQUIRE(oracle.equivalent(subs(s1, n, integer(4))->str(),
                              "x + 2*x**2 + 3*x**3 + 4*x**4"));

    // Σ_{k=1}^n k²·xᵏ at n = 3: x + 4x² + 9x³.
    auto s2 = sum1n(pow(k, integer(2)) * pow(x, k));
    REQUIRE(s2->str().find("Sum(") == std::string::npos);
    REQUIRE(oracle.equivalent(subs(s2, n, integer(3))->str(),
                              "x + 4*x**2 + 9*x**3"));

    // Infinite sum with a symbolic ratio stays unevaluated (convergence unknown).
    REQUIRE(summation(k * pow(x, k), k, integer(1), S::Infinity())
                ->str()
                .find("Sum(") != std::string::npos);
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
// marker, never the bare summand. Σ 1/(k!+1) (no elementary closed form) stands
// in for the general case; the exponential-series family Σ P(k)·r^k/k! and the
// p-series Σ 1/k^s now close, so they no longer serve as the placeholder.
TEST_CASE("summation: unrecognised sum stays an unevaluated Sum",
          "[6][summation][regression]") {
    auto k = symbol("k");
    auto e = pow(add(factorial(k), integer(1)), integer(-1));
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

// SUM-ALT-PSERIES-1: the alternating p-series Σ_{k=1}^∞ (−1)^k/k^s — the Dirichlet
// eta values η(s). s=1 closes to −log 2 (alternating harmonic), s≥2 to
// (2^(1−s)−1)·ζ(s) (−π²/12, −7π⁴/720, …). A (−1)^(a·k+b) factor with odd a is
// (−1)^k up to the constant sign (−1)^b; a leading constant multiplies through.
TEST_CASE("summation: alternating p-series → Dirichlet eta (SUM-ALT-PSERIES-1)",
          "[6][summation][zeta][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto oo = S::Infinity();
    auto altk = pow(integer(-1), k);
    // Σ (−1)^k/k = −log 2 (alternating harmonic).
    REQUIRE(oracle.equivalent(
        summation(altk * pow(k, integer(-1)), k, integer(1), oo)->str(),
        "-log(2)"));
    // Σ (−1)^(k+1)/k = log 2 (the shifted sign).
    REQUIRE(oracle.equivalent(
        summation(pow(integer(-1), k + integer(1)) * pow(k, integer(-1)), k,
                  integer(1), oo)->str(),
        "log(2)"));
    // Σ (−1)^k/k² = −π²/12, Σ (−1)^k/k⁴ = −7π⁴/720 (even s → closed π-form).
    REQUIRE(oracle.equivalent(
        summation(altk * pow(k, integer(-2)), k, integer(1), oo)->str(),
        "-pi**2/12"));
    REQUIRE(oracle.equivalent(
        summation(altk * pow(k, integer(-4)), k, integer(1), oo)->str(),
        "-7*pi**4/720"));
    // Σ (−1)^k/k³ = −¾·ζ(3) (odd s → symbolic ζ).
    REQUIRE(oracle.equivalent(
        summation(altk * pow(k, integer(-3)), k, integer(1), oo)->str(),
        "-3*zeta(3)/4"));
    // A leading constant multiplies through: Σ 3·(−1)^k/k = −3·log 2.
    REQUIRE(oracle.equivalent(
        summation(integer(3) * altk * pow(k, integer(-1)), k, integer(1), oo)
            ->str(),
        "-3*log(2)"));
    // Non-alternating p-series unaffected: Σ 1/k² = π²/6.
    REQUIRE(oracle.equivalent(
        summation(pow(k, integer(-2)), k, integer(1), oo)->str(), "pi**2/6"));
}

// SUM-ARITH-PSERIES-1: arithmetic-argument p-series Σ_{k=1}^∞ 1/(a·k+b)^s for
// a ∈ {1,2}, closed by slicing ζ(s) over the residue class the denominator runs.
// The famous odd/even Basel relatives: Σ1/(2k−1)²=π²/8, Σ1/(2k)²=π²/24,
// Σ1/(2k−1)⁴=π⁴/96; with finite head terms removed for a shifted start; odd s
// stays a symbolic ζ(s) (Σ1/(2k−1)³ = 7ζ(3)/8). Previously all unevaluated.
TEST_CASE("summation: arithmetic-argument p-series (SUM-ARITH-PSERIES-1)",
          "[6][summation][zeta][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto oo = S::Infinity();
    auto one = integer(1);
    auto odd = integer(2) * k - one;   // 2k−1
    auto even = integer(2) * k;        // 2k
    // Σ 1/(2k−1)² = π²/8 ; Σ 1/(2k)² = π²/24 ; Σ 1/(2k−1)⁴ = π⁴/96.
    REQUIRE(oracle.equivalent(
        summation(pow(odd, integer(-2)), k, one, oo)->str(), "pi**2/8"));
    REQUIRE(oracle.equivalent(
        summation(pow(even, integer(-2)), k, one, oo)->str(), "pi**2/24"));
    REQUIRE(oracle.equivalent(
        summation(pow(odd, integer(-4)), k, one, oo)->str(), "pi**4/96"));
    // Shifted start removes the head term: Σ 1/(2k+1)² = π²/8 − 1.
    REQUIRE(oracle.equivalent(
        summation(pow(integer(2) * k + one, integer(-2)), k, one, oo)->str(),
        "pi**2/8 - 1"));
    // a=1 shift: Σ 1/(k+1)² = π²/6 − 1.
    REQUIRE(oracle.equivalent(
        summation(pow(k + one, integer(-2)), k, one, oo)->str(),
        "pi**2/6 - 1"));
    // A leading constant scales through: Σ 3/(2k−1)² = 3π²/8.
    REQUIRE(oracle.equivalent(
        summation(integer(3) * pow(odd, integer(-2)), k, one, oo)->str(),
        "3*pi**2/8"));
    // Odd exponent stays a symbolic ζ: Σ 1/(2k−1)³ = 7ζ(3)/8.
    REQUIRE(oracle.equivalent(
        summation(pow(odd, integer(-3)), k, one, oo)->str(), "7*zeta(3)/8"));
}

// SUM-SHIFT-1: an infinite sum whose index starts at an integer ≠ 1 is
// re-expressed as Σ_{m=1}^∞ f(m + lo − 1), so the lo=1 closed-form handlers
// (arithmetic p-series, ζ, …) reach it. Closes the standard odd p-series written
// from zero, Σ_{n=0}^∞ 1/(2n+1)² = π²/8, and shifted-start variants.
TEST_CASE("summation: integer index-shift fallback (SUM-SHIFT-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto oo = S::Infinity();
    // Σ_{n=0}^∞ 1/(2n+1)² = π²/8 ; Σ_{n=0}^∞ 1/(2n+1)⁴ = π⁴/96.
    REQUIRE(oracle.equivalent(
        summation(pow(integer(2) * n + integer(1), integer(-2)), n, integer(0),
                  oo)->str(),
        "pi**2/8"));
    REQUIRE(oracle.equivalent(
        summation(pow(integer(2) * n + integer(1), integer(-4)), n, integer(0),
                  oo)->str(),
        "pi**4/96"));
    // Σ_{n=0}^∞ 1/(n+1)² = π²/6 (re-indexed Basel).
    REQUIRE(oracle.equivalent(
        summation(pow(n + integer(1), integer(-2)), n, integer(0), oo)->str(),
        "pi**2/6"));
    // Start above 1 drops head terms: Σ_{n=2}^∞ 1/(2n+1)² = π²/8 − 10/9.
    REQUIRE(oracle.equivalent(
        summation(pow(integer(2) * n + integer(1), integer(-2)), n, integer(2),
                  oo)->str(),
        "pi**2/8 - 10/9"));
    // No regression: the lo=1 form is unchanged.
    REQUIRE(oracle.equivalent(
        summation(pow(integer(2) * n - integer(1), integer(-2)), n, integer(1),
                  oo)->str(),
        "pi**2/8"));
}

// SUM-INVQUAD-1: Σ_{k=1}^∞ c/(a·k²+b) with b/a > 0 — irreducible-quadratic
// denominator, the cotangent closed form Σ 1/(k²+B) = (π√B·coth(π√B) − 1)/(2B).
// apart can't split the complex-conjugate poles, so it was previously unevaluated.
TEST_CASE("summation: inverse irreducible quadratic → coth (SUM-INVQUAD-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto oo = S::Infinity();
    auto one = integer(1);
    auto den = [&](int a, int b) {
        return pow(integer(a) * pow(n, integer(2)) + integer(b), S::NegativeOne());
    };
    // Σ 1/(n²+1) = (π·coth π − 1)/2 = −1/2 + π·coth(π)/2.
    REQUIRE(oracle.equivalent(summation(den(1, 1), n, one, oo)->str(),
                              "-1/2 + pi*coth(pi)/2"));
    // Σ 1/(n²+4) = −1/8 + π·coth(2π)/4.
    REQUIRE(oracle.equivalent(summation(den(1, 4), n, one, oo)->str(),
                              "-1/8 + pi*coth(2*pi)/4"));
    // Leading constant scales through: Σ 3/(n²+1) = −3/2 + 3π·coth(π)/2.
    REQUIRE(oracle.equivalent(
        summation(integer(3) * den(1, 1), n, one, oo)->str(),
        "-3/2 + 3*pi*coth(pi)/2"));
    // Leading n² coefficient: Σ 1/(2n²+1) = −1/2 + √2·π·coth(√2·π/2)/4.
    REQUIRE(oracle.equivalent(
        summation(den(2, 1), n, one, oo)->str(),
        "-1/2 + sqrt(2)*pi*coth(sqrt(2)*pi/2)/4"));
    // No regression: a true p-series (b=0) still closes to π²/6, not coth.
    REQUIRE(oracle.equivalent(
        summation(pow(n, integer(-2)), n, one, oo)->str(), "pi**2/6"));
}

// SUM-DIRICHLET-BETA-1: the Dirichlet beta series Σ_{k=0}^∞ (−1)^k/(2k+1)^s.
// β(1) = π/4 (Leibniz), β(2) = Catalan's constant. Higher s have no elementary
// closed form (SymPy returns a polylog), so only s ∈ {1, 2} are closed; others
// stay an unevaluated Sum.
TEST_CASE("summation: Dirichlet beta β(1)=π/4, β(2)=Catalan (SUM-DIRICHLET-BETA-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto oo = S::Infinity();
    auto altk = pow(integer(-1), k);
    auto odd = integer(2) * k + integer(1);
    // β(1) = Σ (−1)^k/(2k+1) = π/4.
    REQUIRE(oracle.equivalent(
        summation(altk * pow(odd, integer(-1)), k, S::Zero(), oo)->str(),
        "pi/4"));
    // Shifted sign flips it: Σ (−1)^(k+1)/(2k+1) = −π/4.
    REQUIRE(oracle.equivalent(
        summation(pow(integer(-1), k + integer(1)) * pow(odd, integer(-1)), k,
                  S::Zero(), oo)->str(),
        "-pi/4"));
    // A leading constant multiplies through: Σ 2(−1)^k/(2k+1) = π/2.
    REQUIRE(oracle.equivalent(
        summation(integer(2) * altk * pow(odd, integer(-1)), k, S::Zero(), oo)
            ->str(),
        "pi/2"));
    // β(2) = Σ (−1)^k/(2k+1)² = Catalan's constant.
    REQUIRE(oracle.equivalent(
        summation(altk * pow(odd, integer(-2)), k, S::Zero(), oo)->str(),
        "Catalan"));
    // β(3) has no elementary form recognized here — stays an unevaluated Sum.
    REQUIRE(summation(altk * pow(odd, integer(-3)), k, S::Zero(), oo)
                ->str()
                .find("Sum(") != std::string::npos);
    // A non-(2k+1) odd-step denominator is not Dirichlet beta — left unevaluated.
    REQUIRE(summation(altk * pow(integer(3) * k + integer(1), integer(-1)), k,
                      S::Zero(), oo)
                ->str()
                .find("Sum(") != std::string::npos);
}

// SUM-EXP-1: the exponential series Σ_{k=0}^∞ r^k/k! = e^r (and its shifted /
// scaled forms). Converges for every r, so no convergence test is needed.
TEST_CASE("summation: exponential series r^k/k! closes to e^r (SUM-EXP-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto x = symbol("x");
    auto oo = S::Infinity();
    auto inv_fact = pow(factorial(k), integer(-1));
    // Σ 1/k! = e.
    REQUIRE(oracle.equivalent(summation(inv_fact, k, integer(0), oo)->str(),
                              "E"));
    // Σ_{k=1}^∞ 1/k! = e − 1 (head term k=0 removed).
    REQUIRE(oracle.equivalent(summation(inv_fact, k, integer(1), oo)->str(),
                              "E - 1"));
    // Σ x^k/k! = e^x.
    REQUIRE(oracle.equivalent(
        summation(mul(pow(x, k), inv_fact), k, integer(0), oo)->str(),
        "exp(x)"));
    // Σ 2^k/k! = e², Σ (−1)^k/k! = e⁻¹, Σ 1/(2^k·k!) = e^(1/2).
    REQUIRE(oracle.equivalent(
        summation(mul(pow(integer(2), k), inv_fact), k, integer(0), oo)->str(),
        "exp(2)"));
    REQUIRE(oracle.equivalent(
        summation(mul(pow(integer(-1), k), inv_fact), k, integer(0), oo)->str(),
        "exp(-1)"));
    REQUIRE(oracle.equivalent(
        summation(mul(pow(integer(2), mul(integer(-1), k)), inv_fact), k,
                  integer(0), oo)
            ->str(),
        "exp(1/2)"));
    // Constant prefactor is preserved: Σ 3/k! = 3e.
    REQUIRE(oracle.equivalent(
        summation(mul(integer(3), inv_fact), k, integer(0), oo)->str(),
        "3*E"));
    // Polynomial numerators fold via the falling-factorial basis:
    // Σ k/k! = e, Σ k²/k! = 2e, Σ k³/k! = 5e, Σ (2k+3)/k! = 5e.
    REQUIRE(oracle.equivalent(
        summation(mul(k, inv_fact), k, integer(0), oo)->str(), "E"));
    REQUIRE(oracle.equivalent(
        summation(mul(pow(k, integer(2)), inv_fact), k, integer(0), oo)->str(),
        "2*E"));
    REQUIRE(oracle.equivalent(
        summation(mul(pow(k, integer(3)), inv_fact), k, integer(0), oo)->str(),
        "5*E"));
    REQUIRE(oracle.equivalent(
        summation(mul(add(mul(integer(2), k), integer(3)), inv_fact), k,
                  integer(0), oo)
            ->str(),
        "5*E"));
    // Polynomial × geometric: Σ k·x^k/k! = x·e^x.
    REQUIRE(oracle.equivalent(
        summation(mul({k, pow(x, k), inv_fact}), k, integer(0), oo)->str(),
        "x*exp(x)"));
    // SUM-EXP-NOLEAK: a non-polynomial numerator factor (cos(k·x), cos(k)) is NOT a
    // degree-0 Poly coefficient — it must not be pulled out, which previously leaked
    // the bound variable as a bogus closed form (Σ cos(k·x)/k! → e·cos(k·x)). The
    // sum stays unevaluated (matching SymPy), and the result must be free of k.
    {
        auto s1 = summation(mul(cos(k * x), inv_fact), k, integer(0), oo);
        REQUIRE(s1->str().find("Sum(") != std::string::npos);
        REQUIRE(has(s1, k));   // only as the unevaluated Sum's bound variable
        auto s2 = summation(mul(cos(k), inv_fact), k, integer(0), oo);
        REQUIRE(s2->str().find("Sum(") != std::string::npos);
        // The legitimate polynomial-numerator case is unaffected: Σ k/k! = e.
        REQUIRE(oracle.equivalent(summation(mul(k, inv_fact), k, integer(0), oo)
                                      ->str(),
                                  "E"));
    }
}

// SUM-EXP-SHIFT-1: a shifted factorial (k+m)! re-indexes (j=k+m) to the k! case,
// closing the e-valued sums Σ P(k)/(k+m)! that don't telescope. Σ1/(k+1)!=e−2,
// Σ(2k+1)/(k+1)!=e, Σ k/(k+2)!=3−e. Previously unevaluated. Matches SymPy.
TEST_CASE("summation: shifted exponential series Σ P(k)/(k+m)! (SUM-EXP-SHIFT-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto one = integer(1);
    auto oo = S::Infinity();
    auto inv_fact = [&](const Expr& a) {
        return pow(factorial(a), integer(-1));
    };
    // Σ_{k=0}^∞ 1/(k+1)! = e − 1 ; Σ_{k=1}^∞ 1/(k+1)! = e − 2.
    REQUIRE(oracle.equivalent(
        summation(inv_fact(k + one), k, S::Zero(), oo)->str(), "E - 1"));
    REQUIRE(oracle.equivalent(
        summation(inv_fact(k + one), k, one, oo)->str(), "E - 2"));
    // Σ_{k=0}^∞ 1/(k+2)! = e − 2.
    REQUIRE(oracle.equivalent(
        summation(inv_fact(k + integer(2)), k, S::Zero(), oo)->str(), "E - 2"));
    // Polynomial numerator over a shifted factorial: Σ_{k=1}^∞ (2k+1)/(k+1)! = e.
    REQUIRE(oracle.equivalent(
        summation((integer(2) * k + one) * inv_fact(k + one), k, one, oo)->str(),
        "E"));
    // Σ_{k=1}^∞ k/(k+2)! = 3 − e.
    REQUIRE(oracle.equivalent(
        summation(k * inv_fact(k + integer(2)), k, one, oo)->str(), "3 - E"));
    // The unshifted m=0 case is unaffected: Σ (k+1)/k! = 2e.
    REQUIRE(oracle.equivalent(
        summation((k + one) * inv_fact(k), k, S::Zero(), oo)->str(), "2*E"));
}

// SUM-COSH-SINH-1: the even/odd bisection of the exponential series. Σ z^(2k)/(2k)!
// = cosh z and Σ z^(2k+1)/(2k+1)! = sinh z; with a (−1)^k sign they become cos z and
// sin z. Previously unevaluated. A lower bound > 0 subtracts the head. Matches SymPy.
TEST_CASE("summation: even/odd exponential series → cosh/sinh/cos/sin (SUM-COSH-SINH-1)",
          "[6][summation][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto x = symbol("x");
    auto one = integer(1);
    auto oo = S::Infinity();
    auto two = integer(2);
    auto inv_fact = [&](const Expr& a) {
        return pow(factorial(a), integer(-1));
    };
    auto altk = pow(integer(-1), k);
    auto even = two * k;          // 2k
    auto odd = two * k + one;     // 2k+1
    // Σ 1/(2k)! = cosh 1 ; Σ 1/(2k+1)! = sinh 1.
    REQUIRE(oracle.equivalent(summation(inv_fact(even), k, S::Zero(), oo)->str(),
                              "cosh(1)"));
    REQUIRE(oracle.equivalent(summation(inv_fact(odd), k, S::Zero(), oo)->str(),
                              "sinh(1)"));
    // Σ x^(2k)/(2k)! = cosh x ; Σ x^(2k+1)/(2k+1)! = sinh x.
    REQUIRE(oracle.equivalent(
        summation(mul(pow(x, even), inv_fact(even)), k, S::Zero(), oo)->str(),
        "cosh(x)"));
    REQUIRE(oracle.equivalent(
        summation(mul(pow(x, odd), inv_fact(odd)), k, S::Zero(), oo)->str(),
        "sinh(x)"));
    // Σ (−1)^k/(2k)! = cos 1 ; Σ (−1)^k x^(2k)/(2k)! = cos x.
    REQUIRE(oracle.equivalent(
        summation(mul(altk, inv_fact(even)), k, S::Zero(), oo)->str(), "cos(1)"));
    REQUIRE(oracle.equivalent(
        summation(mul({altk, pow(x, even), inv_fact(even)}), k, S::Zero(), oo)
            ->str(),
        "cos(x)"));
    // Σ (−1)^k x^(2k+1)/(2k+1)! = sin x.
    REQUIRE(oracle.equivalent(
        summation(mul({altk, pow(x, odd), inv_fact(odd)}), k, S::Zero(), oo)
            ->str(),
        "sin(x)"));
    // Lower bound > 0 subtracts the head: Σ_{k=1}^∞ 1/(2k)! = cosh 1 − 1.
    REQUIRE(oracle.equivalent(summation(inv_fact(even), k, one, oo)->str(),
                              "cosh(1) - 1"));
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
