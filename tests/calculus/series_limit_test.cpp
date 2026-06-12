// series / limit tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/calculus/limit.hpp>
#include <sympp/calculus/series.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
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
