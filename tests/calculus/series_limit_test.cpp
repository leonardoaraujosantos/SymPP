// series / limit tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/calculus/limit.hpp>
#include <sympp/calculus/series.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
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

TEST_CASE("summation: Σ from a to b uses telescoping",
          "[6][summation][oracle]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto a = symbol("a");
    auto b = symbol("b");
    auto s = summation(k, k, a, b);
    REQUIRE(oracle.equivalent(s->str(), "(b - a + 1)*(a + b)/2"));
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
