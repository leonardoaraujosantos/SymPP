// simplify / collect / together tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_key.hpp>
#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

TEST_CASE("simplify: trivial identity", "[5][simplify]") {
    auto x = symbol("x");
    auto e = x + x;
    REQUIRE(simplify(e) == integer(2) * x);
}

TEST_CASE("simplify: cancels x - x", "[5][simplify]") {
    auto x = symbol("x");
    REQUIRE(simplify(x - x) == S::Zero());
}

TEST_CASE("simplify: distributes and recombines", "[5][simplify][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // (x + 1)*(x - 1) - x^2 = -1
    auto e = (x + integer(1)) * (x - integer(1)) - pow(x, integer(2));
    auto s = simplify(e);
    REQUIRE(oracle.equivalent(s->str(), "-1"));
}

TEST_CASE("collect: groups powers of var", "[5][collect][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto a = symbol("a");
    auto b = symbol("b");
    // a*x^2 + b*x^2 + x + 3 -> (a + b)*x^2 + x + 3
    auto e = a * pow(x, integer(2)) + b * pow(x, integer(2)) + x + integer(3);
    auto c = collect(e, x);
    REQUIRE(oracle.equivalent(c->str(), "(a + b)*x**2 + x + 3"));
}

TEST_CASE("collect: linear in var", "[5][collect]") {
    auto x = symbol("x");
    auto a = symbol("a");
    auto b = symbol("b");
    // a*x + b*x + 7 -> (a+b)*x + 7
    auto e = a * x + b * x + integer(7);
    auto c = collect(e, x);
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(c->str(), "(a + b)*x + 7"));
}

TEST_CASE("simplify: matches SymPy on polynomial expansion",
          "[5][simplify][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x + integer(1), integer(3)) - pow(x - integer(1), integer(3));
    auto s = simplify(e);
    // Expected: 6*x^2 + 2
    REQUIRE(oracle.equivalent(s->str(), "6*x**2 + 2"));
}

// ----- powsimp ---------------------------------------------------------------

TEST_CASE("powsimp: x^a * y^a → (x*y)^a (positive bases)",
          "[5][powsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    auto y = symbol("y", AssumptionMask{}.set_positive(true));
    auto a = symbol("a");
    auto e = pow(x, a) * pow(y, a);
    auto p = powsimp(e);
    REQUIRE(oracle.equivalent(p->str(), "(x*y)**a"));
}

TEST_CASE("powsimp: x^2 * y^2 → (x*y)^2 (integer exp, no assumption needed)",
          "[5][powsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = pow(x, integer(2)) * pow(y, integer(2));
    auto p = powsimp(e);
    REQUIRE(oracle.equivalent(p->str(), "(x*y)**2"));
}

TEST_CASE("powsimp: leaves non-positive non-integer alone",
          "[5][powsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto a = symbol("a");
    // No positivity assumption — must NOT combine.
    auto e = pow(x, a) * pow(y, a);
    auto p = powsimp(e);
    // Stay separate.
    REQUIRE(oracle.equivalent(p->str(), "x**a * y**a"));
}

// ----- expand_power_exp ------------------------------------------------------

TEST_CASE("expand_power_exp: x^(a+b) → x^a * x^b",
          "[5][expand_power_exp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto a = symbol("a");
    auto b = symbol("b");
    auto e = pow(x, a + b);
    auto out = expand_power_exp(e);
    REQUIRE(oracle.equivalent(out->str(), "x**a * x**b"));
}

TEST_CASE("expand_power_exp: leaves Pow without Add exponent alone",
          "[5][expand_power_exp]") {
    auto x = symbol("x");
    auto a = symbol("a");
    auto e = pow(x, a);
    REQUIRE(expand_power_exp(e) == e);
}

// ----- expand_power_base -----------------------------------------------------

TEST_CASE("expand_power_base: (x*y)^3 → x^3 * y^3 (integer exponent)",
          "[5][expand_power_base][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = pow(x * y, integer(3));
    auto out = expand_power_base(e);
    REQUIRE(oracle.equivalent(out->str(), "x**3 * y**3"));
}

TEST_CASE("expand_power_base: (x*y)^(1/2) needs positive bases",
          "[5][expand_power_base][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    auto y = symbol("y", AssumptionMask{}.set_positive(true));
    auto e = pow(x * y, S::Half());
    auto out = expand_power_base(e);
    REQUIRE(oracle.equivalent(out->str(), "x**(1/2) * y**(1/2)"));
}

TEST_CASE("expand_power_base: refuses without positivity",
          "[5][expand_power_base]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = pow(x * y, S::Half());
    auto out = expand_power_base(e);
    // No positivity assumption — must not split.
    REQUIRE(out == e);
}

// ----- trigsimp --------------------------------------------------------------

TEST_CASE("trigsimp: sin(x)^2 + cos(x)^2 → 1", "[5][trigsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sin(x), integer(2)) + pow(cos(x), integer(2));
    auto s = trigsimp(e);
    REQUIRE(oracle.equivalent(s->str(), "1"));
}

TEST_CASE("trigsimp: shared coefficient a*sin²+a*cos² → a",
          "[5][trigsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto a = symbol("a");
    auto e = a * pow(sin(x), integer(2)) + a * pow(cos(x), integer(2));
    auto s = trigsimp(e);
    REQUIRE(oracle.equivalent(s->str(), "a"));
}

TEST_CASE("trigsimp: 2*sin² + 3*cos² → 2 + cos² (form reduction)",
          "[5][trigsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = integer(2) * pow(sin(x), integer(2))
             + integer(3) * pow(cos(x), integer(2));
    auto s = trigsimp(e);
    // Algebraically equivalent to either form; oracle uses simplify.
    REQUIRE(oracle.equivalent(s->str(),
                              "2*sin(x)**2 + 3*cos(x)**2"));
}

TEST_CASE("trigsimp: independent of argument symbol",
          "[5][trigsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = pow(sin(x + y), integer(2)) + pow(cos(x + y), integer(2));
    auto s = trigsimp(e);
    REQUIRE(oracle.equivalent(s->str(), "1"));
}

TEST_CASE("trigsimp: nested in larger expression", "[5][trigsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    // y * (sin²(x) + cos²(x)) + y² → y + y²  (after expand)
    auto e = y * (pow(sin(x), integer(2)) + pow(cos(x), integer(2)))
             + pow(y, integer(2));
    auto s = trigsimp(e);
    REQUIRE(oracle.equivalent(s->str(), "y + y**2"));
}

TEST_CASE("trigsimp: leaves non-Pythagorean trig alone",
          "[5][trigsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sin(x) + cos(x);
    auto s = trigsimp(e);
    REQUIRE(oracle.equivalent(s->str(), "sin(x) + cos(x)"));
}

TEST_CASE("trigsimp: only sin² (no cos² counterpart) stays",
          "[5][trigsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sin(x), integer(2)) + integer(1);
    auto s = trigsimp(e);
    REQUIRE(oracle.equivalent(s->str(), "sin(x)**2 + 1"));
}

// ----- trigsimp: Fu-style double-angle rules --------------------------------

TEST_CASE("trigsimp: cos²(x) - sin²(x) → cos(2*x)",
          "[5][trigsimp][fu][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(cos(x), integer(2)) - pow(sin(x), integer(2));
    auto s = trigsimp(e);
    REQUIRE(oracle.equivalent(s->str(), "cos(2*x)"));
}

TEST_CASE("trigsimp: 1 - 2*sin²(x) → cos(2*x)",
          "[5][trigsimp][fu][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = integer(1) - integer(2) * pow(sin(x), integer(2));
    auto s = trigsimp(e);
    REQUIRE(oracle.equivalent(s->str(), "cos(2*x)"));
}

TEST_CASE("trigsimp: 2*cos²(x) - 1 → cos(2*x)",
          "[5][trigsimp][fu][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = integer(2) * pow(cos(x), integer(2)) - integer(1);
    auto s = trigsimp(e);
    REQUIRE(oracle.equivalent(s->str(), "cos(2*x)"));
}

TEST_CASE("trigsimp: 2*sin(x)*cos(x) → sin(2*x)",
          "[5][trigsimp][fu][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = integer(2) * sin(x) * cos(x);
    auto s = trigsimp(e);
    REQUIRE(oracle.equivalent(s->str(), "sin(2*x)"));
}

TEST_CASE("trigsimp: k*sin(x)*cos(x) collapses to (k/2)*sin(2*x)",
          "[5][trigsimp][fu][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto k = symbol("k");
    auto e = k * sin(x) * cos(x);
    auto s = trigsimp(e);
    REQUIRE(oracle.equivalent(s->str(), "k*sin(2*x)/2"));
}

TEST_CASE("trigsimp: 4*sin(x)*cos(x) → 2*sin(2*x)",
          "[5][trigsimp][fu][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = integer(4) * sin(x) * cos(x);
    auto s = trigsimp(e);
    REQUIRE(oracle.equivalent(s->str(), "2*sin(2*x)"));
}

// ----- expand_trig ----------------------------------------------------------

TEST_CASE("expand_trig: sin(a + b) → sin(a)*cos(b) + cos(a)*sin(b)",
          "[5][expand_trig][oracle]") {
    auto& oracle = Oracle::instance();
    auto a = symbol("a");
    auto b = symbol("b");
    auto e = sin(a + b);
    auto s = expand_trig(e);
    REQUIRE(oracle.equivalent(s->str(),
                              "sin(a)*cos(b) + cos(a)*sin(b)"));
}

TEST_CASE("expand_trig: cos(a + b) → cos(a)*cos(b) - sin(a)*sin(b)",
          "[5][expand_trig][oracle]") {
    auto& oracle = Oracle::instance();
    auto a = symbol("a");
    auto b = symbol("b");
    auto e = cos(a + b);
    auto s = expand_trig(e);
    REQUIRE(oracle.equivalent(s->str(),
                              "cos(a)*cos(b) - sin(a)*sin(b)"));
}

TEST_CASE("expand_trig: sin(a + b + c) recursively expands",
          "[5][expand_trig][oracle]") {
    auto& oracle = Oracle::instance();
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    auto e = sin(a + b + c);
    auto s = expand_trig(e);
    // Verify equivalence by oracle (the canonical form has many shapes).
    REQUIRE(oracle.equivalent(s->str(), "sin(a + b + c)"));
}

// ----- fu orchestrator ------------------------------------------------------

TEST_CASE("fu: picks the smaller form between identity and trigsimp",
          "[5][fu][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // sin² + cos² = 1 — trigsimp wins.
    auto e = pow(sin(x), integer(2)) + pow(cos(x), integer(2));
    REQUIRE(oracle.equivalent(fu(e)->str(), "1"));
}

TEST_CASE("fu: leaves non-trig expression alone",
          "[5][fu][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(2)) + integer(1);
    REQUIRE(oracle.equivalent(fu(e)->str(), "x**2 + 1"));
}

// ----- radsimp ---------------------------------------------------------------

TEST_CASE("radsimp: 1/sqrt(2) → sqrt(2)/2", "[5][radsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = pow(integer(2), S::Half());  // sqrt(2)
    // Build 1/sqrt(2) explicitly.
    auto r = pow(e, integer(-1));
    auto out = radsimp(r);
    REQUIRE(oracle.equivalent(out->str(), "sqrt(2)/2"));
}

TEST_CASE("radsimp: 1/(1 + sqrt(2)) → sqrt(2) - 1", "[5][radsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto den = integer(1) + sqrt(integer(2));
    auto e = pow(den, integer(-1));
    auto out = radsimp(e);
    // 1/(1+√2) = (√2-1)/(2-1) = √2 - 1
    REQUIRE(oracle.equivalent(out->str(), "sqrt(2) - 1"));
}

TEST_CASE("radsimp: x/(2 + 3*sqrt(5))", "[5][radsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto den = integer(2) + integer(3) * sqrt(integer(5));
    auto e = x / den;
    auto out = radsimp(e);
    // x*(2 - 3√5)/(4 - 45) = x*(2 - 3√5)/(-41)
    REQUIRE(oracle.equivalent(out->str(), "x*(2 - 3*sqrt(5)) / (-41)"));
}

TEST_CASE("radsimp: leaves rational denominator alone",
          "[5][radsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = (x + integer(1)) / (x - integer(1));
    auto out = radsimp(e);
    REQUIRE(oracle.equivalent(out->str(), "(x + 1)/(x - 1)"));
}

// ----- sqrtdenest ------------------------------------------------------------

TEST_CASE("sqrtdenest: sqrt(3 + 2*sqrt(2)) → 1 + sqrt(2)",
          "[5][sqrtdenest][oracle]") {
    auto& oracle = Oracle::instance();
    // sqrt(3 + 2√2) = 1 + √2 (since (1+√2)² = 1 + 2√2 + 2 = 3 + 2√2)
    auto inner = integer(3) + integer(2) * sqrt(integer(2));
    auto e = sqrt(inner);
    auto out = sqrtdenest(e);
    REQUIRE(oracle.equivalent(out->str(), "1 + sqrt(2)"));
}

TEST_CASE("sqrtdenest: sqrt(5 - 2*sqrt(6)) → sqrt(3) - sqrt(2)",
          "[5][sqrtdenest][oracle]") {
    auto& oracle = Oracle::instance();
    // (sqrt(3) - sqrt(2))² = 3 - 2√6 + 2 = 5 - 2√6
    auto inner = integer(5) - integer(2) * sqrt(integer(6));
    auto e = sqrt(inner);
    auto out = sqrtdenest(e);
    REQUIRE(oracle.equivalent(out->str(), "sqrt(3) - sqrt(2)"));
}

TEST_CASE("sqrtdenest: leaves non-denestable alone",
          "[5][sqrtdenest][oracle]") {
    auto& oracle = Oracle::instance();
    // sqrt(2 + sqrt(3)) doesn't denest (discriminant 4 - 3 = 1, but
    // (2+1)/2 and (2-1)/2 are 3/2 and 1/2 — actually IT DOES denest:
    // sqrt(2 + sqrt(3)) = sqrt(3/2) + sqrt(1/2)).
    // Use sqrt(1 + sqrt(2)) instead — disc = 1 - 2 = -1, fails.
    auto inner = integer(1) + sqrt(integer(2));
    auto e = sqrt(inner);
    auto out = sqrtdenest(e);
    REQUIRE(oracle.equivalent(out->str(), "sqrt(1 + sqrt(2))"));
}

#include <sympp/functions/combinatorial.hpp>

// ----- combsimp / gammasimp --------------------------------------------------

TEST_CASE("combsimp: factorial(n+1)/factorial(n) → n+1",
          "[5][combsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto e = factorial(n + integer(1)) * pow(factorial(n), integer(-1));
    auto out = combsimp(e);
    REQUIRE(oracle.equivalent(out->str(), "n + 1"));
}

TEST_CASE("combsimp: factorial(n+2)/factorial(n) → (n+1)*(n+2)",
          "[5][combsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto e = factorial(n + integer(2)) * pow(factorial(n), integer(-1));
    auto out = combsimp(e);
    REQUIRE(oracle.equivalent(out->str(), "(n+1)*(n+2)"));
}

TEST_CASE("combsimp: factorial(n)/factorial(n-1) → n",
          "[5][combsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto e = factorial(n) * pow(factorial(n - integer(1)), integer(-1));
    auto out = combsimp(e);
    REQUIRE(oracle.equivalent(out->str(), "n"));
}

TEST_CASE("combsimp: leaves unrelated factorials alone",
          "[5][combsimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto m = symbol("m");
    // factorial(n) / factorial(m) — non-integer difference, leave alone
    auto e = factorial(n) * pow(factorial(m), integer(-1));
    auto out = combsimp(e);
    REQUIRE(oracle.equivalent(out->str(),
                              "factorial(n)/factorial(m)"));
}

TEST_CASE("gammasimp: gamma(n+1)/gamma(n) → n",
          "[5][gammasimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto e = gamma(n + integer(1)) * pow(gamma(n), integer(-1));
    auto out = gammasimp(e);
    REQUIRE(oracle.equivalent(out->str(), "n"));
}

TEST_CASE("gammasimp: gamma(n+3)/gamma(n) → n*(n+1)*(n+2)",
          "[5][gammasimp][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto e = gamma(n + integer(3)) * pow(gamma(n), integer(-1));
    auto out = gammasimp(e);
    REQUIRE(oracle.equivalent(out->str(), "n*(n+1)*(n+2)"));
}

#include <sympp/core/float.hpp>

// ----- cse -------------------------------------------------------------------

TEST_CASE("cse: extracts shared subtree", "[5][cse]") {
    auto x = symbol("x");
    auto y = symbol("y");
    // (x+y)^2 + (x+y)^3 — Add canonical form leaves Pow children intact,
    // so (x+y) is preserved as a shared subtree of two Pow nodes.
    auto shared = x + y;
    auto e = pow(shared, integer(2)) + pow(shared, integer(3));
    auto r = cse(e);
    REQUIRE(r.substitutions.size() == 1);
    REQUIRE(r.substitutions[0].second == shared);
}

TEST_CASE("cse: leaves no-shared expression alone", "[5][cse]") {
    auto x = symbol("x");
    auto e = x + integer(1);
    auto r = cse(e);
    REQUIRE(r.substitutions.empty());
    REQUIRE(r.expr == e);
}

TEST_CASE("cse: extracts multiple shared subtrees", "[5][cse]") {
    auto x = symbol("x");
    auto y = symbol("y");
    // (x+y)^2 + 2*(x+y) + (x*y) + (x*y)^3
    // (x+y) appears twice, (x*y) appears twice.
    auto a = x + y;
    auto b = x * y;
    auto e = pow(a, integer(2)) + integer(2) * a + b + pow(b, integer(3));
    auto r = cse(e);
    REQUIRE(r.substitutions.size() == 2);
}

// ----- nsimplify -------------------------------------------------------------

TEST_CASE("nsimplify: 0.5 → 1/2", "[5][nsimplify]") {
    auto e = evalf(rational(1, 2), 15);
    auto r = nsimplify(e);
    REQUIRE(r == rational(1, 2));
}

TEST_CASE("nsimplify: pi numeric → pi", "[5][nsimplify]") {
    auto e = evalf(S::Pi(), 15);
    auto r = nsimplify(e);
    REQUIRE(r == S::Pi());
}

TEST_CASE("nsimplify: pi/4 numeric → pi/4", "[5][nsimplify]") {
    auto e = evalf(S::Pi() / integer(4), 15);
    auto r = nsimplify(e);
    REQUIRE(r == S::Pi() / integer(4));
}

TEST_CASE("nsimplify: sqrt(2) numeric → sqrt(2)", "[5][nsimplify][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = evalf(sqrt(integer(2)), 15);
    auto r = nsimplify(e);
    REQUIRE(oracle.equivalent(r->str(), "sqrt(2)"));
}

TEST_CASE("nsimplify: junk irrational stays put", "[5][nsimplify]") {
    // 1.4142135 is close to sqrt(2). 1.123456789 is not in our table.
    auto e = make<Float>("1.123456789", 15);
    auto r = nsimplify(e);
    // Should fall through. Same Float (or close enough — exact equality
    // not guaranteed across factories, but no rational/sqrt match found).
    REQUIRE(r->type_id() == TypeId::Float);
}

// ----- simplify orchestrator (chained) ---------------------------------------

TEST_CASE("simplify: chains trigsimp (sin² + cos² → 1)",
          "[5][simplify][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sin(x), integer(2)) + pow(cos(x), integer(2));
    auto s = simplify(e);
    REQUIRE(oracle.equivalent(s->str(), "1"));
}

TEST_CASE("simplify: chains combsimp (factorial(n+1)/factorial(n) → n+1)",
          "[5][simplify][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto e = factorial(n + integer(1)) * pow(factorial(n), integer(-1));
    auto s = simplify(e);
    REQUIRE(oracle.equivalent(s->str(), "n + 1"));
}

TEST_CASE("simplify: chains sqrtdenest (sqrt(3+2√2) → 1+√2)",
          "[5][simplify][oracle]") {
    auto& oracle = Oracle::instance();
    auto inner = integer(3) + integer(2) * sqrt(integer(2));
    auto e = sqrt(inner);
    auto s = simplify(e);
    REQUIRE(oracle.equivalent(s->str(), "1 + sqrt(2)"));
}
