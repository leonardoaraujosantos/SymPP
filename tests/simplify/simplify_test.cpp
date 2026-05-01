// simplify / collect / together tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_key.hpp>
#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
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
