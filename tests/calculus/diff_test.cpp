// diff() tests.
//
// References:
//   sympy/core/tests/test_diff.py
//   sympy/functions/elementary/tests/test_trigonometric.py — d/dx sin(x) = cos(x)

#include <catch2/catch_test_macros.hpp>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/derivative.hpp>
#include <sympp/core/undefined_function.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Atomic rules ----------------------------------------------------------

TEST_CASE("diff: x w.r.t. x is 1", "[6a][diff]") {
    auto x = symbol("x");
    REQUIRE(diff(x, x) == S::One());
}

TEST_CASE("diff: y w.r.t. x is 0", "[6a][diff]") {
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(diff(y, x) == S::Zero());
}

TEST_CASE("diff: number is 0", "[6a][diff]") {
    auto x = symbol("x");
    REQUIRE(diff(integer(7), x) == S::Zero());
    REQUIRE(diff(rational(3, 4), x) == S::Zero());
    REQUIRE(diff(S::Pi(), x) == S::Zero());
}

// ----- Polynomial rules ------------------------------------------------------

TEST_CASE("diff: x^2 = 2*x", "[6a][diff][pow]") {
    auto x = symbol("x");
    auto e = pow(x, integer(2));
    REQUIRE(diff(e, x) == integer(2) * x);
}

TEST_CASE("diff: x^3 = 3*x^2", "[6a][diff][pow]") {
    auto x = symbol("x");
    auto e = pow(x, integer(3));
    auto d = diff(e, x);
    auto expected = integer(3) * pow(x, integer(2));
    REQUIRE(d == expected);
}

TEST_CASE("diff: 2*x^3 + x = 6*x^2 + 1", "[6a][diff][add]") {
    auto x = symbol("x");
    auto e = integer(2) * pow(x, integer(3)) + x;
    auto d = diff(e, x);
    auto expected = integer(6) * pow(x, integer(2)) + integer(1);
    REQUIRE(d == expected);
}

TEST_CASE("diff: product rule (x * y) w.r.t. x = y", "[6a][diff][mul]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = x * y;
    REQUIRE(diff(e, x) == y);
}

TEST_CASE("diff: x*sin(x)", "[6a][diff][mul][trig]") {
    auto x = symbol("x");
    auto e = x * sin(x);
    // d/dx (x*sin(x)) = sin(x) + x*cos(x)
    auto d = diff(e, x);
    auto expected = sin(x) + x * cos(x);
    REQUIRE(d == expected);
}

// ----- Function chain rule ---------------------------------------------------

TEST_CASE("diff: sin(x) = cos(x)", "[6a][diff][trig]") {
    auto x = symbol("x");
    REQUIRE(diff(sin(x), x) == cos(x));
}

TEST_CASE("diff: cos(x) = -sin(x)", "[6a][diff][trig]") {
    auto x = symbol("x");
    REQUIRE(diff(cos(x), x) == mul(S::NegativeOne(), sin(x)));
}

TEST_CASE("diff: exp(x) = exp(x)", "[6a][diff][exp]") {
    auto x = symbol("x");
    REQUIRE(diff(exp(x), x) == exp(x));
}

TEST_CASE("diff: log(x) = 1/x", "[6a][diff][log]") {
    auto x = symbol("x");
    auto d = diff(log(x), x);
    auto expected = pow(x, S::NegativeOne());
    REQUIRE(d == expected);
}

// Regression (issue #13 / DIFF-1): diff(Abs(x)) returned 0 because Abs had
// no diff_arg override and fell through to the default. d/dx|x| = sign(x).
TEST_CASE("diff: Abs(x) = sign(x)", "[6a][diff][abs][regression]") {
    auto x = symbol("x");
    REQUIRE(diff(abs(x), x) == sign(x));
}

TEST_CASE("diff: Abs(2x+1) = 2*sign(2x+1) (chain rule)",
          "[6a][diff][abs][chain][regression]") {
    auto x = symbol("x");
    auto d = diff(abs(integer(2) * x + integer(1)), x);
    REQUIRE(d == integer(2) * sign(integer(2) * x + integer(1)));
}

TEST_CASE("diff: x*Abs(x) = x*sign(x) + Abs(x) (product rule)",
          "[6a][diff][abs][mul][regression]") {
    auto x = symbol("x");
    auto d = diff(x * abs(x), x);
    REQUIRE(d == x * sign(x) + abs(x));
}

// Chain rule
TEST_CASE("diff: sin(2x) = 2*cos(2x)", "[6a][diff][trig][chain]") {
    auto x = symbol("x");
    auto e = sin(integer(2) * x);
    auto d = diff(e, x);
    auto expected = integer(2) * cos(integer(2) * x);
    REQUIRE(d == expected);
}

TEST_CASE("diff: exp(x^2) = 2*x*exp(x^2)", "[6a][diff][exp][chain]") {
    auto x = symbol("x");
    auto e = exp(pow(x, integer(2)));
    auto d = diff(e, x);
    auto expected = integer(2) * x * exp(pow(x, integer(2)));
    REQUIRE(d == expected);
}

// Higher-order
TEST_CASE("diff: second derivative of x^4", "[6a][diff][higher-order]") {
    auto x = symbol("x");
    auto e = pow(x, integer(4));
    auto d2 = diff(e, x, 2);
    // 12*x^2
    auto expected = integer(12) * pow(x, integer(2));
    REQUIRE(d2 == expected);
}

TEST_CASE("diff: third derivative of sin(x) = -cos(x)",
          "[6a][diff][higher-order]") {
    auto x = symbol("x");
    auto d3 = diff(sin(x), x, 3);
    REQUIRE(d3 == mul(S::NegativeOne(), cos(x)));
}

// ----- Oracle parity --------------------------------------------------------

TEST_CASE("diff: matches SymPy on polynomial", "[6a][diff][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(3)) + integer(2) * pow(x, integer(2)) - integer(7) * x + integer(5);
    auto d = diff(e, x);
    REQUIRE(oracle.equivalent(d->str(), "3*x**2 + 4*x - 7"));
}

TEST_CASE("diff: matches SymPy on trig composition", "[6a][diff][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sin(cos(x));
    auto d = diff(e, x);
    REQUIRE(oracle.equivalent(d->str(), "-sin(x)*cos(cos(x))"));
}

TEST_CASE("diff: matches SymPy on log of polynomial",
          "[6a][diff][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = log(pow(x, integer(2)) + integer(1));
    auto d = diff(e, x);
    REQUIRE(oracle.equivalent(d->str(), "2*x/(x**2 + 1)"));
}

TEST_CASE("diff: tanh derivative", "[6a][diff][hyperbolic]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto d = diff(tanh(x), x);
    // 1 - tanh(x)^2
    REQUIRE(oracle.equivalent(d->str(), "1 - tanh(x)**2"));
}

// ----- Unevaluated Derivative for functions with no closed-form rule ---------
//
// DERIV-1: diff() must NOT silently collapse the derivative of an undefined /
// untabulated function to 0 (the old Function::diff_arg default). It now keeps
// an unevaluated Derivative node, exactly as SymPy does.

TEST_CASE("diff: undefined function keeps an unevaluated Derivative",
          "[6a][diff][derivative][regression]") {
    auto x = symbol("x");
    auto f = function_symbol("f");
    // d/dx f(x) = Derivative(f(x), x), NOT 0.
    auto d = diff(f(x), x);
    REQUIRE(d != S::Zero());
    REQUIRE(d->type_id() == TypeId::Derivative);
    REQUIRE(d->str() == "Derivative(f(x), x)");
}

TEST_CASE("diff: product / power / sum rules carry the Derivative node",
          "[6a][diff][derivative][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto g = function_symbol("g");
    auto f = function_symbol("f");
    // Product rule: d/dx (x*g(x)) = x*Derivative(g(x), x) + g(x).
    REQUIRE(oracle.equivalent(diff(x * g(x), x)->str(),
                              "x*Derivative(g(x), x) + g(x)"));
    // Power rule: d/dx f(x)**2 = 2*f(x)*Derivative(f(x), x).
    REQUIRE(oracle.equivalent(diff(pow(f(x), integer(2)), x)->str(),
                              "2*f(x)*Derivative(f(x), x)"));
    // A function of an independent variable still differentiates to 0.
    auto y = symbol("y");
    REQUIRE(diff(f(y), x) == S::Zero());
}

TEST_CASE("diff: second derivative of an undefined function bumps the order",
          "[6a][diff][derivative][regression]") {
    auto x = symbol("x");
    auto f = function_symbol("f");
    auto d2 = diff(f(x), x, 2);
    REQUIRE(d2->str() == "Derivative(f(x), (x, 2))");
}
