// lambdify — compiled double-precision evaluation, cross-checked against SymPy.

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/lambdify.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/special.hpp>
#include <sympp/functions/trigonometric.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using Catch::Matchers::WithinAbs;

TEST_CASE("lambdify: arithmetic and powers", "[lambdify]") {
    auto x = symbol("x");
    auto f = lambdify(x, pow(x, integer(3)) + integer(2) * x);  // x³ + 2x
    REQUIRE_THAT(f({2.0}), WithinAbs(12.0, 1e-12));
    REQUIRE_THAT(f({0.0}), WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(f({-1.0}), WithinAbs(-3.0, 1e-12));
}

TEST_CASE("lambdify: elementary functions match libm", "[lambdify]") {
    auto x = symbol("x");
    auto f = lambdify(x, sin(x) + cos(x));
    for (double xv : {0.0, 0.3, 1.7, -2.1}) {
        REQUIRE_THAT(f({xv}), WithinAbs(std::sin(xv) + std::cos(xv), 1e-12));
    }
    auto g = lambdify(x, exp(x) * sqrt(x));
    REQUIRE_THAT(g({2.0}), WithinAbs(std::exp(2.0) * std::sqrt(2.0), 1e-12));
}

TEST_CASE("lambdify: multivariate and constants", "[lambdify]") {
    auto x = symbol("x"), y = symbol("y");
    auto f = lambdify({x, y}, x * y + exp(y));
    REQUIRE_THAT(f({3.0, 0.0}), WithinAbs(1.0, 1e-12));  // 0 + e^0
    REQUIRE_THAT(f({2.0, 1.0}), WithinAbs(2.0 + std::exp(1.0), 1e-12));

    auto h = lambdify(x, integer(2) * S::Pi() * x);
    REQUIRE_THAT(h({1.0}), WithinAbs(2.0 * M_PI, 1e-12));
}

TEST_CASE("lambdify: cross-checked against SymPy", "[lambdify][oracle]") {
    auto& oracle = sympp::testing::Oracle::instance();
    auto x = symbol("x");
    // A compound real expression; compare lambdify to SymPy's evalf of the
    // numerically-substituted form at several sample points.
    auto body = sin(x) * exp(x / integer(2)) + log(x + integer(2));
    auto f = lambdify(x, body);
    for (double xv : {0.5, 1.3, 2.7, 4.0}) {
        auto sub = subs(body, x, float_value(xv));
        double expected = std::stod(oracle.evalf(sub->str(), 15));
        REQUIRE_THAT(f({xv}), WithinAbs(expected, 1e-9));
    }
}

TEST_CASE("lambdify: rejects what it cannot evaluate in ℝ", "[lambdify]") {
    auto x = symbol("x"), y = symbol("y");
    // Imaginary unit.
    REQUIRE_THROWS_AS(lambdify(x, x + S::I()), std::runtime_error);
    // Free symbol not in the variable list.
    REQUIRE_THROWS_AS(lambdify(x, x + y), std::runtime_error);
    // Unsupported special function.
    REQUIRE_THROWS_AS(lambdify(x, zeta(x)), std::runtime_error);
    // Wrong arity at call time.
    auto f = lambdify(x, x);
    REQUIRE_THROWS_AS(f({1.0, 2.0}), std::runtime_error);
}
