// autowrap — compile expressions to native code via the system C compiler and
// check the result against lambdify / closed-form values.

#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <sympp/codegen/autowrap.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/lambdify.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/trigonometric.hpp>

using namespace sympp;
using namespace sympp::codegen;
using Catch::Matchers::WithinAbs;

TEST_CASE("autowrap: single-variable polynomial matches lambdify", "[autowrap]") {
    if (!autowrap_available()) {
        SKIP("no C compiler available");
    }
    auto x = symbol("x");
    auto body = pow(x, integer(2)) + integer(1);  // x² + 1

    auto f = autowrap(x, body);
    auto ref = lambdify(x, body);

    for (double xv : {-3.0, -0.5, 0.0, 1.0, 2.5, 7.3}) {
        REQUIRE_THAT(f({xv}), WithinAbs(ref({xv}), 1e-9));
        REQUIRE_THAT(f({xv}), WithinAbs(xv * xv + 1.0, 1e-9));
    }
}

TEST_CASE("autowrap: transcendental sin(x)*exp(x)", "[autowrap]") {
    if (!autowrap_available()) {
        SKIP("no C compiler available");
    }
    auto x = symbol("x");
    auto body = sin(x) * exp(x);

    auto f = autowrap(x, body);
    auto ref = lambdify(x, body);

    for (double xv : {-2.0, -0.3, 0.0, 0.7, 1.9, 3.1}) {
        REQUIRE_THAT(f({xv}), WithinAbs(ref({xv}), 1e-9));
        REQUIRE_THAT(f({xv}), WithinAbs(std::sin(xv) * std::exp(xv), 1e-9));
    }
}

TEST_CASE("autowrap: multivariate x*y - y", "[autowrap]") {
    if (!autowrap_available()) {
        SKIP("no C compiler available");
    }
    auto x = symbol("x"), y = symbol("y");
    auto body = x * y - y;

    auto f = autowrap({x, y}, body);
    auto ref = lambdify({x, y}, body);

    const std::vector<std::vector<double>> pts = {
        {0.0, 0.0}, {1.0, 2.0}, {-3.0, 4.0}, {2.5, -1.5}, {10.0, 0.25}};
    for (const auto& p : pts) {
        REQUIRE_THAT(f(p), WithinAbs(ref(p), 1e-9));
        REQUIRE_THAT(f(p), WithinAbs(p[0] * p[1] - p[1], 1e-9));
    }
}
