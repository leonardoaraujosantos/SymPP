// MATLAB facade — PDE tests.
//
// Validates pdsolve / pdsolve_heat / pdsolve_wave wrappers.

#include <catch2/catch_test_macros.hpp>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/matlab/matlab.hpp>

using namespace sympp;

TEST_CASE("matlab::pdsolve constant-coefficient: 2 u_x + 3 u_y = 0",
          "[15m][matlab][pdsolve]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto sol = matlab::pdsolve(integer(2), integer(3), S::Zero(), x, y);
    // u(x, y) = F(3x - 2y) — opaque function call.
    REQUIRE(sol->type_id() == TypeId::Function);
}

TEST_CASE("matlab::pdsolve inhomogeneous: u_x + u_y = 5",
          "[15m][matlab][pdsolve]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto sol = matlab::pdsolve(integer(1), integer(1), integer(5), x, y);
    // u(x, y) = 5x + F(x - y) — Add of particular + homogeneous.
    REQUIRE(sol->type_id() == TypeId::Add);
}

TEST_CASE("matlab::pdsolve_heat satisfies u_t = k u_xx",
          "[15m][matlab][pdsolve_heat]") {
    auto x = symbol("x");
    auto t = symbol("t");
    auto k = symbol("k");
    auto lambda = symbol("lambda");
    auto u = matlab::pdsolve_heat(k, lambda, x, t);
    auto u_t = sympp::diff(u, t);
    auto u_xx = sympp::diff(sympp::diff(u, x), x);
    auto residual = matlab::simplify(u_t - k * u_xx);
    REQUIRE(residual == S::Zero());
}

TEST_CASE("matlab::pdsolve_wave d'Alembert form",
          "[15m][matlab][pdsolve_wave]") {
    auto x = symbol("x");
    auto t = symbol("t");
    auto c = symbol("c");
    auto u = matlab::pdsolve_wave(c, x, t);
    REQUIRE(u->type_id() == TypeId::Add);
}
