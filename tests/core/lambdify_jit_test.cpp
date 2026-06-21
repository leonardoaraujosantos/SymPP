// lambdify LLVM-JIT backend — native-compiled evaluation, cross-checked against
// the portable interpreter and known values. Skips gracefully when SymPP was
// built without LLVM.

#include <cmath>
#include <stdexcept>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/lambdify.hpp>
#include <sympp/core/lambdify_llvm.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>

using namespace sympp;
using Catch::Matchers::WithinAbs;

TEST_CASE("lambdify_jit availability and fallback", "[lambdify][jit]") {
    auto x = symbol("x");
    if (!lambdify_jit_available()) {
        REQUIRE_THROWS_AS(lambdify_jit(x, x), std::runtime_error);
        SUCCEED("built without LLVM — JIT correctly unavailable");
        return;
    }
    auto f = lambdify_jit(x, x);
    REQUIRE_THAT(f({3.0}), WithinAbs(3.0, 1e-12));
}

TEST_CASE("lambdify_jit matches the interpreter backend", "[lambdify][jit]") {
    if (!lambdify_jit_available()) return;
    auto x = symbol("x"), y = symbol("y");
    struct Case { std::vector<Expr> vars; Expr body; std::vector<double> at; };
    std::vector<Case> cases = {
        {{x}, pow(x, integer(3)) + integer(2) * x, {2.0}},
        {{x}, sin(x) * exp(x) + cos(x), {1.3}},
        {{x}, sqrt(x) + log(x), {2.5}},
        {{x}, tan(x) + abs(x), {0.6}},
        {{x, y}, x * y + exp(y) + pow(x, integer(2)), {1.5, -0.7}},
        {{x}, integer(2) * S::Pi() * x + tanh(x), {0.9}},
    };
    for (const auto& c : cases) {
        auto jit = lambdify_jit(c.vars, c.body);
        auto interp = lambdify(c.vars, c.body);
        REQUIRE_THAT(jit(c.at), WithinAbs(interp(c.at), 1e-9));
    }
}

TEST_CASE("lambdify_jit known values", "[lambdify][jit]") {
    if (!lambdify_jit_available()) return;
    auto x = symbol("x");
    auto f = lambdify_jit(x, sin(x) + pow(x, integer(2)));
    REQUIRE_THAT(f({0.0}), WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(f({1.0}), WithinAbs(std::sin(1.0) + 1.0, 1e-12));
}

TEST_CASE("lambdify_jit rejects complex / arity errors", "[lambdify][jit]") {
    if (!lambdify_jit_available()) return;
    auto x = symbol("x");
    REQUIRE_THROWS_AS(lambdify_jit(x, x + S::I()), std::runtime_error);
    auto f = lambdify_jit(x, x);
    REQUIRE_THROWS_AS(f({1.0, 2.0}), std::runtime_error);
}
