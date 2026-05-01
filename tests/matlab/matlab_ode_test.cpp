// MATLAB facade — ODE tests.
//
// Validates dsolve overloads (first-order, second-order with
// auto-classification, IVP variants, system).

#include <catch2/catch_test_macros.hpp>

#include <utility>
#include <vector>

#include <sympp/core/function.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/matlab/matlab.hpp>
#include <sympp/matrices/matrix.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {

bool is_unevaluated_marker(const Expr& sol) {
    if (sol->type_id() != TypeId::Function) return false;
    const auto& fn = static_cast<const Function&>(*sol);
    return std::string_view(fn.name()) == "Dsolve";
}

}  // namespace

// ----- First-order ----------------------------------------------------------

TEST_CASE("matlab::dsolve first-order linear: y' + y = 0",
          "[15m][matlab][dsolve]") {
    auto [x, y, yp] = matlab::syms("x", "y", "yp");
    auto sol = matlab::dsolve(yp + y, y, yp, x);
    REQUIRE(sol);
    REQUIRE_FALSE(is_unevaluated_marker(sol));
}

TEST_CASE("matlab::dsolve_ivp first-order applies condition",
          "[15m][matlab][dsolve][ivp]") {
    auto [x, y, yp] = matlab::syms("x", "y", "yp");
    // y' + y = 0, y(0) = 1 → y = exp(-x)
    std::vector<std::pair<Expr, Expr>> conds = {{S::Zero(), integer(1)}};
    auto sol = matlab::dsolve_ivp(yp + y, y, yp, x, conds);
    REQUIRE(sol);
    REQUIRE_FALSE(is_unevaluated_marker(sol));
}

// ----- Second-order: constant-coefficient -----------------------------------

TEST_CASE("matlab::dsolve second-order const-coef: y'' + y = 0",
          "[15m][matlab][dsolve]") {
    auto [x, y, yp, ypp] = matlab::syms("x", "y", "yp", "ypp");
    auto sol = matlab::dsolve(ypp + y, y, yp, ypp, x);
    REQUIRE(sol);
    REQUIRE_FALSE(is_unevaluated_marker(sol));
    // y = C0 cos x + C1 sin x — Add of two terms.
    REQUIRE(sol->type_id() == TypeId::Add);
}

TEST_CASE("matlab::dsolve second-order const-coef: y'' - 3y' + 2y = 0",
          "[15m][matlab][dsolve]") {
    auto [x, y, yp, ypp] = matlab::syms("x", "y", "yp", "ypp");
    auto sol = matlab::dsolve(
        ypp - integer(3) * yp + integer(2) * y, y, yp, ypp, x);
    REQUIRE(sol);
    REQUIRE_FALSE(is_unevaluated_marker(sol));
    // Verify residual: y'' - 3y' + 2y = 0.
    auto yp_val = matlab::diff(sol, x);
    auto ypp_val = matlab::diff(yp_val, x);
    auto residual = matlab::simplify(
        ypp_val - integer(3) * yp_val + integer(2) * sol);
    REQUIRE(residual == S::Zero());
}

// ----- Second-order: Cauchy-Euler -------------------------------------------

TEST_CASE("matlab::dsolve second-order Cauchy-Euler: x²y'' + xy' - y = 0",
          "[15m][matlab][dsolve][cauchy_euler]") {
    auto [x, y, yp, ypp] = matlab::syms("x", "y", "yp", "ypp");
    // Coefficient pattern: a₂(x) = x², a₁(x) = x, a₀(x) = -1.
    auto eq = pow(x, integer(2)) * ypp + x * yp - y;
    auto sol = matlab::dsolve(eq, y, yp, ypp, x);
    REQUIRE(sol);
    REQUIRE_FALSE(is_unevaluated_marker(sol));
}

// ----- Second-order: nonhomogeneous → variation of parameters --------------

TEST_CASE("matlab::dsolve nonhomogeneous: y'' - 3y' + 2y = x routes through varparams",
          "[15m][matlab][dsolve][varparams][oracle]") {
    auto& oracle = sympp::testing::Oracle::instance();
    auto [x, y, yp, ypp] = matlab::syms("x", "y", "yp", "ypp");
    // y'' - 3y' + 2y - x = 0  →  RHS g(x) = x.
    auto eq = ypp - integer(3) * yp + integer(2) * y - x;
    auto sol = matlab::dsolve(eq, y, yp, ypp, x);
    REQUIRE(sol);
    // The unevaluated marker would have type_id Function with name "Dsolve".
    bool is_marker = false;
    if (sol->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*sol);
        is_marker = (std::string_view(fn.name()) == "Dsolve");
    }
    REQUIRE_FALSE(is_marker);
    // Verify residual at a numeric point.
    auto yp_val = matlab::diff(sol, x);
    auto ypp_val = matlab::diff(yp_val, x);
    auto residual = matlab::simplify(
        ypp_val - integer(3) * yp_val + integer(2) * sol - x);
    auto resp = oracle.send(
        {{"op", "evalf_is_zero"},
         {"expr", matlab::simplify(matlab::subs(residual,
                       x, rational(7, 5)))->str()},
         {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

// ----- Second-order: nonlinear → marker -------------------------------------

TEST_CASE("matlab::dsolve nonlinear ODE returns unevaluated Dsolve marker",
          "[15m][matlab][dsolve]") {
    auto [x, y, yp, ypp] = matlab::syms("x", "y", "yp", "ypp");
    // y·y'' = 1 — the y² term makes it nonlinear in y.
    auto eq = y * ypp - integer(1);
    auto sol = matlab::dsolve(eq, y, yp, ypp, x);
    REQUIRE(sol);
    REQUIRE(is_unevaluated_marker(sol));
}

// ----- dsolve_constant_coeff re-export --------------------------------------

TEST_CASE("matlab::dsolve_constant_coeff direct coefficient list",
          "[15m][matlab][dsolve]") {
    auto x = symbol("x");
    // y'' - 3y' + 2y = 0 → coeffs [2, -3, 1]
    auto sol = matlab::dsolve_constant_coeff(
        {integer(2), integer(-3), integer(1)}, x);
    REQUIRE(sol);
}

// ----- Linear ODE system ----------------------------------------------------

TEST_CASE("matlab::dsolve(A, x) linear system",
          "[15m][matlab][dsolve][system]") {
    auto x = symbol("x");
    Matrix A = {{integer(2), integer(0)}, {integer(0), integer(3)}};
    auto sol = matlab::dsolve(A, x);
    REQUIRE(sol.rows() == 2);
    REQUIRE(sol.cols() == 1);
}
