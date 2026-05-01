// dsolve / pdsolve tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/ode/dsolve.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Separable ODE --------------------------------------------------------

TEST_CASE("dsolve_separable: y' = y → y = C·e^x",
          "[11][dsolve][separable][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    // eq: yp - y = 0
    auto eq = yp - y;
    auto sol = dsolve_separable(eq, y, yp, x);
    // Verify by checkodesol — substitute y → sol, yp → diff(sol, x).
    auto residual = checkodesol(eq, sol, y, yp, x);
    REQUIRE(oracle.equivalent(residual->str(), "0"));
}

TEST_CASE("dsolve_separable: y' = x·y → y = C·exp(x²/2)",
          "[11][dsolve][separable][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    auto eq = yp - x * y;
    auto sol = dsolve_separable(eq, y, yp, x);
    auto residual = checkodesol(eq, sol, y, yp, x);
    REQUIRE(oracle.equivalent(residual->str(), "0"));
}

TEST_CASE("dsolve_separable: y' = x² → y = x³/3 + C",
          "[11][dsolve][separable][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    auto eq = yp - pow(x, integer(2));
    auto sol = dsolve_separable(eq, y, yp, x);
    auto residual = checkodesol(eq, sol, y, yp, x);
    REQUIRE(oracle.equivalent(residual->str(), "0"));
}

// ----- Linear first-order ----------------------------------------------------

TEST_CASE("dsolve_linear_first_order: y' + y = e^x",
          "[11][dsolve][linear1][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    auto eq = yp + y - exp(x);
    auto sol = dsolve_linear_first_order(eq, y, yp, x);
    auto residual = checkodesol(eq, sol, y, yp, x);
    REQUIRE(oracle.equivalent(residual->str(), "0"));
}

TEST_CASE("dsolve_linear_first_order: y' + 2y = 0 → y = C·e^(-2x)",
          "[11][dsolve][linear1][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    auto eq = yp + integer(2) * y;
    auto sol = dsolve_linear_first_order(eq, y, yp, x);
    auto residual = checkodesol(eq, sol, y, yp, x);
    REQUIRE(oracle.equivalent(residual->str(), "0"));
}

// ----- Bernoulli -------------------------------------------------------------

TEST_CASE("dsolve_bernoulli: y' + y = y² (n=2)",
          "[11][dsolve][bernoulli][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    auto eq = yp + y - pow(y, integer(2));
    auto sol = dsolve_bernoulli(eq, y, yp, x);
    auto residual = checkodesol(eq, sol, y, yp, x);
    // Substitute __C0 = 1 and x = 1 to get a numeric residual; should be 0.
    auto C0 = symbol("__C0");
    Expr fixed = subs(residual, C0, integer(1));
    Expr at_1 = simplify(subs(fixed, x, integer(1)));
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", at_1->str()},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

// ----- Exact ODE -------------------------------------------------------------

TEST_CASE("dsolve_exact: 2xy + (x² + 3y²) y' = 0",
          "[11][dsolve][exact][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    // M = 2xy, N = x² + 3y² — exact since ∂M/∂y = 2x = ∂N/∂x.
    auto M = integer(2) * x * y;
    auto N = pow(x, integer(2)) + integer(3) * pow(y, integer(2));
    auto sol = dsolve_exact(M, N, y, x);
    // Solution is x²·y + y³ - C. Check by differentiating implicitly:
    //   ∂ψ/∂x + ∂ψ/∂y · y' = 0  → 2xy + (x² + 3y²) y' = 0.
    auto psi_x = diff(sol, x);
    auto psi_y = diff(sol, y);
    Expr lhs = psi_x + psi_y * (mul(S::NegativeOne(), M) / N);
    REQUIRE(oracle.equivalent(simplify(lhs)->str(), "0"));
}

// ----- Constant-coefficient linear ODE --------------------------------------

TEST_CASE("dsolve_constant_coeff: y'' + y = 0 → C₀ cos x + C₁ sin x form",
          "[11][dsolve][constcoeff]") {
    auto x = symbol("x");
    // y'' + y = 0 → coeffs [1, 0, 1] (a₀=1, a₁=0, a₂=1).
    auto sol = dsolve_constant_coeff({integer(1), integer(0), integer(1)}, x);
    // Char poly r² + 1 = 0 → r = ±I.
    // General solution: C₀ exp(I x) + C₁ exp(-I x). The shape is fine;
    // we check it's an Add with two terms.
    REQUIRE(sol->type_id() == TypeId::Add);
    REQUIRE(sol->args().size() == 2);
}

TEST_CASE("dsolve_constant_coeff: y'' - 3y' + 2y = 0 → C₀ e^x + C₁ e^(2x)",
          "[11][dsolve][constcoeff][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // coeffs [2, -3, 1] for 2 y - 3 y' + y'' = 0.
    auto sol = dsolve_constant_coeff(
        {integer(2), integer(-3), integer(1)}, x);
    // Verify by differentiating: y'' - 3y' + 2y should be 0.
    auto yp = diff(sol, x);
    auto ypp = diff(yp, x);
    Expr residual = simplify(ypp - integer(3) * yp + integer(2) * sol);
    REQUIRE(oracle.equivalent(residual->str(), "0"));
}

TEST_CASE("dsolve_constant_coeff: repeated root y'' - 2y' + y = 0",
          "[11][dsolve][constcoeff][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // (r-1)² = r² - 2r + 1, coeffs [1, -2, 1]
    auto sol = dsolve_constant_coeff(
        {integer(1), integer(-2), integer(1)}, x);
    auto yp = diff(sol, x);
    auto ypp = diff(yp, x);
    Expr residual = simplify(ypp - integer(2) * yp + sol);
    REQUIRE(oracle.equivalent(residual->str(), "0"));
}

// ----- Cauchy-Euler ----------------------------------------------------------

TEST_CASE("dsolve_cauchy_euler: x²y'' - x y' + y = 0",
          "[11][dsolve][cauchy_euler][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // a₀=1 (y term), a₁=-1 (x y' term), a₂=1 (x² y'' term)
    auto sol = dsolve_cauchy_euler({integer(1), integer(-1), integer(1)}, x);
    auto yp = diff(sol, x);
    auto ypp = diff(yp, x);
    Expr residual = simplify(pow(x, integer(2)) * ypp - x * yp + sol);
    // Should be zero. May include log(x) terms from repeated roots.
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", simplify(subs(residual,
                                  x, integer(2)))->str()},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

// ----- ODE systems -----------------------------------------------------------

TEST_CASE("dsolve_system: diagonal A → independent exponentials",
          "[11][dsolve][system]") {
    auto x = symbol("x");
    Matrix A = {{integer(2), integer(0)}, {integer(0), integer(3)}};
    auto sol = dsolve_system(A, x);
    REQUIRE(sol.rows() == 2);
    REQUIRE(sol.cols() == 1);
    // Solution should be C₀ e^(2x) and C₁ e^(3x) in some order.
}

TEST_CASE("dsolve_system: 2x2 satisfies y' = A·y",
          "[11][dsolve][system][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    Matrix A = {{integer(2), integer(0)}, {integer(0), integer(3)}};
    auto sol = dsolve_system(A, x);
    // Compute y' and A·y, verify equal.
    Matrix yp(2, 1);
    yp.set(0, 0, diff(sol.at(0, 0), x));
    yp.set(1, 0, diff(sol.at(1, 0), x));
    auto Ay = A * sol;
    for (std::size_t i = 0; i < 2; ++i) {
        Expr d = simplify(yp.at(i, 0) - Ay.at(i, 0));
        REQUIRE(oracle.equivalent(d->str(), "0"));
    }
}

// ----- IVP -------------------------------------------------------------------

TEST_CASE("apply_ivp: y' = y, y(0) = 5 → y = 5·e^x",
          "[11][dsolve][ivp][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    auto sol = dsolve_separable(yp - y, y, yp, x);
    auto fixed = apply_ivp(sol, x, {{S::Zero(), integer(5)}});
    // At x=0, fixed should equal 5.
    Expr at_0 = simplify(subs(fixed, x, S::Zero()));
    REQUIRE(oracle.equivalent(at_0->str(), "5"));
}

// ----- checkodesol -----------------------------------------------------------

TEST_CASE("checkodesol: rejects wrong solution", "[11][dsolve][check]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    auto eq = yp - y;
    // Wrong: y = x is not a solution to y' = y.
    auto wrong = x;
    auto residual = checkodesol(eq, wrong, y, yp, x);
    // Residual should be 1 - x ≠ 0.
    REQUIRE(!(residual == S::Zero()));
}

// ----- PDE: first-order linear -----------------------------------------------

TEST_CASE("pdsolve_first_order_linear: 2 u_x + 3 u_y = 0",
          "[11][pdsolve]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto sol = pdsolve_first_order_linear(integer(2), integer(3), S::Zero(),
                                            x, y);
    // u(x, y) = F(3x - 2y). Check it has the F-call structure.
    REQUIRE(sol->type_id() == TypeId::Function);
}

TEST_CASE("pdsolve_first_order_linear: 1 u_x + 1 u_y = 5",
          "[11][pdsolve]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto sol = pdsolve_first_order_linear(integer(1), integer(1), integer(5),
                                            x, y);
    // u(x, y) = 5x + F(x - y).
    REQUIRE(sol->type_id() == TypeId::Add);
}
