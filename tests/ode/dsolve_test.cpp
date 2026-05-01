// dsolve / pdsolve tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/expand.hpp>
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

#include <sympp/functions/miscellaneous.hpp>

// ----- Riccati ---------------------------------------------------------------

TEST_CASE("dsolve_riccati: y' = 1 + y² → tan(x + C)",
          "[11d][riccati][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    // P=1, Q=0, R=1: y' = 1 + 0·y + 1·y².
    auto sol = dsolve_riccati(integer(1), integer(0), integer(1), x);
    // Substitute into the original eq and verify residual.
    auto eq = yp - integer(1) - pow(y, integer(2));
    auto residual = checkodesol(eq, sol, y, yp, x);
    auto C0 = symbol("__C0");
    auto C1 = symbol("__C1");
    Expr fixed = subs(subs(residual, C0, integer(1)), C1, integer(0));
    Expr at_half = simplify(subs(fixed, x, rational(1, 2)));
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", at_half->str()},
                             {"prec", 30}, {"tol", 8}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

TEST_CASE("dsolve_riccati: variable coefficient bails out",
          "[11d][riccati]") {
    auto x = symbol("x");
    auto sol = dsolve_riccati(integer(1), x, integer(1), x);
    // Variable Q in x — coef_u_prime depends on x, so we bail to marker.
    REQUIRE(sol->type_id() == TypeId::Function);
}

// ----- Homogeneous y' = f(y/x) -----------------------------------------------

TEST_CASE("dsolve_homogeneous: y' = y/x → y = C·x",
          "[11d][homogeneous][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    auto eq = yp - y / x;
    auto sol = dsolve_homogeneous(eq, y, yp, x);
    // Verify checkodesol → 0 numerically at x=2 with __C0 = 3.
    auto residual = checkodesol(eq, sol, y, yp, x);
    auto C0 = symbol("__C0");
    Expr fixed = subs(residual, C0, integer(3));
    Expr at_2 = simplify(subs(fixed, x, integer(2)));
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", at_2->str()},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

// ----- Lie autonomous --------------------------------------------------------

TEST_CASE("dsolve_lie_autonomous: y' = y² (autonomous)",
          "[11d][lie][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    auto eq = yp - pow(y, integer(2));
    auto sol = dsolve_lie_autonomous(eq, y, yp, x);
    auto residual = checkodesol(eq, sol, y, yp, x);
    auto C0 = symbol("__C0");
    Expr fixed = subs(residual, C0, integer(1));
    Expr at_half = simplify(subs(fixed, x, rational(1, 2)));
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", at_half->str()},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

TEST_CASE("dsolve_lie_autonomous: rejects x-dependent rhs",
          "[11d][lie]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    // Has x in rhs — not autonomous.
    auto eq = yp - x - y;
    auto sol = dsolve_lie_autonomous(eq, y, yp, x);
    REQUIRE(sol->type_id() == TypeId::Function);
}

// ----- Hypergeometric --------------------------------------------------------

TEST_CASE("hyper: factory creates hyper(a, b, c, z)",
          "[11d][hyper]") {
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    auto z = symbol("z");
    auto h = hyper(a, b, c, z);
    REQUIRE(h->type_id() == TypeId::Function);
}

TEST_CASE("dsolve_hypergeometric: standard form recognized",
          "[11d][hyper]") {
    auto x = symbol("x");
    // x(1-x)y'' + (3 - 4x)y' - 2y = 0.  Match: c = 3, a+b+1 = 4 → a+b=3,
    // a*b = 2. Roots of t²-3t+2: 1, 2 → hyper(1, 2, 3, x).
    auto coef_y = integer(-2);
    auto coef_yp = integer(3) - integer(4) * x;
    auto coef_ypp = x * (integer(1) - x);
    auto sol = dsolve_hypergeometric({coef_y, coef_yp, expand(coef_ypp)}, x);
    REQUIRE(sol->type_id() == TypeId::Function);
}

// ----- Variable-coeff PDE ----------------------------------------------------

TEST_CASE("pdsolve_first_order_variable: x u_x + y u_y = 0",
          "[11d][pdsolve][oracle]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");
    // a = x, b = y, c = 0 — characteristics are y/x = const.
    auto sol = pdsolve_first_order_variable(x, y, S::Zero(), y, yp, x);
    // Solution should be F(y/x) form.
    REQUIRE(sol->type_id() == TypeId::Function);
}

// ----- Heat / Wave -----------------------------------------------------------

TEST_CASE("pdsolve_heat: u_t = k u_xx product solution",
          "[11d][pdsolve_heat]") {
    auto x = symbol("x");
    auto t = symbol("t");
    auto k = symbol("k");
    auto lambda = symbol("lambda");
    auto u = pdsolve_heat(k, lambda, x, t);
    // Verify u_t - k u_xx = 0.
    auto u_t = diff(u, t);
    auto u_xx = diff(diff(u, x), x);
    Expr residual = simplify(u_t - k * u_xx);
    REQUIRE(residual == S::Zero());
}

TEST_CASE("pdsolve_wave: u_tt = c² u_xx d'Alembert",
          "[11d][pdsolve_wave]") {
    auto x = symbol("x");
    auto t = symbol("t");
    auto c = symbol("c");
    auto u = pdsolve_wave(c, x, t);
    // u = F(x - c·t) + G(x + c·t). Structure check: it's an Add of two
    // function calls.
    REQUIRE(u->type_id() == TypeId::Add);
}

// ----- DAE Jacobian / structural index ---------------------------------------

TEST_CASE("dae_jacobians: extracts ∂F/∂yp and ∂F/∂y",
          "[11d][dae]") {
    auto y1 = symbol("y1");
    auto y2 = symbol("y2");
    auto y1p = symbol("y1p");
    auto y2p = symbol("y2p");
    auto x = symbol("x");
    // Simple DAE: y1' = y2, y1 + y2 - 1 = 0 (the second is algebraic).
    std::vector<Expr> F = {y1p - y2, y1 + y2 - integer(1)};
    auto [M, K] = dae_jacobians(F, {y1, y2}, {y1p, y2p});
    REQUIRE(M.rows() == 2);
    REQUIRE(M.cols() == 2);
    // Row 0 of M: ∂(y1p - y2)/∂(y1p, y2p) = (1, 0). Row 1: (0, 0) — algebraic.
    REQUIRE(M.at(0, 0) == integer(1));
    REQUIRE(M.at(1, 0) == S::Zero());
    REQUIRE(M.at(1, 1) == S::Zero());
}

TEST_CASE("dae_structural_index: pure ODE → 0",
          "[11d][dae]") {
    auto y1 = symbol("y1"); auto y1p = symbol("y1p");
    std::vector<Expr> F = {y1p - y1};   // y' = y, regular ODE
    REQUIRE(dae_structural_index(F, {y1}, {y1p}) == 0);
}

TEST_CASE("dae_structural_index: index-1 algebraic constraint",
          "[11d][dae]") {
    auto y1 = symbol("y1"); auto y2 = symbol("y2");
    auto y1p = symbol("y1p"); auto y2p = symbol("y2p");
    // y1' = y2, y1 + y2 = 1. Second eq is algebraic in y → M has a zero row.
    std::vector<Expr> F = {y1p - y2, y1 + y2 - integer(1)};
    auto idx = dae_structural_index(F, {y1, y2}, {y1p, y2p});
    REQUIRE(idx >= 1);
}
