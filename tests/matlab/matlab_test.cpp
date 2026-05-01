// MATLAB facade tests.
//
// Validates that the sympp::matlab namespace exposes MATLAB-named
// wrappers that do the same thing as the underlying SymPP API.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/matlab/matlab.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- syms / sym -----------------------------------------------------------

TEST_CASE("matlab::syms single", "[15m][matlab][syms]") {
    auto x = matlab::syms("x");
    REQUIRE(x == symbol("x"));
}

TEST_CASE("matlab::syms multi via structured binding",
          "[15m][matlab][syms]") {
    auto [x, y, z] = matlab::syms("x", "y", "z");
    REQUIRE(x == symbol("x"));
    REQUIRE(y == symbol("y"));
    REQUIRE(z == symbol("z"));
}

TEST_CASE("matlab::sym from string and integer",
          "[15m][matlab][sym]") {
    REQUIRE(matlab::sym("alpha") == symbol("alpha"));
    REQUIRE(matlab::sym(42L) == integer(42));
}

// ----- vpa ------------------------------------------------------------------

TEST_CASE("matlab::vpa evalfs at requested precision",
          "[15m][matlab][vpa][oracle]") {
    auto& oracle = Oracle::instance();
    auto v = matlab::vpa(S::Pi(), 20);
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", v->str() + " - 3.14159265358979323846"},
                             {"prec", 30}, {"tol", 18}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

// ----- diff -----------------------------------------------------------------

TEST_CASE("matlab::diff first derivative",
          "[15m][matlab][diff][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = matlab::syms("x");
    auto df = matlab::diff(sin(x) * exp(x), x);
    REQUIRE(oracle.equivalent(df->str(), "exp(x)*sin(x) + exp(x)*cos(x)"));
}

TEST_CASE("matlab::diff n-th derivative",
          "[15m][matlab][diff][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = matlab::syms("x");
    // d²/dx² (x⁴) = 12x²
    auto d2 = matlab::diff(pow(x, integer(4)), x, 2);
    REQUIRE(oracle.equivalent(d2->str(), "12*x**2"));
}

// ----- Int ------------------------------------------------------------------

TEST_CASE("matlab::Int indefinite", "[15m][matlab][int][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = matlab::syms("x");
    auto F = matlab::Int(integer(2) * x, x);
    REQUIRE(oracle.equivalent(F->str(), "x**2"));
}

TEST_CASE("matlab::Int definite", "[15m][matlab][int]") {
    auto x = matlab::syms("x");
    auto v = matlab::Int(sin(x), x, S::Zero(), S::Pi());
    REQUIRE(v == integer(2));
}

// ----- Algebra wrappers -----------------------------------------------------

TEST_CASE("matlab::simplify chains the orchestrator",
          "[15m][matlab][simplify][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = matlab::syms("x");
    auto e = pow(sin(x), integer(2)) + pow(cos(x), integer(2));
    REQUIRE(oracle.equivalent(matlab::simplify(e)->str(), "1"));
}

TEST_CASE("matlab::factor", "[15m][matlab][factor][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = matlab::syms("x");
    auto f = matlab::factor(pow(x, integer(2)) - integer(1), x);
    REQUIRE(oracle.equivalent(f->str(), "(x - 1)*(x + 1)"));
}

TEST_CASE("matlab::expand", "[15m][matlab][expand][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = matlab::syms("x");
    auto e = matlab::expand(pow(x + integer(1), integer(2)));
    REQUIRE(oracle.equivalent(e->str(), "x**2 + 2*x + 1"));
}

TEST_CASE("matlab::subs", "[15m][matlab][subs]") {
    auto x = matlab::syms("x");
    auto s = matlab::subs(pow(x, integer(2)) + x, x, integer(3));
    REQUIRE(s == integer(12));
}

// ----- Solvers --------------------------------------------------------------

TEST_CASE("matlab::solve quadratic", "[15m][matlab][solve]") {
    auto x = matlab::syms("x");
    auto roots = matlab::solve(pow(x, integer(2)) - integer(4), x);
    REQUIRE(roots.size() == 2);
}

TEST_CASE("matlab::dsolve linear first-order",
          "[15m][matlab][dsolve]") {
    auto [x, y, yp] = matlab::syms("x", "y", "yp");
    auto sol = matlab::dsolve(yp + y - exp(x), y, yp, x);
    REQUIRE(sol);
    // Verify the solver returned something that isn't the Dsolve(...)
    // unevaluated marker. Wrap the comparisons in parens because
    // Catch2 disallows chained comparisons inside REQUIRE.
    bool is_marker = false;
    if (sol->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*sol);
        is_marker = (std::string_view(fn.name()) == "Dsolve");
    }
    REQUIRE(!is_marker);
}

// ----- Series -> taylor ------------------------------------------------------

TEST_CASE("matlab::taylor exp(x) about 0",
          "[15m][matlab][taylor][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = matlab::syms("x");
    auto s = matlab::taylor(exp(x), x, S::Zero(), 5);
    REQUIRE(oracle.equivalent(s->str(),
                              "1 + x + x**2/2 + x**3/6 + x**4/24"));
}

// ----- Transforms -----------------------------------------------------------

TEST_CASE("matlab::laplace and ilaplace round-trip on exp(a*t)",
          "[15m][matlab][laplace]") {
    auto [a, t, s] = matlab::syms("a", "t", "s");
    auto F = matlab::laplace(exp(a * t), t, s);
    auto f = matlab::ilaplace(F, s, t);
    REQUIRE(f == exp(a * t));
}

TEST_CASE("matlab::ztrans on a^n", "[15m][matlab][ztrans][oracle]") {
    auto& oracle = Oracle::instance();
    auto [a, n, z] = matlab::syms("a", "n", "z");
    auto F = matlab::ztrans(pow(a, n), n, z);
    REQUIRE(oracle.equivalent(F->str(), "z/(z - a)"));
}

// ----- limit ----------------------------------------------------------------

TEST_CASE("matlab::limit sin(x)/x at 0 → 1",
          "[15m][matlab][limit]") {
    auto x = matlab::syms("x");
    REQUIRE(matlab::limit(sin(x) / x, x, S::Zero()) == integer(1));
}

// ----- Code generation ------------------------------------------------------

TEST_CASE("matlab::ccode", "[15m][matlab][ccode]") {
    auto x = matlab::syms("x");
    REQUIRE(matlab::ccode(pow(x, integer(2)) + sin(x))
             == "pow(x, 2) + sin(x)");
}

TEST_CASE("matlab::latex", "[15m][matlab][latex]") {
    auto x = matlab::syms("x");
    REQUIRE(matlab::latex(rational(1, 2)) == "\\frac{1}{2}");
}

TEST_CASE("matlab::fortran", "[15m][matlab][fortran]") {
    auto x = matlab::syms("x");
    auto out = matlab::fortran(pow(x, integer(2)));
    REQUIRE(out == "x**2.0d0");
}

TEST_CASE("matlab::matlabFunction emits valid Octave",
          "[15m][matlab][matlabFunction]") {
    auto x = matlab::syms("x");
    auto src = matlab::matlabFunction("f",
                                        pow(x, integer(2)) + integer(3) * x,
                                        {x});
    REQUIRE(src.find("function y = f(x)") != std::string::npos);
    REQUIRE(src.find(".^") != std::string::npos);
}

// ----- Worked example: porting a tiny MATLAB script -------------------------

TEST_CASE("matlab: full toolbox-style script — symbolic kinematics",
          "[15m][matlab][example][oracle]") {
    // Mirror a small MATLAB workflow using only the matlab:: namespace
    // (no `using namespace sympp` — that would shadow-collide with
    // matlab:: names of identical signatures).
    auto& oracle = Oracle::instance();

    // syms x t
    auto [x, t] = matlab::syms("x", "t");
    // pos = sin(t) * exp(t)
    auto pos = sin(t) * exp(t);
    // velocity = diff(x, t)
    auto velocity = matlab::diff(pos, t);
    // acceleration = diff(velocity, t)
    auto acceleration = matlab::diff(velocity, t);
    // simplify the expressions
    auto v_simp = matlab::simplify(velocity);
    auto a_simp = matlab::simplify(acceleration);

    // sanity: position evaluated at t = 0 is 0
    REQUIRE(matlab::subs(pos, t, S::Zero()) == S::Zero());
    // velocity at t = 0 = sin(0)·1 + cos(0)·1 = 1
    REQUIRE(oracle.equivalent(
        matlab::subs(v_simp, t, S::Zero())->str(), "1"));
    // acceleration at t = 0 = 2 (from d²/dt² (sin(t)·exp(t)) at 0)
    REQUIRE(oracle.equivalent(
        matlab::subs(a_simp, t, S::Zero())->str(), "2"));

    // codegen: emit C function for velocity
    auto src = matlab::ccode(v_simp);
    REQUIRE(src.find("exp(t)") != std::string::npos);
}
