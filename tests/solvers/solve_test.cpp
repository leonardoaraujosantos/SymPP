// solve / linsolve tests.
//
// References:
//   sympy/solvers/tests/test_solveset.py
//   MATLAB Symbolic Math Toolbox §3 — solve(eqn, var), linsolve

#include <stdexcept>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/solvers/solve.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- solve (univariate polynomial) ----------------------------------------

TEST_CASE("solve: linear equation", "[10][solve]") {
    auto x = symbol("x");
    // 3x - 7 = 0 -> x = 7/3
    auto roots = solve(integer(3) * x - integer(7), x);
    REQUIRE(roots.size() == 1);
    REQUIRE(roots[0] == rational(7, 3));
}

TEST_CASE("solve: lhs == rhs form", "[10][solve]") {
    auto x = symbol("x");
    // 2x + 1 == 5 -> x = 2
    auto roots = solve(integer(2) * x + integer(1), integer(5), x);
    REQUIRE(roots.size() == 1);
    REQUIRE(roots[0] == integer(2));
}

TEST_CASE("solve: quadratic with rational roots", "[10][solve]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^2 - 5x + 6 = 0 -> x ∈ {2, 3}
    auto roots = solve(pow(x, integer(2)) - integer(5) * x + integer(6), x);
    REQUIRE(roots.size() == 2);
    REQUIRE(oracle.equivalent(roots[0]->str(), "3"));
    REQUIRE(oracle.equivalent(roots[1]->str(), "2"));
}

TEST_CASE("solve: quadratic with irrational roots", "[10][solve][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^2 - 2 = 0 -> ±sqrt(2)
    auto roots = solve(pow(x, integer(2)) - integer(2), x);
    REQUIRE(roots.size() == 2);
    REQUIRE(oracle.equivalent(roots[0]->str(), "sqrt(2)"));
    REQUIRE(oracle.equivalent(roots[1]->str(), "-sqrt(2)"));
}

TEST_CASE("solve: quadratic with parametric coefficient",
          "[10][solve][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    // a*x^2 + b*x + c = 0
    auto roots = solve(a * pow(x, integer(2)) + b * x + c, x);
    REQUIRE(roots.size() == 2);
    INFO("root[0]: " << roots[0]->str());
    INFO("root[1]: " << roots[1]->str());
    REQUIRE(oracle.equivalent(roots[0]->str(),
                              "(-b + sqrt(b**2 - 4*a*c)) / (2*a)"));
    REQUIRE(oracle.equivalent(roots[1]->str(),
                              "(-b - sqrt(b**2 - 4*a*c)) / (2*a)"));
}

TEST_CASE("solve: degree 3+ returns empty (deferred)", "[10][solve]") {
    auto x = symbol("x");
    auto roots = solve(pow(x, integer(3)) - integer(8), x);
    REQUIRE(roots.empty());  // Phase 10 minimal: no Cardano yet
}

// ----- linsolve --------------------------------------------------------------

TEST_CASE("linsolve: 2x2 system", "[10][linsolve]") {
    // x + 2y = 5
    // 3x + 4y = 11
    // -> x = 1, y = 2
    Matrix A = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    Matrix b = {{integer(5)}, {integer(11)}};
    Matrix x = linsolve(A, b);
    REQUIRE(x.rows() == 2);
    REQUIRE(x.cols() == 1);
    REQUIRE(x.at(0, 0) == integer(1));
    REQUIRE(x.at(1, 0) == integer(2));
}

TEST_CASE("linsolve: singular A throws", "[10][linsolve]") {
    Matrix A = {{integer(1), integer(2)}, {integer(2), integer(4)}};
    Matrix b = {{integer(1)}, {integer(2)}};
    REQUIRE_THROWS(linsolve(A, b));
}

TEST_CASE("linsolve: identity returns b unchanged", "[10][linsolve]") {
    Matrix A = Matrix::identity(3);
    Matrix b = {{integer(7)}, {integer(11)}, {integer(13)}};
    Matrix x = linsolve(A, b);
    REQUIRE(x.equals(b));
}
