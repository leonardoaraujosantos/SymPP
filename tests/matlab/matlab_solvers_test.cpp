// MATLAB facade — solvers tests.
//
// Validates linsolve / multi-eq solve / nsolve / vpasolve wrappers.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/matlab/matlab.hpp>
#include <sympp/matrices/matrix.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- linsolve --------------------------------------------------------------

TEST_CASE("matlab::linsolve 2x2 system",
          "[15m][matlab][linsolve]") {
    Matrix A = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    Matrix b = {{integer(5)}, {integer(11)}};
    Matrix x = matlab::linsolve(A, b);
    REQUIRE(x.rows() == 2);
    REQUIRE(x.cols() == 1);
    REQUIRE(x.at(0, 0) == integer(1));
    REQUIRE(x.at(1, 0) == integer(2));
}

TEST_CASE("matlab::linsolve singular A throws",
          "[15m][matlab][linsolve]") {
    Matrix A = {{integer(1), integer(2)}, {integer(2), integer(4)}};
    Matrix b = {{integer(1)}, {integer(2)}};
    REQUIRE_THROWS(matlab::linsolve(A, b));
}

// ----- multi-equation solve --------------------------------------------------

TEST_CASE("matlab::solve 2x2 polynomial system via nonlinsolve",
          "[15m][matlab][solve]") {
    auto x = symbol("x");
    auto y = symbol("y");
    // x² + y² = 1, x = y → ±(√½, √½)
    auto eq1 = pow(x, integer(2)) + pow(y, integer(2)) - integer(1);
    auto eq2 = x - y;
    auto sols = matlab::solve({eq1, eq2}, {x, y});
    REQUIRE(sols.size() >= 1);
    for (const auto& s : sols) {
        REQUIRE(s.size() == 2);
    }
}

TEST_CASE("matlab::solve 3-var system routes to Groebner path",
          "[15m][matlab][solve]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    // Trivial triangular system: x = 1, y = 2, z = 3
    std::vector<Expr> eqs = {x - integer(1), y - integer(2), z - integer(3)};
    auto sols = matlab::solve(eqs, {x, y, z});
    REQUIRE(sols.size() == 1);
    REQUIRE(sols[0].size() == 3);
    REQUIRE(sols[0][0] == integer(1));
    REQUIRE(sols[0][1] == integer(2));
    REQUIRE(sols[0][2] == integer(3));
}

// ----- nsolve / vpasolve -----------------------------------------------------

TEST_CASE("matlab::nsolve x³ - 2 ≈ 1.2599 at default 15 dps",
          "[15m][matlab][nsolve][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto root = matlab::nsolve(pow(x, integer(3)) - integer(2),
                                  x, integer(1));
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", root->str() + " - 1.2599210498948732"},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

TEST_CASE("matlab::vpasolve cos(x) - x ≈ 0.7390851 at default 32 dps",
          "[15m][matlab][vpasolve][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto root = matlab::vpasolve(cos(x) - x, x, integer(1));
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", root->str() + " - 0.7390851332151607"},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}
