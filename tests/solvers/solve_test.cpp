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

TEST_CASE("solve: cubic via rational-root deflation", "[10][solve][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^3 - 8 = (x - 2)(x² + 2x + 4); rational root 2 + two from quadratic.
    auto roots = solve(pow(x, integer(3)) - integer(8), x);
    REQUIRE(roots.size() == 3);
    // Each root must satisfy x³ - 8 = 0.
    for (const auto& r : roots) {
        auto val = pow(r, integer(3)) - integer(8);
        REQUIRE(oracle.equivalent(val->str(), "0"));
    }
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

#include <sympp/core/boolean.hpp>
#include <sympp/sets/sets.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/trigonometric.hpp>

// ----- solveset --------------------------------------------------------------

TEST_CASE("solveset: x^2 - 4 = 0 → {2, -2}", "[10][solveset]") {
    auto x = symbol("x");
    auto s = solveset(pow(x, integer(2)) - integer(4), x);
    REQUIRE(s->kind() == SetKind::FiniteSet);
    auto fs = std::static_pointer_cast<const FiniteSet>(s);
    REQUIRE(fs->size() == 2);
}

TEST_CASE("solveset: no solutions → EmptySet",
          "[10][solveset]") {
    auto x = symbol("x");
    // Constant 5 = 0 has no solution.
    auto s = solveset(integer(5), x);
    REQUIRE(s->kind() == SetKind::Empty);
}

TEST_CASE("solveset: filters via domain",
          "[10][solveset]") {
    auto x = symbol("x");
    auto eq = pow(x, integer(2)) - integer(4);
    auto domain = interval(integer(0), integer(10));
    auto s = solveset(eq, x, domain);
    // Only x=2 falls in [0, 10]; x=-2 is filtered out.
    REQUIRE(s->kind() == SetKind::FiniteSet);
    auto fs = std::static_pointer_cast<const FiniteSet>(s);
    REQUIRE(fs->size() == 1);
    REQUIRE(fs->elements()[0] == integer(2));
}

// ----- nsolve ----------------------------------------------------------------

TEST_CASE("nsolve: x^3 - 2 ≈ 1.2599",
          "[10][nsolve][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto root = nsolve(pow(x, integer(3)) - integer(2), x, integer(1), 15);
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", root->str() + " - 1.2599210498948732"},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

TEST_CASE("nsolve: cos(x) - x ≈ 0.7390851332151607",
          "[10][nsolve][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto root = nsolve(cos(x) - x, x, integer(1), 15);
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", root->str() + " - 0.7390851332151607"},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

// ----- solve_univariate_inequality -------------------------------------------

TEST_CASE("inequality: x^2 - 4 < 0 → (-2, 2)",
          "[10][inequality]") {
    auto x = symbol("x");
    auto op = symbol("<");
    auto s = solve_univariate_inequality(
        pow(x, integer(2)) - integer(4), op, S::Zero(), x);
    // The solution should be a single open interval roughly (-2, 2).
    REQUIRE(s->kind() != SetKind::Empty);
    REQUIRE(s->contains(integer(0)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(3)) == std::optional<bool>{false});
}

TEST_CASE("inequality: x^2 - 4 > 0 outside (-2, 2)",
          "[10][inequality]") {
    auto x = symbol("x");
    auto op = symbol(">");
    auto s = solve_univariate_inequality(
        pow(x, integer(2)) - integer(4), op, S::Zero(), x);
    REQUIRE(s->contains(integer(3)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(0)) == std::optional<bool>{false});
}

// ----- rsolve ----------------------------------------------------------------

TEST_CASE("rsolve: y(n+1) - 2*y(n) = 0 → y = C * 2^n",
          "[10][rsolve][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    // Coeffs: c_0 = -2, c_1 = 1 (i.e., -2 y(n) + y(n+1) = 0)
    auto sol = rsolve({integer(-2), integer(1)}, n);
    // Should be C * 2^n for some symbolic C.
    REQUIRE(oracle.equivalent(sol->str(), "__C0 * 2**n"));
}

TEST_CASE("rsolve: Fibonacci recurrence y(n+2) = y(n+1) + y(n)",
          "[10][rsolve]") {
    auto n = symbol("n");
    // y(n+2) - y(n+1) - y(n) = 0  → coeffs [-1, -1, 1]
    auto sol = rsolve({integer(-1), integer(-1), integer(1)}, n);
    // Solution: C0 * φ^n + C1 * (1-φ)^n with φ = (1+√5)/2.
    REQUIRE(sol->type_id() == TypeId::Add);
}

// ----- nonlinsolve -----------------------------------------------------------

TEST_CASE("nonlinsolve: x^2 + y^2 = 1 ∧ x = y",
          "[10][nonlinsolve][oracle]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto eq1 = pow(x, integer(2)) + pow(y, integer(2)) - integer(1);
    auto eq2 = x - y;
    auto sols = nonlinsolve({eq1, eq2}, {x, y});
    // Two solutions: (sqrt(2)/2, sqrt(2)/2) and (-sqrt(2)/2, -sqrt(2)/2).
    REQUIRE(sols.size() == 2);
}

// ----- linear_diophantine ----------------------------------------------------

TEST_CASE("diophantine: 3x + 5y = 1 has solution",
          "[10][diophantine]") {
    auto sol = linear_diophantine(integer(3), integer(5), integer(1));
    REQUIRE(sol.has_value());
    auto [x, y] = *sol;
    // Verify 3*x + 5*y = 1.
    auto check = integer(3) * x + integer(5) * y;
    REQUIRE(check == integer(1));
}

TEST_CASE("diophantine: 6x + 9y = 5 has no solution (gcd = 3 ∤ 5)",
          "[10][diophantine]") {
    auto sol = linear_diophantine(integer(6), integer(9), integer(5));
    REQUIRE(!sol.has_value());
}

TEST_CASE("diophantine: 4x + 6y = 10 (gcd = 2 ∣ 10)",
          "[10][diophantine]") {
    auto sol = linear_diophantine(integer(4), integer(6), integer(10));
    REQUIRE(sol.has_value());
    auto [x, y] = *sol;
    auto check = integer(4) * x + integer(6) * y;
    REQUIRE(check == integer(10));
}

// ----- Trig solveset → ImageSet ----------------------------------------------

TEST_CASE("solveset: sin(x) = 0 → ImageSet of n*pi",
          "[10d][solveset][trig]") {
    auto x = symbol("x");
    auto s = solveset(sin(x), x);
    REQUIRE(s->kind() == SetKind::ImageSet);
}

TEST_CASE("solveset: cos(x) = 0 → ImageSet of (2n+1)*pi/2",
          "[10d][solveset][trig]") {
    auto x = symbol("x");
    auto s = solveset(cos(x), x);
    REQUIRE(s->kind() == SetKind::ImageSet);
}

TEST_CASE("solveset: sin(2*x) = 0 → ImageSet of n*pi/2",
          "[10d][solveset][trig]") {
    auto x = symbol("x");
    auto s = solveset(sin(integer(2) * x), x);
    REQUIRE(s->kind() == SetKind::ImageSet);
}

TEST_CASE("solveset: tan(x) = 0 → ImageSet of n*pi",
          "[10d][solveset][trig]") {
    auto x = symbol("x");
    auto s = solveset(tan(x), x);
    REQUIRE(s->kind() == SetKind::ImageSet);
}

// ----- reduce_inequalities ---------------------------------------------------

TEST_CASE("reduce_inequalities: single x < 5",
          "[10d][reduce_inequalities]") {
    auto x = symbol("x");
    auto rel = lt(x, integer(5));
    auto s = reduce_inequalities(rel, x);
    REQUIRE(s->contains(integer(0)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(10)) == std::optional<bool>{false});
}

TEST_CASE("reduce_inequalities: AND combines with intersection",
          "[10d][reduce_inequalities]") {
    auto x = symbol("x");
    // x > 0 AND x < 5
    std::vector<Expr> rels = {gt(x, integer(0)), lt(x, integer(5))};
    auto s = reduce_inequalities(rels, x, true);
    REQUIRE(s->contains(integer(3)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(-1)) == std::optional<bool>{false});
    REQUIRE(s->contains(integer(10)) == std::optional<bool>{false});
}

TEST_CASE("reduce_inequalities: OR combines with union",
          "[10d][reduce_inequalities]") {
    auto x = symbol("x");
    // x < -5 OR x > 5 — outside [-5, 5]
    std::vector<Expr> rels = {lt(x, integer(-5)), gt(x, integer(5))};
    auto s = reduce_inequalities(rels, x, false);
    REQUIRE(s->contains(integer(-10)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(10)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(0)) == std::optional<bool>{false});
}

// ----- Pythagorean Diophantine -----------------------------------------------

TEST_CASE("pythagorean_triples: max_z=20 contains (3,4,5), (5,12,13)",
          "[10d][pythagorean]") {
    auto triples = pythagorean_triples(20);
    REQUIRE(!triples.empty());
    bool has_345 = false, has_5_12_13 = false;
    for (auto& t : triples) {
        if (t[0] == integer(3) && t[1] == integer(4) && t[2] == integer(5)) {
            has_345 = true;
        }
        if (t[0] == integer(5) && t[1] == integer(12) && t[2] == integer(13)) {
            has_5_12_13 = true;
        }
    }
    REQUIRE(has_345);
    REQUIRE(has_5_12_13);
    // Verify x² + y² = z² for each.
    for (auto& t : triples) {
        REQUIRE(pow(t[0], integer(2)) + pow(t[1], integer(2))
                == pow(t[2], integer(2)));
    }
}

// ----- Gröbner basis ---------------------------------------------------------

TEST_CASE("groebner: simple linear system",
          "[10d][groebner]") {
    auto x = symbol("x");
    auto y = symbol("y");
    // x - y, x + y - 2 → unique solution (1, 1).
    auto basis = groebner({x - y, x + y - integer(2)}, {x, y});
    REQUIRE(!basis.empty());
}

TEST_CASE("groebner: x²+y²=1 ∧ x=y returns triangular form",
          "[10d][groebner]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto basis = groebner({pow(x, integer(2)) + pow(y, integer(2)) - integer(1),
                            x - y}, {x, y});
    REQUIRE(!basis.empty());
}

TEST_CASE("nonlinsolve_groebner: 2-variable polynomial system",
          "[10d][nonlinsolve_groebner]") {
    auto x = symbol("x");
    auto y = symbol("y");
    // x*y = 1, x + y = 3 → roots of t² - 3t + 1, two solutions.
    auto sols = nonlinsolve_groebner(
        {x * y - integer(1), x + y - integer(3)}, {x, y});
    REQUIRE(!sols.empty());
}
