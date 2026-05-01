// Poly tests.
//
// References:
//   sympy/polys/tests/test_polytools.py — Poly.degree, .coeffs, .as_expr
//   roots: sympy/polys/polyroots.py

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/polys/poly.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Construction & accessors ----------------------------------------------

TEST_CASE("Poly: degree from expression", "[4][poly]") {
    auto x = symbol("x");
    auto e = pow(x, integer(3)) + integer(2) * pow(x, integer(2)) - integer(7) * x + integer(5);
    Poly p(e, x);
    REQUIRE(p.degree() == 3);
}

TEST_CASE("Poly: coefficients in increasing degree", "[4][poly]") {
    auto x = symbol("x");
    auto e = pow(x, integer(3)) + integer(2) * pow(x, integer(2))
             - integer(7) * x + integer(5);
    Poly p(e, x);
    REQUIRE(p.coeffs().size() == 4);
    REQUIRE(p.coeffs()[0] == integer(5));
    REQUIRE(p.coeffs()[1] == integer(-7));
    REQUIRE(p.coeffs()[2] == integer(2));
    REQUIRE(p.coeffs()[3] == integer(1));
}

TEST_CASE("Poly: leading_coeff and constant_term", "[4][poly]") {
    auto x = symbol("x");
    auto e = integer(3) * pow(x, integer(2)) + integer(7);
    Poly p(e, x);
    REQUIRE(p.leading_coeff() == integer(3));
    REQUIRE(p.constant_term() == integer(7));
}

// ----- Arithmetic ------------------------------------------------------------

TEST_CASE("Poly: addition", "[4][poly]") {
    auto x = symbol("x");
    Poly a(pow(x, integer(2)) + integer(1), x);
    Poly b(integer(2) * x, x);
    Poly c = a + b;
    REQUIRE(c.degree() == 2);
    REQUIRE(c.coeffs()[0] == integer(1));
    REQUIRE(c.coeffs()[1] == integer(2));
    REQUIRE(c.coeffs()[2] == integer(1));
}

TEST_CASE("Poly: multiplication", "[4][poly]") {
    auto x = symbol("x");
    // (x + 1) * (x - 1) = x^2 - 1
    Poly a(x + integer(1), x);
    Poly b(x - integer(1), x);
    Poly c = a * b;
    REQUIRE(c.degree() == 2);
    REQUIRE(c.coeffs()[0] == integer(-1));
    REQUIRE(c.coeffs()[1] == integer(0));
    REQUIRE(c.coeffs()[2] == integer(1));
}

// ----- eval ------------------------------------------------------------------

TEST_CASE("Poly: eval at integer", "[4][poly]") {
    auto x = symbol("x");
    Poly p(pow(x, integer(2)) + integer(2) * x + integer(3), x);
    // p(2) = 4 + 4 + 3 = 11
    REQUIRE(p.eval(integer(2)) == integer(11));
}

// ----- as_expr roundtrip -----------------------------------------------------

TEST_CASE("Poly: as_expr round-trips through SymPy", "[4][poly][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(3)) - integer(7) * x + integer(2);
    Poly p(e, x);
    REQUIRE(oracle.equivalent(p.as_expr()->str(), "x**3 - 7*x + 2"));
}

// ----- Roots -----------------------------------------------------------------

TEST_CASE("Poly: linear root", "[4][poly][roots]") {
    auto x = symbol("x");
    Poly p(integer(3) * x - integer(7), x);
    auto roots = p.roots();
    REQUIRE(roots.size() == 1);
    REQUIRE(roots[0] == rational(7, 3));
}

TEST_CASE("Poly: quadratic roots", "[4][poly][roots]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^2 - 5x + 6 = (x - 2)(x - 3) -> roots 2, 3
    Poly p(pow(x, integer(2)) - integer(5) * x + integer(6), x);
    auto roots = p.roots();
    REQUIRE(roots.size() == 2);
    // The roots come in (-b ± sqrt(disc))/(2a) order; check both numeric values
    // collapse via the oracle to the expected set.
    REQUIRE(oracle.equivalent(roots[0]->str(), "3"));
    REQUIRE(oracle.equivalent(roots[1]->str(), "2"));
}

TEST_CASE("Poly: quadratic with irrational roots", "[4][poly][roots][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^2 - 2 -> roots ±sqrt(2)
    Poly p(pow(x, integer(2)) - integer(2), x);
    auto roots = p.roots();
    REQUIRE(roots.size() == 2);
    REQUIRE(oracle.equivalent(roots[0]->str(), "sqrt(2)"));
    REQUIRE(oracle.equivalent(roots[1]->str(), "-sqrt(2)"));
}

TEST_CASE("Poly: degree 3+ roots return empty (deferred)", "[4][poly][roots]") {
    auto x = symbol("x");
    Poly p(pow(x, integer(3)) - integer(1), x);
    auto roots = p.roots();
    REQUIRE(roots.empty());
}
