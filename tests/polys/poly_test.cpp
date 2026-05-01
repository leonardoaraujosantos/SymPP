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

// ----- Polynomial division ---------------------------------------------------

namespace {
struct DivOracle { std::string quotient, remainder; };
[[nodiscard]] DivOracle oracle_div(std::string_view p, std::string_view q,
                                   std::string_view var) {
    auto& oracle = Oracle::instance();
    auto resp = oracle.send({{"op", "div"}, {"p", p}, {"q", q}, {"var", var}});
    REQUIRE(resp.ok);
    return DivOracle{resp.raw.at("quotient").get<std::string>(),
                     resp.raw.at("remainder").get<std::string>()};
}
}  // namespace

TEST_CASE("Poly: divmod exact (x^2 - 1) / (x - 1)", "[4][poly][divmod][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    Poly p(pow(x, integer(2)) - integer(1), x);
    Poly d(x - integer(1), x);
    auto [q, r] = p.divmod(d);
    auto ref = oracle_div("x**2 - 1", "x - 1", "x");
    REQUIRE(oracle.equivalent(q.as_expr()->str(), ref.quotient));
    REQUIRE(oracle.equivalent(r.as_expr()->str(), ref.remainder));
    REQUIRE(r.is_zero());
}

TEST_CASE("Poly: divmod with remainder (x^2 + 1) / (x + 1)", "[4][poly][divmod][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    Poly p(pow(x, integer(2)) + integer(1), x);
    Poly d(x + integer(1), x);
    auto [q, r] = p.divmod(d);
    auto ref = oracle_div("x**2 + 1", "x + 1", "x");
    REQUIRE(oracle.equivalent(q.as_expr()->str(), ref.quotient));
    REQUIRE(oracle.equivalent(r.as_expr()->str(), ref.remainder));
}

TEST_CASE("Poly: divmod cubic by quadratic", "[4][poly][divmod][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // (x^3 - 6x^2 + 11x - 6) / (x^2 - 3x + 2) = x - 3, remainder 0
    Poly p(pow(x, integer(3)) - integer(6) * pow(x, integer(2))
               + integer(11) * x - integer(6), x);
    Poly d(pow(x, integer(2)) - integer(3) * x + integer(2), x);
    auto [q, r] = p.divmod(d);
    auto ref = oracle_div("x**3 - 6*x**2 + 11*x - 6",
                          "x**2 - 3*x + 2", "x");
    REQUIRE(oracle.equivalent(q.as_expr()->str(), ref.quotient));
    REQUIRE(oracle.equivalent(r.as_expr()->str(), ref.remainder));
}

TEST_CASE("Poly: divmod with rational coefficients", "[4][poly][divmod][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // (2x^3 + 3x^2 - x + 5) / (2x + 1)
    Poly p(integer(2) * pow(x, integer(3)) + integer(3) * pow(x, integer(2))
               - x + integer(5), x);
    Poly d(integer(2) * x + integer(1), x);
    auto [q, r] = p.divmod(d);
    auto ref = oracle_div("2*x**3 + 3*x**2 - x + 5", "2*x + 1", "x");
    REQUIRE(oracle.equivalent(q.as_expr()->str(), ref.quotient));
    REQUIRE(oracle.equivalent(r.as_expr()->str(), ref.remainder));
}

TEST_CASE("Poly: divmod degree(divisor) > degree(dividend)", "[4][poly][divmod]") {
    auto x = symbol("x");
    // (x + 1) / (x^2 + 1)  -> q=0, r=x+1
    Poly p(x + integer(1), x);
    Poly d(pow(x, integer(2)) + integer(1), x);
    auto [q, r] = p.divmod(d);
    REQUIRE(q.is_zero());
    REQUIRE(r.as_expr() == p.as_expr());
}

TEST_CASE("Poly: divmod by zero throws", "[4][poly][divmod]") {
    auto x = symbol("x");
    Poly p(pow(x, integer(2)) + x, x);
    Poly z(std::vector<Expr>{S::Zero()}, x);
    REQUIRE_THROWS_AS(p.divmod(z), std::invalid_argument);
}

TEST_CASE("Poly: monic normalization", "[4][poly][divmod]") {
    auto x = symbol("x");
    Poly p(integer(3) * pow(x, integer(2)) + integer(6) * x + integer(9), x);
    auto m = p.monic();
    // Now leading coeff is 1, others divided by 3
    REQUIRE(m.leading_coeff() == S::One());
    REQUIRE(m.coeffs()[1] == integer(2));
    REQUIRE(m.coeffs()[0] == integer(3));
}

TEST_CASE("Poly: diff lowers degree", "[4][poly][diff]") {
    auto x = symbol("x");
    // (x^3 + 2x^2 - 7x + 5)' = 3x^2 + 4x - 7
    Poly p(pow(x, integer(3)) + integer(2) * pow(x, integer(2))
               - integer(7) * x + integer(5), x);
    auto dp = p.diff();
    REQUIRE(dp.degree() == 2);
    REQUIRE(dp.coeffs()[2] == integer(3));
    REQUIRE(dp.coeffs()[1] == integer(4));
    REQUIRE(dp.coeffs()[0] == integer(-7));
}
