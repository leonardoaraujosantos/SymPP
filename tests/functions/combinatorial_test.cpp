// factorial / binomial / gamma / loggamma tests.
//
// References:
//   sympy/functions/combinatorial/tests/test_comb_factorials.py
//   sympy/functions/special/tests/test_gamma_functions.py

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/functions/combinatorial.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- factorial -------------------------------------------------------------

TEST_CASE("factorial: small integer values", "[3i][factorial]") {
    REQUIRE(factorial(integer(0)) == integer(1));
    REQUIRE(factorial(integer(1)) == integer(1));
    REQUIRE(factorial(integer(5)) == integer(120));
    REQUIRE(factorial(integer(10)) == integer(3628800));
}

TEST_CASE("factorial: large integer", "[3i][factorial][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = factorial(integer(20));
    REQUIRE(oracle.equivalent(e->str(), "factorial(20)"));
    // 20! = 2432902008176640000
    REQUIRE(e == integer("2432902008176640000"));
}

TEST_CASE("factorial: stays unevaluated on a symbol", "[3i][factorial]") {
    auto x = symbol("x");
    REQUIRE(factorial(x)->type_id() == TypeId::Function);
}

// ----- binomial --------------------------------------------------------------

TEST_CASE("binomial: classic values", "[3i][binomial]") {
    REQUIRE(binomial(integer(5), integer(0)) == integer(1));
    REQUIRE(binomial(integer(5), integer(1)) == integer(5));
    REQUIRE(binomial(integer(5), integer(2)) == integer(10));
    REQUIRE(binomial(integer(5), integer(3)) == integer(10));
    REQUIRE(binomial(integer(5), integer(5)) == integer(1));
}

TEST_CASE("binomial: k > n returns 0", "[3i][binomial]") {
    REQUIRE(binomial(integer(3), integer(5)) == integer(0));
}

TEST_CASE("binomial: stays symbolic for symbolic args", "[3i][binomial]") {
    auto n = symbol("n");
    auto k = symbol("k");
    auto e = binomial(n, k);
    REQUIRE(e->type_id() == TypeId::Function);
    REQUIRE(e->str() == "binomial(n, k)");
}

TEST_CASE("binomial(n, 0) = 1", "[3i][binomial]") {
    auto n = symbol("n");
    REQUIRE(binomial(n, S::Zero()) == S::One());
}

// ----- gamma -----------------------------------------------------------------

TEST_CASE("gamma: positive integer reduces to factorial", "[3i][gamma]") {
    // gamma(n) = (n-1)!
    REQUIRE(gamma(integer(1)) == integer(1));
    REQUIRE(gamma(integer(2)) == integer(1));
    REQUIRE(gamma(integer(3)) == integer(2));
    REQUIRE(gamma(integer(4)) == integer(6));
    REQUIRE(gamma(integer(5)) == integer(24));
}

TEST_CASE("gamma: numeric Float matches MPFR", "[3i][gamma]") {
    auto e = gamma(float_value(1.5));
    REQUIRE(e->type_id() == TypeId::Float);
}

TEST_CASE("gamma: high-precision Float matches SymPy",
          "[3i][gamma][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = gamma(float_value("1.5", 30));
    auto sympy_value = oracle.evalf("gamma(Rational(3, 2))", 30);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

TEST_CASE("gamma: stays symbolic for symbolic arg", "[3i][gamma]") {
    auto x = symbol("x");
    REQUIRE(gamma(x)->str() == "gamma(x)");
}

// ----- loggamma --------------------------------------------------------------

TEST_CASE("loggamma: classic values", "[3i][loggamma]") {
    REQUIRE(loggamma(integer(1)) == integer(0));   // log(0!) = 0
    REQUIRE(loggamma(integer(2)) == integer(0));   // log(1!) = 0
}

TEST_CASE("loggamma: numeric Float matches mpfr_lngamma", "[3i][loggamma]") {
    auto e = loggamma(float_value(5.0));
    REQUIRE(e->type_id() == TypeId::Float);
}

// ----- Substitution ----------------------------------------------------------

TEST_CASE("factorial: subs collapses", "[3i][subs]") {
    auto x = symbol("x");
    REQUIRE(subs(factorial(x), x, integer(6)) == integer(720));
}

TEST_CASE("binomial: subs collapses", "[3i][subs]") {
    auto n = symbol("n");
    REQUIRE(subs(binomial(n, integer(2)), n, integer(10)) == integer(45));
}

// ----- Oracle structural -----------------------------------------------------

TEST_CASE("combinatorial: structural output matches SymPy",
          "[3i][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(factorial(integer(5))->str(), "factorial(5)"));
    REQUIRE(oracle.equivalent(binomial(integer(10), integer(3))->str(),
                              "binomial(10, 3)"));
    REQUIRE(oracle.equivalent(gamma(integer(5))->str(), "gamma(5)"));
}
