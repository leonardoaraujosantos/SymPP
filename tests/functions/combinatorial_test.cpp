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
#include <sympp/parsing/parser.hpp>

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

// Regression (BINOM-1): binomial(n, 1) = n for any n (SymPy auto-evaluates it;
// SymPP previously kept it symbolic).
TEST_CASE("binomial(n, 1) = n", "[3i][binomial][regression]") {
    auto n = symbol("n");
    REQUIRE(binomial(n, S::One()) == n);
    REQUIRE(binomial(integer(5), S::One()) == integer(5));
    // binomial(n, 2) is not a special identity — must stay symbolic.
    REQUIRE(binomial(n, integer(2))->type_id() == TypeId::Function);
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

// GAMMA-1: half-integer arguments stayed symbolic; SymPy gives closed forms
// with sqrt(pi). gamma(p/2) = C·sqrt(pi) via the duplication recurrence.
TEST_CASE("gamma: half-integer reduces to a multiple of sqrt(pi)",
          "[3i][gamma][oracle][regression]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(gamma(rational(1, 2))->str(), "sqrt(pi)"));
    REQUIRE(oracle.equivalent(gamma(rational(3, 2))->str(), "sqrt(pi)/2"));
    REQUIRE(oracle.equivalent(gamma(rational(5, 2))->str(), "3*sqrt(pi)/4"));
    REQUIRE(oracle.equivalent(gamma(rational(7, 2))->str(), "15*sqrt(pi)/8"));
}

TEST_CASE("gamma: negative half-integer reduces to a multiple of sqrt(pi)",
          "[3i][gamma][oracle][regression]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(gamma(rational(-1, 2))->str(), "-2*sqrt(pi)"));
    REQUIRE(oracle.equivalent(gamma(rational(-3, 2))->str(), "4*sqrt(pi)/3"));
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

// ----- fibonacci / catalan (FIB-CAT) -----------------------------------------

TEST_CASE("fibonacci: integer values", "[3i][fibonacci]") {
    REQUIRE(fibonacci(integer(0)) == integer(0));
    REQUIRE(fibonacci(integer(1)) == integer(1));
    REQUIRE(fibonacci(integer(2)) == integer(1));
    REQUIRE(fibonacci(integer(10)) == integer(55));
    REQUIRE(fibonacci(integer(20)) == integer(6765));
    // Symbolic / negative arguments stay unevaluated.
    auto x = symbol("x");
    REQUIRE(fibonacci(x)->type_id() == TypeId::Function);
    REQUIRE(fibonacci(integer(-1))->type_id() == TypeId::Function);
}

TEST_CASE("catalan: integer values", "[3i][catalan]") {
    REQUIRE(catalan(integer(0)) == integer(1));
    REQUIRE(catalan(integer(1)) == integer(1));
    REQUIRE(catalan(integer(2)) == integer(2));
    REQUIRE(catalan(integer(3)) == integer(5));
    REQUIRE(catalan(integer(5)) == integer(42));
    REQUIRE(catalan(integer(10)) == integer(16796));
    auto x = symbol("x");
    REQUIRE(catalan(x)->type_id() == TypeId::Function);
}

TEST_CASE("fibonacci/catalan: parse round-trip and subs",
          "[3i][fibonacci][catalan][parser]") {
    auto x = symbol("x");
    REQUIRE(parsing::parse("fibonacci(x)") == fibonacci(x));
    REQUIRE(parsing::parse("catalan(x)") == catalan(x));
    REQUIRE(fibonacci(x)->str() == "fibonacci(x)");
    REQUIRE(catalan(x)->str() == "catalan(x)");
    REQUIRE(subs(fibonacci(x), x, integer(12)) == integer(144));
    REQUIRE(subs(catalan(x), x, integer(4)) == integer(14));
}

TEST_CASE("fibonacci/catalan: values match SymPy", "[3i][fibonacci][catalan][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(fibonacci(integer(15))->str(), "fibonacci(15)"));
    REQUIRE(oracle.equivalent(catalan(integer(7))->str(), "catalan(7)"));
}

// ----- gcd / lcm (GCD-LCM) ---------------------------------------------------

TEST_CASE("gcd: integer values", "[3i][gcd]") {
    REQUIRE(gcd(integer(12), integer(18)) == integer(6));
    REQUIRE(gcd(integer(17), integer(5)) == integer(1));
    REQUIRE(gcd(integer(-12), integer(8)) == integer(4));  // non-negative result
    REQUIRE(gcd(integer(0), integer(5)) == integer(5));
    REQUIRE(gcd(integer(0), integer(0)) == integer(0));
    // Symbolic args stay unevaluated.
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(gcd(x, y)->type_id() == TypeId::Function);
}

TEST_CASE("lcm: integer values", "[3i][lcm]") {
    REQUIRE(lcm(integer(4), integer(6)) == integer(12));
    REQUIRE(lcm(integer(21), integer(6)) == integer(42));
    REQUIRE(lcm(integer(0), integer(5)) == integer(0));
    REQUIRE(lcm(integer(7), integer(5)) == integer(35));
    auto x = symbol("x");
    REQUIRE(lcm(x, integer(4))->type_id() == TypeId::Function);
}

TEST_CASE("gcd/lcm: parse round-trip and subs", "[3i][gcd][lcm][parser]") {
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(parsing::parse("gcd(x, y)") == gcd(x, y));
    REQUIRE(parsing::parse("lcm(x, y)") == lcm(x, y));
    REQUIRE(gcd(x, y)->str() == "gcd(x, y)");
    REQUIRE(lcm(x, y)->str() == "lcm(x, y)");
    REQUIRE(subs(gcd(x, integer(18)), x, integer(12)) == integer(6));
    REQUIRE(subs(lcm(x, integer(6)), x, integer(4)) == integer(12));
}
