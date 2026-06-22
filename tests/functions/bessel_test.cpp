// Bessel functions J/Y/I/K — special values, half-integer closed forms and
// derivative recurrences, cross-checked against SymPy.

#include <catch2/catch_test_macros.hpp>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/bessel.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {
auto& oracle() { return sympp::testing::Oracle::instance(); }
void eq(const Expr& mine, const std::string& sympy) {
    REQUIRE(oracle().equivalent(mine->str(), sympy));
}
}  // namespace

TEST_CASE("bessel: special values at 0", "[bessel]") {
    auto z = symbol("z");
    REQUIRE(besselj(integer(0), S::Zero()) == integer(1));
    REQUIRE(besselj(integer(2), S::Zero()) == S::Zero());
    REQUIRE(besseli(integer(0), S::Zero()) == integer(1));
    REQUIRE(besseli(integer(3), S::Zero()) == S::Zero());
    // Generic order stays symbolic.
    REQUIRE(besselj(symbol("n"), z)->str().rfind("besselj", 0) == 0);
}

TEST_CASE("bessel: half-integer closed forms match SymPy", "[bessel][oracle]") {
    auto z = symbol("z");
    auto h = rational(1, 2), mh = rational(-1, 2);
    eq(besselj(h, z), "sqrt(2)*sin(z)/(sqrt(pi)*sqrt(z))");
    eq(besselj(mh, z), "sqrt(2)*cos(z)/(sqrt(pi)*sqrt(z))");
    eq(besseli(h, z), "sqrt(2)*sinh(z)/(sqrt(pi)*sqrt(z))");
    eq(besseli(mh, z), "sqrt(2)*cosh(z)/(sqrt(pi)*sqrt(z))");
    eq(besselk(h, z), "sqrt(2)*sqrt(pi)*exp(-z)/(2*sqrt(z))");
    eq(besselk(mh, z), "sqrt(2)*sqrt(pi)*exp(-z)/(2*sqrt(z))");  // K even in ν
    eq(bessely(h, z), "-sqrt(2)*cos(z)/(sqrt(pi)*sqrt(z))");
    eq(bessely(mh, z), "sqrt(2)*sin(z)/(sqrt(pi)*sqrt(z))");
}

TEST_CASE("bessel: z-derivative recurrences match SymPy", "[bessel][oracle]") {
    auto z = symbol("z");
    auto n = symbol("n");
    eq(diff(besselj(n, z), z), "besselj(n - 1, z)/2 - besselj(n + 1, z)/2");
    eq(diff(bessely(n, z), z), "bessely(n - 1, z)/2 - bessely(n + 1, z)/2");
    eq(diff(besseli(n, z), z), "besseli(n - 1, z)/2 + besseli(n + 1, z)/2");
    eq(diff(besselk(n, z), z), "-besselk(n - 1, z)/2 - besselk(n + 1, z)/2");
}
