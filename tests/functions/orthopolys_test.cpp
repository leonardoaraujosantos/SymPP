// Classical orthogonal polynomials — cross-checked against SymPy.

#include <stdexcept>
#include <string>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/orthopolys.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {
// True iff SymPP's polynomial equals SymPy's reference `sympy_name(n, x)`.
bool matches(const Expr& got, const std::string& sympy_name, long n) {
    auto& oracle = sympp::testing::Oracle::instance();
    std::string ref = sympy_name + "(" + std::to_string(n) + ", x)";
    return oracle.equivalent(got->str(), ref);
}
}  // namespace

TEST_CASE("orthogonal polynomials match SymPy", "[orthopolys][oracle]") {
    auto x = symbol("x");
    for (long n = 0; n <= 8; ++n) {
        REQUIRE(matches(legendre(integer(n), x), "legendre", n));
        REQUIRE(matches(chebyshevt(integer(n), x), "chebyshevt", n));
        REQUIRE(matches(chebyshevu(integer(n), x), "chebyshevu", n));
        REQUIRE(matches(hermite(integer(n), x), "hermite", n));
        REQUIRE(matches(laguerre(integer(n), x), "laguerre", n));
    }
}

TEST_CASE("orthogonal polynomials: low-order closed forms", "[orthopolys]") {
    auto x = symbol("x");
    REQUIRE(legendre(integer(0), x) == integer(1));
    REQUIRE(chebyshevt(integer(1), x) == x);
    // chebyshevu(1, x) = 2x
    REQUIRE(chebyshevu(integer(1), x)->str() == "2*x");
}

TEST_CASE("orthogonal polynomials reject bad degree", "[orthopolys]") {
    auto x = symbol("x");
    REQUIRE_THROWS_AS(legendre(integer(-1), x), std::invalid_argument);
    REQUIRE_THROWS_AS(hermite(symbol("n"), x), std::invalid_argument);
}
