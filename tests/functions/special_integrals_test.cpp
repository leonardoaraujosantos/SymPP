// Fresnel integrals & generalized exponential integral — derivatives and
// special values cross-checked against SymPy.

#include <string>

#include <catch2/catch_test_macros.hpp>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/special.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {
bool deriv_matches(const Expr& f, const std::string& expr, const std::string& var) {
    auto& oracle = sympp::testing::Oracle::instance();
    std::string ref = oracle.diff(expr, var);  // SymPy d/dvar
    return oracle.equivalent(diff(f, symbol(var))->str(), ref);
}
}  // namespace

TEST_CASE("Fresnel integrals: derivatives and special values", "[special][oracle]") {
    auto x = symbol("x");
    REQUIRE(deriv_matches(fresnels(x), "fresnels(x)", "x"));  // S'(x)=sin(pi x^2/2)
    REQUIRE(deriv_matches(fresnelc(x), "fresnelc(x)", "x"));  // C'(x)=cos(pi x^2/2)

    REQUIRE(fresnels(integer(0)) == integer(0));
    REQUIRE(fresnelc(integer(0)) == integer(0));
    REQUIRE(fresnels(S::Infinity()) == rational(1, 2));
    REQUIRE(fresnelc(S::Infinity()) == rational(1, 2));
    // Odd symmetry: S(-x) = -S(x).
    REQUIRE(fresnels(mul(integer(-1), x)) == mul(integer(-1), fresnels(x)));
}

TEST_CASE("generalized exponential integral", "[special][oracle]") {
    auto z = symbol("z");
    auto n = symbol("n");
    // d/dz Eₙ(z) = -E_{n-1}(z).
    REQUIRE(deriv_matches(expint(n, z), "expint(n, z)", "z"));
    // E₀(z) = e^{-z}/z.
    REQUIRE(expint(integer(0), z) == mul(exp(mul(integer(-1), z)), pow(z, integer(-1))));
    // Eₙ(∞) = 0.
    REQUIRE(expint(integer(2), S::Infinity()) == integer(0));
}
