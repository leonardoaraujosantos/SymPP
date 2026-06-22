// Berlekamp–Zassenhaus univariate ℤ factorization — product reconstruction and
// factor counts cross-checked against SymPy.

#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <sympp/core/expand.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/polys/zassenhaus.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {
auto& oracle() { return sympp::testing::Oracle::instance(); }

int sympy_factor_count(const Expr& e, const Expr& var) {
    nlohmann::json req = {{"op", "polyfactor"}, {"expr", e->str()}, {"var", var->str()}};
    auto r = oracle().send(req);
    REQUIRE(r.ok);
    return std::stoi(r.result_str());
}
// The factors multiply back to the input, each is irreducible, and the count
// matches SymPy's factorization.
void check_factors(const Expr& p, const Expr& x) {
    auto facs = factor_zassenhaus(p, x);
    Expr prod = integer(1);
    for (const auto& f : facs) prod = mul(prod, f);
    REQUIRE(oracle().equivalent(prod->str(), p->str()));
    REQUIRE(static_cast<int>(facs.size()) == sympy_factor_count(p, x));
    for (const auto& f : facs) REQUIRE(sympy_factor_count(f, x) == 1);  // irreducible
}
}  // namespace

TEST_CASE("Zassenhaus factorization matches SymPy", "[zassenhaus][oracle]") {
    auto x = symbol("x");
    check_factors(pow(x, integer(2)) - integer(1), x);            // (x-1)(x+1)
    check_factors(pow(x, integer(4)) - integer(1), x);            // (x-1)(x+1)(x²+1)
    check_factors(pow(x, integer(3)) - x, x);                     // x(x-1)(x+1)
    check_factors(pow(x, integer(2)) - integer(5) * x + integer(6), x);   // (x-2)(x-3)
    check_factors(pow(x, integer(4)) + pow(x, integer(2)) + integer(1), x);  // (x²+x+1)(x²-x+1)
    check_factors(pow(x, integer(6)) - integer(1), x);
    check_factors(pow(x, integer(2)) - integer(9), x);  // (x-3)(x+3)
    // repeated factors
    check_factors(expand(pow(x + integer(1), integer(2)) * (x - integer(2))), x);
    check_factors(expand(pow(x - integer(1), integer(3))), x);
}

// Higher-degree / larger-coefficient inputs whose Landau–Mignotte modulus needs
// several Hensel-lift steps from a small prime (p^l ≫ p) before recombination.
TEST_CASE("Zassenhaus: Hensel lifting on larger inputs", "[zassenhaus][oracle][hensel]") {
    auto x = symbol("x");
    auto X = [&](int k) { return pow(x, integer(k)); };
    // (x²+1)(x²+2)(x²+3): three irreducible quadratics, coefficients up to 11.
    check_factors(expand(mul(mul(X(2) + integer(1), X(2) + integer(2)), X(2) + integer(3))), x);
    // (x²+x+1)(x²−x+1)(x²+1) = x⁶+x⁴+x²+1 region — distinct cyclotomic factors.
    check_factors(expand(mul(mul(X(2) + x + integer(1), X(2) - x + integer(1)),
                             X(2) + integer(1))),
                  x);
    // x⁸ − 1: splits into (x−1)(x+1)(x²+1)(x⁴+1).
    check_factors(X(8) - integer(1), x);
    // (3x²+2)(x³−x+5): large leading/content interplay, degree 5.
    check_factors(expand(mul(integer(3) * X(2) + integer(2), X(3) - x + integer(5))), x);
    // (x²−2x+3)² · (x+4): repeated irreducible quadratic with multiplicity.
    check_factors(expand(mul(pow(X(2) - integer(2) * x + integer(3), integer(2)),
                             x + integer(4))),
                  x);
    // x¹⁰ − 1: ten-th cyclotomic split, exercises a wider recombination search.
    check_factors(X(10) - integer(1), x);
}

TEST_CASE("Zassenhaus: irreducible stays whole", "[zassenhaus]") {
    auto x = symbol("x");
    auto facs = factor_zassenhaus(pow(x, integer(2)) - integer(2), x);  // irreducible over ℤ
    REQUIRE(facs.size() == 1);
    auto fc = factor_zassenhaus(pow(x, integer(2)) + integer(1), x);
    REQUIRE(fc.size() == 1);
}
