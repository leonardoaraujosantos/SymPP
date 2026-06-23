// product tests — closed-form symbolic products ∏_{k=lo}^{hi} f(k).

#include <catch2/catch_test_macros.hpp>

#include <string>

#include <sympp/calculus/product.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// PROD-1: the core closed forms, each cross-checked against SymPy's product().
// The Gamma-ratio / factorial spellings may differ textually from SymPy's, so
// equivalence is checked numerically by the oracle rather than by string match.
TEST_CASE("product: closed forms (PROD-1)", "[6][product][oracle]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto x = symbol("x");
    auto I = [](int i) { return integer(i); };
    auto chk = [&](const Expr& got, const std::string& sympy) {
        REQUIRE(oracle.equivalent(got->str(), sympy));
    };

    // ∏_{k=1}^{n} k = n!  (canonical factorial spelling).
    REQUIRE(product(k, k, I(1), n)->str() == "factorial(n)");
    chk(product(k, k, I(1), n), "product(k,(k,1,n))");
    // Shifted lower bound → Gamma ratio.
    chk(product(k, k, I(2), n), "product(k,(k,2,n))");
    // Constant w.r.t. k: x → x^n.
    chk(product(x, k, I(1), n), "product(x,(k,1,n))");
    // Multiplicative split: ∏ 2k = 2^n·n!.
    chk(product(mul(I(2), k), k, I(1), n), "product(2*k,(k,1,n))");
    // Constant exponent: ∏ k² = (n!)².
    chk(product(pow(k, I(2)), k, I(1), n), "product(k**2,(k,1,n))");
    // Geometric exponent: ∏ 2^k = 2^(n(n+1)/2).
    chk(product(pow(I(2), k), k, I(1), n), "product(2**k,(k,1,n))");
    // ∏ x^k = x^(n(n+1)/2) for a symbolic base.
    chk(product(pow(x, k), k, I(1), n), "product(x**k,(k,1,n))");
    // Linear factor ∏ (k+1) and ∏ (2k+1).
    chk(product(k + I(1), k, I(1), n), "product(k+1,(k,1,n))");
    chk(product(mul(I(2), k) + I(1), k, I(1), n), "product(2*k+1,(k,1,n))");
    // Quadratic that factors: ∏ (k²−1) over k=2..n.
    chk(product(pow(k, I(2)) - I(1), k, I(2), n), "product(k**2-1,(k,2,n))");
    // Telescoping rational: ∏ (k+1)/(k+2) = 2/(n+2).
    chk(product(mul(k + I(1), pow(k + I(2), I(-1))), k, I(1), n),
        "product((k+1)/(k+2),(k,1,n))");
    chk(product(mul(k + I(1), pow(k + I(2), I(-1))), k, I(1), n), "2/(n+2)");
    // Product of two factors: ∏ k(k+1).
    chk(product(mul(k, k + I(1)), k, I(1), n), "product(k*(k+1),(k,1,n))");
}

// PROD-2: numeric ranges fully evaluate, and forms without a closed product
// stay as an unevaluated Product(...) marker (mirroring summation's Sum).
TEST_CASE("product: numeric evaluation and unevaluated marker (PROD-2)",
          "[6][product][oracle]") {
    auto& oracle = Oracle::instance();
    auto k = symbol("k");
    auto n = symbol("n");
    auto I = [](int i) { return integer(i); };

    // ∏_{k=1}^{5} 3 = 3^5 = 243; ∏_{k=1}^{5} k = 120.
    REQUIRE(oracle.equivalent(product(I(3), k, I(1), I(5))->str(), "243"));
    REQUIRE(oracle.equivalent(product(k, k, I(1), I(5))->str(), "120"));
    // Single-term range collapses to the term.
    REQUIRE(oracle.equivalent(product(pow(k, I(2)), k, n, n)->str(), "n**2"));

    // ∏ k^k has no elementary closed form → unevaluated marker, never dropped
    // to its bare term.
    Expr p = product(pow(k, k), k, I(1), n);
    REQUIRE(p->str().find("Product(") != std::string::npos);
}
