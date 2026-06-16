// expand() tests.
//
// References:
//   sympy/core/tests/test_expand.py        — distribution patterns
//   sympy/core/tests/test_arit.py          — (x+1)**2 expansion via Pow
//   MATLAB Symbolic Math Toolbox §3 — expand()

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/add.hpp>
#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/miscellaneous.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Distribution ----------------------------------------------------------

TEST_CASE("expand: x*(y + z) = x*y + x*z", "[1i][expand]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    auto e = expand(x * (y + z));
    REQUIRE(e->type_id() == TypeId::Add);
    REQUIRE(e->args().size() == 2);
}

TEST_CASE("expand: (a + b)*(c + d) = a*c + a*d + b*c + b*d", "[1i][expand]") {
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    auto d = symbol("d");
    auto e = expand((a + b) * (c + d));
    REQUIRE(e->type_id() == TypeId::Add);
    REQUIRE(e->args().size() == 4);
}

TEST_CASE("expand: leaves atomic and Number unchanged", "[1i][expand]") {
    auto x = symbol("x");
    REQUIRE(expand(x) == x);
    REQUIRE(expand(integer(5)) == integer(5));
    REQUIRE(expand(S::Pi()) == S::Pi());
}

TEST_CASE("expand: leaves expanded Add unchanged", "[1i][expand]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = x + y;
    REQUIRE(expand(e) == e);
}

// ----- Integer-power expansion -----------------------------------------------

TEST_CASE("expand: (x + 1)^2 = x^2 + 2*x + 1", "[1i][expand][pow]") {
    auto x = symbol("x");
    auto e = expand(pow(x + integer(1), integer(2)));
    REQUIRE(e->type_id() == TypeId::Add);
    REQUIRE(e->args().size() == 3);
}

TEST_CASE("expand: (x - 1)^2 = x^2 - 2*x + 1", "[1i][expand][pow]") {
    auto x = symbol("x");
    auto e = expand(pow(x - integer(1), integer(2)));
    REQUIRE(e->type_id() == TypeId::Add);
    REQUIRE(e->args().size() == 3);
}

TEST_CASE("expand: (x + 1)^3 = x^3 + 3*x^2 + 3*x + 1", "[1i][expand][pow]") {
    auto x = symbol("x");
    auto e = expand(pow(x + integer(1), integer(3)));
    REQUIRE(e->type_id() == TypeId::Add);
    REQUIRE(e->args().size() == 4);
}

TEST_CASE("expand: (a + b)^0 = 1", "[1i][expand][pow]") {
    auto a = symbol("a");
    auto b = symbol("b");
    auto e = expand(pow(a + b, integer(0)));
    REQUIRE(e == S::One());
}

TEST_CASE("expand: (a + b)^1 = a + b", "[1i][expand][pow]") {
    auto a = symbol("a");
    auto b = symbol("b");
    auto e = expand(pow(a + b, integer(1)));
    REQUIRE(e == a + b);
}

TEST_CASE("expand: leaves negative-int Pow unevaluated", "[1i][expand][pow]") {
    // Phase 1i doesn't expand 1/(x+1)^2 into a partial fraction. Verify it stays.
    auto x = symbol("x");
    auto e = expand(pow(x + integer(1), integer(-2)));
    REQUIRE(e->type_id() == TypeId::Pow);
}

// ----- Recursion -------------------------------------------------------------

TEST_CASE("expand: descends into nested Mul", "[1i][expand]") {
    auto x = symbol("x");
    auto y = symbol("y");
    // (x + y) * (x*(x + 1))  → expand inner first, then distribute outer
    auto inner = x * (x + integer(1));   // x*(x+1) = x^2 + x after inner expand
    auto e = expand((x + y) * inner);
    REQUIRE(e->type_id() == TypeId::Add);
}

// ----- Difference of squares & cubes ---------------------------------------

TEST_CASE("expand: (x + y)*(x - y) = x^2 - y^2", "[1i][expand][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = expand((x + y) * (x - y));
    REQUIRE(oracle.equivalent(e->str(), "x**2 - y**2"));
}

TEST_CASE("expand: matches SymPy on classic identities", "[1i][expand][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");

    SECTION("(x + 1)^2") {
        auto e = expand(pow(x + integer(1), integer(2)));
        REQUIRE(oracle.equivalent(e->str(), "x**2 + 2*x + 1"));
    }
    SECTION("(x + 1)^3") {
        auto e = expand(pow(x + integer(1), integer(3)));
        REQUIRE(oracle.equivalent(e->str(), "x**3 + 3*x**2 + 3*x + 1"));
    }
    SECTION("(x + 1)*(x + 2)*(x + 3)") {
        auto e = expand((x + integer(1)) * (x + integer(2)) * (x + integer(3)));
        REQUIRE(oracle.equivalent(e->str(), "x**3 + 6*x**2 + 11*x + 6"));
    }
    SECTION("(x + y)^4 binomial") {
        auto e = expand(pow(x + y, integer(4)));
        REQUIRE(oracle.equivalent(e->str(),
                                  "x**4 + 4*x**3*y + 6*x**2*y**2 + 4*x*y**3 + y**4"));
    }
    SECTION("(a + b + c)^2 multinomial") {
        auto a = symbol("a");
        auto b = symbol("b");
        auto c = symbol("c");
        auto e = expand(pow(a + b + c, integer(2)));
        REQUIRE(oracle.equivalent(
            e->str(), "a**2 + 2*a*b + 2*a*c + b**2 + 2*b*c + c**2"));
    }
}

TEST_CASE("expand: matches SymPy on polynomial multiplication", "[1i][expand][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");

    auto e = expand((pow(x, integer(2)) + integer(2) * x + integer(3))
                   * (x - integer(1)));
    REQUIRE(oracle.equivalent(e->str(), "x**3 + x**2 + x - 3"));
}

TEST_CASE("expand: subs after expand matches SymPy", "[1i][expand][subs][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");

    // (x + 1)^3 expanded, then x = 2 -> 27
    auto e = expand(pow(x + integer(1), integer(3)));
    REQUIRE(oracle.equivalent(e->str(), "x**3 + 3*x**2 + 3*x + 1"));
    // Substitution-after-expansion still works:
    auto sub_result = subs(e, x, integer(2));
    REQUIRE(sub_result == integer(27));
}

// EXPAND-POWBASE-1: expand distributes a power over a product base,
// (a·b)^n → a^n·b^n, folding numeric factors — e.g. (2x)² → 4x². Previously the
// power-of-product was left intact (only powers of an Add expanded).
TEST_CASE("expand: power of a product distributes (EXPAND-POWBASE-1)",
          "[1i][expand][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    // (2x)² → 4x², (2x)³ → 8x³.
    REQUIRE(oracle.equivalent(
        expand(pow(integer(2) * x, integer(2)))->str(), "4*x**2"));
    REQUIRE(oracle.equivalent(
        expand(pow(integer(2) * x, integer(3)))->str(), "8*x**3"));
    // (xy)² → x²y², (3xy)² → 9x²y².
    REQUIRE(oracle.equivalent(expand(pow(x * y, integer(2)))->str(),
                              "x**2*y**2"));
    REQUIRE(oracle.equivalent(
        expand(pow(integer(3) * x * y, integer(2)))->str(), "9*x**2*y**2"));
    // (2xy²)³ → 8x³y⁶.
    REQUIRE(oracle.equivalent(
        expand(pow(integer(2) * x * pow(y, integer(2)), integer(3)))->str(),
        "8*x**3*y**6"));
    // Negative integer exponent distributes too: (2x)⁻² → 1/(4x²).
    REQUIRE(oracle.equivalent(
        expand(pow(integer(2) * x, integer(-2)))->str(), "1/(4*x**2)"));
    // (−x)² → x².
    REQUIRE(oracle.equivalent(
        expand(pow(mul(S::NegativeOne(), x), integer(2)))->str(), "x**2"));
    // Composed with an Add factor: ((x+1)y)² → x²y² + 2xy² + y².
    REQUIRE(oracle.equivalent(
        expand(pow((x + integer(1)) * y, integer(2)))->str(),
        "x**2*y**2 + 2*x*y**2 + y**2"));
    // A non-integer exponent over a possibly-negative factor is NOT split.
    REQUIRE(expand(pow(integer(2) * x, rational(1, 2)))->str()
            == pow(integer(2) * x, rational(1, 2))->str());
}

TEST_CASE("expand: log of a positive product/power splits (ASSUME-4)",
          "[1i][expand][assumptions][regression]") {
    auto p = symbol("p", AssumptionMask{}.set_positive(true));
    auto q = symbol("q", AssumptionMask{}.set_positive(true));
    auto z = symbol("z");  // generic
    // log(p·q) → log(p)+log(q); log(p³) → 3·log(p) — only when positive.
    REQUIRE(expand(log(p * q)) == log(p) + log(q));
    REQUIRE(expand(log(pow(p, integer(3)))) == integer(3) * log(p));
    // A non-positive (generic) factor blocks the split.
    REQUIRE(expand(log(p * z))->str() == "log(p*z)");
    REQUIRE(expand(log(z * symbol("w")))->str() == "log(w*z)");
}

// EXPAND-LOG-FRACPOW-1: log(bᵉ) = e·log(b) is branch-safe for a generic (not
// necessarily positive) base when e is rational with −1 < e < 1, so it is
// extracted even without a positivity assumption — log(√x) → log(x)/2. Exponents
// with |e| ≥ 1 (log(x²), log(1/x), log(x^(3/2))) stay intact, matching SymPy.
TEST_CASE("expand: log of a fractional power of a generic base (EXPAND-LOG-FRACPOW-1)",
          "[1i][expand][regression]") {
    auto x = symbol("x");  // generic, no assumptions
    auto half = rational(1, 2);
    // |e| < 1 extracts on a generic base.
    REQUIRE(expand(log(pow(x, half))) == half * log(x));
    REQUIRE(expand(log(sqrt(x))) == half * log(x));
    REQUIRE(expand(log(pow(x, rational(1, 3)))) == rational(1, 3) * log(x));
    REQUIRE(expand(log(pow(x, rational(2, 3)))) == rational(2, 3) * log(x));
    REQUIRE(expand(log(pow(x, rational(-1, 2)))) == rational(-1, 2) * log(x));
    // |e| ≥ 1 (and symbolic e) stay intact — not branch-safe in general.
    REQUIRE(expand(log(pow(x, integer(2))))->str() == "log(x**2)");
    REQUIRE(expand(log(pow(x, integer(-1))))->str() == "log(x**(-1))");
    REQUIRE(expand(log(pow(x, rational(3, 2))))->str() == "log(x**(3/2))");
    REQUIRE(expand(log(pow(x, symbol("n"))))->str() == "log(x**n)");
}
