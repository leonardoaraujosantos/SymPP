// expand() tests.
//
// References:
//   sympy/core/tests/test_expand.py        — distribution patterns
//   sympy/core/tests/test_arit.py          — (x+1)**2 expansion via Pow
//   MATLAB Symbolic Math Toolbox §3 — expand()

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/add.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>

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
