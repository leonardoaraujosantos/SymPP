// Add / Mul / Pow / operator-overload tests.
//
// Reference inputs adapted from:
//   sympy/core/tests/test_arit.py::test_arit0
//   sympy/core/tests/test_arit.py::test_div
//   sympy/core/tests/test_arit.py::test_pow_arit
//   sympy/core/tests/test_arit.py::test_mul_add_identity
//
// Phase 1f.1 doesn't yet do same-term collection (x + x → 2*x); the oracle
// equivalence check handles that on the comparison side.

#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/add.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Numeric arithmetic short-circuits -------------------------------------

TEST_CASE("add: pure-numeric collapses to a single Number", "[1][add]") {
    auto e = add(integer(2), integer(3));
    REQUIRE(e->type_id() == TypeId::Integer);
    REQUIRE(e->str() == "5");

    auto e2 = add({integer(1), integer(2), integer(3), integer(4)});
    REQUIRE(e2->str() == "10");

    // Integer + Rational -> Rational
    auto e3 = add(integer(1), rational(1, 2));
    REQUIRE(e3->str() == "3/2");
}

TEST_CASE("mul: pure-numeric collapses to a single Number", "[1][mul]") {
    auto e = mul(integer(2), integer(3));
    REQUIRE(e->type_id() == TypeId::Integer);
    REQUIRE(e->str() == "6");

    auto e2 = mul({integer(2), rational(1, 4)});
    REQUIRE(e2->str() == "1/2");
}

TEST_CASE("add: identity rules", "[1][add]") {
    auto x = symbol("x");
    REQUIRE(add(x, S::Zero()) == x);   // x + 0 = x
    REQUIRE(add(S::Zero(), x) == x);   // 0 + x = x
    REQUIRE(add(S::Zero(), S::Zero()) == S::Zero());
}

TEST_CASE("mul: identity rules", "[1][mul]") {
    auto x = symbol("x");
    REQUIRE(mul(x, S::One()) == x);    // x * 1 = x
    REQUIRE(mul(S::One(), x) == x);    // 1 * x = x
    REQUIRE(mul(x, S::Zero()) == S::Zero());     // x * 0 = 0
    REQUIRE(mul(S::Zero(), x) == S::Zero());     // 0 * x = 0
}

TEST_CASE("add: flattens nested Adds", "[1][add]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    auto e = add(add(x, y), z);
    REQUIRE(e->type_id() == TypeId::Add);
    REQUIRE(e->args().size() == 3);
}

TEST_CASE("mul: flattens nested Muls", "[1][mul]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    auto e = mul(mul(x, y), z);
    REQUIRE(e->type_id() == TypeId::Mul);
    REQUIRE(e->args().size() == 3);
}

TEST_CASE("add: combines all numerics inside flattened Add", "[1][add]") {
    auto x = symbol("x");
    // Add(2, x, 3) should combine 2+3 and keep x
    auto e = add({integer(2), x, integer(3)});
    REQUIRE(e->type_id() == TypeId::Add);
    // expect: 2 args (x, 5)
    REQUIRE(e->args().size() == 2);
}

// ----- Pow auto-eval ----------------------------------------------------------

TEST_CASE("pow: x^0 = 1", "[1][pow]") {
    auto x = symbol("x");
    REQUIRE(pow(x, S::Zero()) == S::One());
}

TEST_CASE("pow: x^1 = x", "[1][pow]") {
    auto x = symbol("x");
    REQUIRE(pow(x, S::One()) == x);
}

TEST_CASE("pow: 1^x = 1", "[1][pow]") {
    auto x = symbol("x");
    REQUIRE(pow(S::One(), x) == S::One());
}

TEST_CASE("pow: 0^positive = 0", "[1][pow]") {
    REQUIRE(pow(S::Zero(), integer(5)) == S::Zero());
}

TEST_CASE("pow: Integer^Integer is computed exactly", "[1][pow]") {
    REQUIRE(pow(integer(2), integer(10))->str() == "1024");
    REQUIRE(pow(integer(3), integer(4))->str() == "81");
    REQUIRE(pow(integer(-2), integer(3))->str() == "-8");
}

TEST_CASE("pow: Integer^negative-Integer becomes Rational", "[1][pow]") {
    auto e = pow(integer(2), integer(-3));
    REQUIRE(e->type_id() == TypeId::Rational);
    REQUIRE(e->str() == "1/8");
}

// Regression (POW-RAT): a^(p/q) of a perfect q-th power folds to an exact
// value — 8^(2/3)=4, not just the 1/q roots. A non-exact root stays symbolic.
TEST_CASE("pow: perfect rational power a^(p/q)", "[1][pow][regression]") {
    REQUIRE(pow(integer(8), rational(2, 3)) == integer(4));
    REQUIRE(pow(integer(16), rational(3, 4)) == integer(8));
    REQUIRE(pow(integer(4), rational(3, 2)) == integer(8));
    REQUIRE(pow(integer(32), rational(2, 5)) == integer(4));
    REQUIRE(pow(rational(8, 27), rational(2, 3)) == rational(4, 9));
    REQUIRE(pow(integer(8), rational(-2, 3)) == rational(1, 4));
    // Unit-numerator roots unchanged; a prime base with p<q stays an
    // irreducible Pow (nothing pulls out).
    REQUIRE(pow(integer(27), rational(1, 3)) == integer(3));
    REQUIRE(pow(integer(2), rational(2, 3))->type_id() == TypeId::Pow);
}

// Non-unit rational power of a non-perfect base (POW-RAT-2): pull integer powers
// of each prime factor out of the radical, keeping primes separate to match
// SymPy's canonical form. 16^(2/3)=4·2^(2/3), 2^(5/2)=4√2, 12^(2/3)=2·2^(1/3)·3^(2/3).
TEST_CASE("pow: non-unit rational power extraction", "[1][pow][regression][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(pow(integer(16), rational(2, 3))->str(),
                              "4*2**Rational(2,3)"));
    REQUIRE(oracle.equivalent(pow(integer(2), rational(5, 2))->str(),
                              "4*sqrt(2)"));
    REQUIRE(oracle.equivalent(pow(integer(12), rational(2, 3))->str(),
                              "2*2**Rational(1,3)*3**Rational(2,3)"));
    REQUIRE(oracle.equivalent(pow(integer(72), rational(2, 3))->str(),
                              "12*3**Rational(1,3)"));
    REQUIRE(oracle.equivalent(pow(integer(7), rational(3, 2))->str(),
                              "7*sqrt(7)"));
    // No integer part pulls out → stays the irreducible Pow.
    REQUIRE(pow(integer(2), rational(2, 3))->type_id() == TypeId::Pow);
}

// Rational base, non-unit power (POW-RAT-2): (a/b)^(p/q) = (a·b^(q-1))^(p/q)/b^p
// rationalises the denominator. (2/3)^(2/3) = 2^(2/3)·3^(1/3)/3, (1/2)^(3/2) =
// √2/4, and a perfect rational power collapses fully ((8/27)^(2/3) = 4/9).
TEST_CASE("pow: non-unit rational power of a rational base",
          "[1][pow][regression][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(pow(rational(2, 3), rational(2, 3))->str(),
                              "2**Rational(2,3)*3**Rational(1,3)/3"));
    REQUIRE(oracle.equivalent(pow(rational(1, 2), rational(3, 2))->str(),
                              "sqrt(2)/4"));
    REQUIRE(oracle.equivalent(pow(rational(2, 9), rational(3, 2))->str(),
                              "2*sqrt(2)/27"));
    REQUIRE(oracle.equivalent(pow(rational(3, 2), rational(5, 2))->str(),
                              "9*sqrt(6)/8"));
    REQUIRE(pow(rational(8, 27), rational(2, 3)) == rational(4, 9));
}

// ----- Operator overloads -----------------------------------------------------

TEST_CASE("operator+ produces Add", "[1][operators]") {
    auto x = symbol("x");
    auto e = x + integer(1);
    REQUIRE(e->type_id() == TypeId::Add);
}

TEST_CASE("operator- = a + (-b)", "[1][operators]") {
    auto x = symbol("x");
    auto e = x - integer(1);
    // Phase 1f.1: this is Add(x, -1).
    REQUIRE(e->type_id() == TypeId::Add);
}

TEST_CASE("unary -x = -1*x", "[1][operators]") {
    auto x = symbol("x");
    auto e = -x;
    REQUIRE(e->type_id() == TypeId::Mul);
}

TEST_CASE("operator/ = a * b^-1", "[1][operators]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = x / y;
    // Phase 1f.1: x/y = Mul(x, Pow(y, -1)).
    REQUIRE(e->type_id() == TypeId::Mul);
}

// ----- Oracle-validated tests ------------------------------------------------
// These verify that whatever canonical form we produce is *equivalent* to the
// expected SymPy expression. SymPy will simplify(a - b) == 0 even when our
// representation hasn't done same-term collection yet.

TEST_CASE("Add / Mul / Pow expressions match SymPy", "[1][arithmetic][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");

    // Reference: test_arit0
    SECTION("a + b") {
        auto e = x + y;
        REQUIRE(oracle.equivalent(e->str(), "x + y"));
    }
    SECTION("a + a (canonical: 2*a after collection — equivalent now)") {
        auto e = x + x;
        REQUIRE(oracle.equivalent(e->str(), "2*x"));
    }
    SECTION("x*y + y*x") {
        auto e = mul(x, y) + mul(y, x);
        REQUIRE(oracle.equivalent(e->str(), "2*x*y"));
    }
    SECTION("(x + 1)^2 expanded, structural value equal") {
        auto e = pow(x + integer(1), integer(2));
        REQUIRE(oracle.equivalent(e->str(), "(x + 1)**2"));
    }
    SECTION("integer power") {
        auto e = pow(x, integer(3));
        REQUIRE(oracle.equivalent(e->str(), "x**3"));
    }
    SECTION("Rational coefficient") {
        auto e = mul(rational(1, 2), x);
        REQUIRE(oracle.equivalent(e->str(), "x/2"));
    }
    SECTION("subtraction") {
        auto e = x - integer(2);
        REQUIRE(oracle.equivalent(e->str(), "x - 2"));
    }
    SECTION("division") {
        auto e = x / integer(3);
        REQUIRE(oracle.equivalent(e->str(), "x/3"));
    }
    SECTION("polynomial-like sum") {
        auto e = pow(x, integer(2)) + integer(3) * x - integer(7);
        REQUIRE(oracle.equivalent(e->str(), "x**2 + 3*x - 7"));
    }
    SECTION("nested operations") {
        auto e = (x + integer(1)) * (x + integer(2));
        REQUIRE(oracle.equivalent(e->str(), "(x + 1)*(x + 2)"));
    }
}

TEST_CASE("Operator overloads match SymPy", "[1][operators][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");

    auto e1 = x + y - integer(2) * x + integer(3);
    REQUIRE(oracle.equivalent(e1->str(), "x + y - 2*x + 3"));

    auto e2 = x * y / pow(x, integer(2));
    REQUIRE(oracle.equivalent(e2->str(), "x*y/x**2"));

    auto e3 = -(x + y);
    REQUIRE(oracle.equivalent(e3->str(), "-(x + y)"));
}

TEST_CASE("Pure-number arithmetic results match SymPy", "[1][arithmetic][oracle]") {
    auto& oracle = Oracle::instance();

    // Reference: test_arit0 lines 67–105
    REQUIRE(oracle.equivalent(add(integer(5), integer(3))->str(), "8"));
    REQUIRE(oracle.equivalent(mul(integer(7), integer(6))->str(), "42"));
    REQUIRE(oracle.equivalent(pow(integer(2), integer(10))->str(), "1024"));
    REQUIRE(oracle.equivalent(add({integer(1), integer(2), integer(3), integer(4)})->str(), "10"));
    REQUIRE(oracle.equivalent(mul(rational(1, 2), rational(2, 3))->str(), "1/3"));
    REQUIRE(oracle.equivalent(pow(rational(2, 3), integer(3))->str(), "8/27"));
    REQUIRE(oracle.equivalent(pow(integer(2), integer(-5))->str(), "1/32"));
}
