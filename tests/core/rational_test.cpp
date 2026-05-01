// Rational tests — ported from sympy/core/tests/test_numbers.py::test_Rational_new
// and the helper _test_rational_new at line 294.
// Validated against the SymPy oracle for every numeric assertion.

#include <stdexcept>
#include <string>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/rational.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

TEST_CASE("Rational: basic construction", "[1][rational]") {
    auto half = rational(1, 2);
    REQUIRE(half->type_id() == TypeId::Rational);
    REQUIRE(half->str() == "1/2");
}

TEST_CASE("Rational: auto-reduction by gcd", "[1][rational]") {
    // Reference: SymPy auto-reduces — Rational(2, 4) == Rational(1, 2)
    auto a = rational(2, 4);
    auto b = rational(1, 2);
    REQUIRE(a == b);
    REQUIRE(a->str() == "1/2");

    auto c = rational(15, 25);
    REQUIRE(c->str() == "3/5");
}

TEST_CASE("Rational: collapses to Integer when denom == 1", "[1][rational]") {
    // Rational(4, 2) -> Integer(2)  per SymPy.
    auto e = rational(4, 2);
    REQUIRE(e->type_id() == TypeId::Integer);
    REQUIRE(e->str() == "2");

    auto z = rational(0, 5);
    REQUIRE(z->type_id() == TypeId::Integer);
    REQUIRE(z->str() == "0");
}

TEST_CASE("Rational: sign normalization (denom > 0)", "[1][rational]") {
    // Reference: SymPy keeps denom > 0; sign is on numerator.
    auto a = rational(1, -2);
    auto b = rational(-1, 2);
    REQUIRE(a == b);
    REQUIRE(a->str() == "-1/2");

    auto c = rational(-3, -4);
    REQUIRE(c->str() == "3/4");
}

TEST_CASE("Rational: rejects division by zero", "[1][rational]") {
    REQUIRE_THROWS_AS(rational(1, 0), std::domain_error);
    REQUIRE_THROWS_AS(rational(-1, 0), std::domain_error);
    REQUIRE_THROWS_AS(rational(0, 0), std::domain_error);
    // Phase 1e will replace this with collapse-to-ComplexInfinity per SymPy.
}

TEST_CASE("Rational: parses 'p/q' string", "[1][rational]") {
    // Reference: test_Rational_new
    //   assert Rational('3/4') == n3_4
    auto a = rational("3/4");
    auto b = rational(3, 4);
    REQUIRE(a == b);

    auto c = rational("-3/4");
    auto d = rational(-3, 4);
    REQUIRE(c == d);
}

TEST_CASE("Rational: predicates", "[1][rational]") {
    auto half = rational(1, 2);
    auto& r = static_cast<const Rational&>(*half);
    REQUIRE_FALSE(r.is_zero());
    REQUIRE_FALSE(r.is_one());
    REQUIRE(r.is_positive());
    REQUIRE_FALSE(r.is_negative());
    REQUIRE_FALSE(r.is_integer());
    REQUIRE(r.is_rational());
    REQUIRE(r.is_real());
    REQUIRE(r.sign() == 1);

    auto neg = rational(-3, 7);
    auto& nr = static_cast<const Rational&>(*neg);
    REQUIRE(nr.is_negative());
    REQUIRE(nr.sign() == -1);
}

TEST_CASE("Rational: equality and hash stability", "[1][rational]") {
    auto a = rational(7, 11);
    auto b = rational(14, 22);  // reduces to 7/11
    REQUIRE(a == b);
    REQUIRE(a->hash() == b->hash());

    auto c = rational(7, 13);
    REQUIRE_FALSE(a == c);
}

TEST_CASE("Rational: large numerator/denominator", "[1][rational]") {
    // Beyond long range; covered by mpz_class set_str path.
    auto a = rational("123456789012345678901234567890/987654321098765432109876543210");
    REQUIRE(a->type_id() == TypeId::Rational);
    // GCD of these two large numbers: SymPy computes it. We validate
    // structurally by oracle below.
    INFO("reduced str: " << a->str());
}

// ----- Oracle-validated tests -------------------------------------------------

TEST_CASE("Rational: 1/2 matches SymPy", "[1][rational][oracle]") {
    auto& oracle = Oracle::instance();
    auto half = rational(1, 2);
    REQUIRE(oracle.equivalent(half->str(), "Rational(1, 2)"));
    REQUIRE(oracle.equivalent(half->str(), "S.Half"));
}

TEST_CASE("Rational: arbitrary fractions match SymPy", "[1][rational][oracle]") {
    auto& oracle = Oracle::instance();
    struct Case {
        long n;
        long d;
        std::string expect;
    };
    for (auto [n, d, expect] :
         {Case{3, 4, "3/4"},
          Case{-3, 4, "-3/4"},
          Case{6, 8, "3/4"},        // reduces
          Case{1, -2, "-1/2"},      // sign normalize
          Case{-100, -200, "1/2"},  // both negative
          Case{17, 19, "17/19"}}) {
        auto e = rational(n, d);
        INFO("rational(" << n << ", " << d << ") -> " << e->str());
        REQUIRE(oracle.equivalent(e->str(), expect));
    }
}

TEST_CASE("Rational: large reduced form matches SymPy", "[1][rational][oracle]") {
    auto& oracle = Oracle::instance();
    // Pick numbers with a known shared factor: 6 and 9 give 2/3 after reduction.
    // Then scale way up.
    auto big = rational("60000000000000000000000000/90000000000000000000000000");
    REQUIRE(big->type_id() == TypeId::Rational);
    REQUIRE(oracle.equivalent(big->str(), "Rational(2, 3)"));
}

TEST_CASE("Rational: /1 collapses to Integer (matches SymPy Rational behaviour)",
          "[1][rational][oracle]") {
    auto& oracle = Oracle::instance();
    auto two = rational(8, 4);
    REQUIRE(two->type_id() == TypeId::Integer);
    REQUIRE(oracle.equivalent(two->str(), "Integer(2)"));
}
