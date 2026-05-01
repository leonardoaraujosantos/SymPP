// Integer tests — ported and adapted from
//   sympy/core/tests/test_numbers.py::test_Integer_new
//   sympy/core/tests/test_numbers.py::test_IntegerInteger
// Validated against the SymPy oracle for every numeric assertion.

#include <climits>
#include <stdexcept>
#include <string>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

TEST_CASE("Integer: construction from int", "[1][integer]") {
    auto z = integer(0);
    auto one = integer(1);
    auto neg = integer(-7);
    auto big = integer(1234567890);

    REQUIRE(z->str() == "0");
    REQUIRE(one->str() == "1");
    REQUIRE(neg->str() == "-7");
    REQUIRE(big->str() == "1234567890");
}

TEST_CASE("Integer: construction from long long beyond 32-bit range", "[1][integer]") {
    long long huge = 9'223'372'036'854'775'807LL;  // INT64_MAX
    auto z = integer(huge);
    REQUIRE(z->str() == "9223372036854775807");
}

TEST_CASE("Integer: construction from arbitrary-precision string", "[1][integer]") {
    // Beyond any builtin integer width — must round-trip exactly.
    auto z = integer("123456789012345678901234567890");
    REQUIRE(z->str() == "123456789012345678901234567890");
}

TEST_CASE("Integer: rejects non-integer strings", "[1][integer]") {
    // Reference: sympy/core/tests/test_numbers.py::test_Integer_new
    //   raises(ValueError, lambda: Integer("10.5"))
    REQUIRE_THROWS_AS(integer("10.5"), std::invalid_argument);
    REQUIRE_THROWS_AS(integer("1e3"), std::invalid_argument);
    REQUIRE_THROWS_AS(integer("not a number"), std::invalid_argument);
}

TEST_CASE("Integer: equality and hash", "[1][integer]") {
    // Reference: sympy/core/tests/test_numbers.py::test_IntegerInteger
    auto a = integer(42);
    auto b = integer(42);
    auto c = integer(43);

    REQUIRE(a == b);
    REQUIRE_FALSE(a == c);
    REQUIRE(a->hash() == b->hash());
    // Different values *very probably* hash differently. Not guaranteed but
    // a regression here would be suspicious.
    REQUIRE(a->hash() != c->hash());
}

TEST_CASE("Integer: predicates", "[1][integer]") {
    auto z = integer(0);
    auto p = integer(7);
    auto n = integer(-3);

    auto& zr = static_cast<const Integer&>(*z);
    auto& pr = static_cast<const Integer&>(*p);
    auto& nr = static_cast<const Integer&>(*n);

    REQUIRE(zr.is_zero());
    REQUIRE_FALSE(zr.is_one());
    REQUIRE_FALSE(zr.is_negative());
    REQUIRE_FALSE(zr.is_positive());

    REQUIRE_FALSE(pr.is_zero());
    REQUIRE(pr.is_positive());
    REQUIRE_FALSE(pr.is_negative());

    REQUIRE_FALSE(nr.is_zero());
    REQUIRE(nr.is_negative());
    REQUIRE_FALSE(nr.is_positive());

    REQUIRE(zr.sign() == 0);
    REQUIRE(pr.sign() == 1);
    REQUIRE(nr.sign() == -1);

    REQUIRE(pr.is_integer());
    REQUIRE(pr.is_rational());
    REQUIRE(pr.is_real());
}

TEST_CASE("Integer: type tag", "[1][integer]") {
    auto z = integer(5);
    REQUIRE(z->type_id() == TypeId::Integer);
    REQUIRE(z->is_atomic());
    REQUIRE(z->args().empty());
}

TEST_CASE("Integer: fits_long / to_long", "[1][integer]") {
    auto small = integer(123);
    auto& sr = static_cast<const Integer&>(*small);
    REQUIRE(sr.fits_long());
    REQUIRE(sr.to_long() == 123L);

    auto huge = integer("123456789012345678901234567890");
    auto& hr = static_cast<const Integer&>(*huge);
    REQUIRE_FALSE(hr.fits_long());
}

// ----- Oracle-validated tests -------------------------------------------------

TEST_CASE("Integer: str() round-trips through SymPy", "[1][integer][oracle]") {
    auto& oracle = Oracle::instance();
    for (int v : {0, 1, -1, 42, -42, 1'000'000}) {
        auto z = integer(v);
        REQUIRE(oracle.equivalent(z->str(), std::to_string(v)));
        REQUIRE(oracle.srepr(z->str()) == "Integer(" + std::to_string(v) + ")");
    }
}

TEST_CASE("Integer: large value matches SymPy", "[1][integer][oracle]") {
    auto& oracle = Oracle::instance();
    // Reference: sympy/core/tests/test_numbers.py — Python int arbitrary size
    auto z = integer("123456789012345678901234567890");
    REQUIRE(oracle.equivalent(z->str(), "123456789012345678901234567890"));
}

TEST_CASE("Integer: negative large value matches SymPy", "[1][integer][oracle]") {
    auto& oracle = Oracle::instance();
    auto z = integer("-99999999999999999999999999999999");
    REQUIRE(oracle.equivalent(z->str(), "-99999999999999999999999999999999"));
}
