// Infinity atoms and their arithmetic — +oo, -oo, zoo, nan (LIM-1).

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/add.hpp>
#include <sympp/core/infinity.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/parsing/parser.hpp>

using namespace sympp;

TEST_CASE("infinity: atoms print canonically", "[1][infinity]") {
    REQUIRE(S::Infinity()->str() == "oo");
    REQUIRE(S::NegativeInfinity()->str() == "-oo");
    REQUIRE(S::ComplexInfinity()->str() == "zoo");
    REQUIRE(S::NaN()->str() == "nan");
}

TEST_CASE("infinity: parser recognises oo / zoo / nan and -oo", "[1][infinity][parser]") {
    using parsing::parse;
    REQUIRE(parse("oo") == S::Infinity());
    REQUIRE(parse("zoo") == S::ComplexInfinity());
    REQUIRE(parse("nan") == S::NaN());
    REQUIRE(parse("-oo") == S::NegativeInfinity());
}

TEST_CASE("infinity: addition", "[1][infinity]") {
    auto oo = S::Infinity();
    REQUIRE(oo + integer(1) == oo);
    REQUIRE(oo - integer(1) == oo);
    REQUIRE(oo + oo == oo);
    REQUIRE(oo - oo == S::NaN());                         // oo - oo
    REQUIRE(S::NegativeInfinity() - oo == S::NegativeInfinity());
    REQUIRE(S::ComplexInfinity() + integer(1) == S::ComplexInfinity());
    REQUIRE(oo + S::NaN() == S::NaN());
}

TEST_CASE("infinity: multiplication", "[1][infinity]") {
    auto oo = S::Infinity();
    REQUIRE(oo * integer(2) == oo);
    REQUIRE(oo * integer(-2) == S::NegativeInfinity());
    REQUIRE(oo * oo == oo);
    REQUIRE(oo * integer(0) == S::NaN());                 // oo * 0
    REQUIRE(oo * S::NegativeInfinity() == S::NegativeInfinity());
    REQUIRE(S::NegativeInfinity() * S::NegativeInfinity() == oo);
}

TEST_CASE("infinity: powers", "[1][infinity]") {
    auto oo = S::Infinity();
    REQUIRE((integer(1) / oo) == S::Zero());              // 1/oo
    REQUIRE(pow(oo, integer(2)) == oo);
    REQUIRE(pow(oo, integer(-1)) == S::Zero());
    REQUIRE(pow(oo, integer(0)) == S::One());
    REQUIRE(pow(integer(2), oo) == oo);                   // 2^oo
    REQUIRE(pow(S::Half(), oo) == S::Zero());             // (1/2)^oo
    REQUIRE(pow(integer(1), oo) == S::NaN());             // 1^oo
    REQUIRE(pow(integer(0), oo) == S::Zero());            // 0^oo
    REQUIRE(pow(S::NegativeInfinity(), integer(2)) == oo);
    REQUIRE(pow(S::NegativeInfinity(), integer(3)) == S::NegativeInfinity());
}
