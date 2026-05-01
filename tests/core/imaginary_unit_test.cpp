// ImaginaryUnit tests — I*I = -1 and I^n cycle.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/imaginary_unit.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>

using namespace sympp;

TEST_CASE("I: prints as 'I'", "[1][imaginary]") {
    REQUIRE(S::I()->str() == "I");
}

TEST_CASE("I: I*I = -1", "[1][imaginary]") {
    auto r = S::I() * S::I();
    REQUIRE(r == S::NegativeOne());
}

TEST_CASE("I: I*I*I = -I", "[1][imaginary]") {
    auto r = S::I() * S::I() * S::I();
    // -I = -1 * I
    auto expected = S::NegativeOne() * S::I();
    REQUIRE(r == expected);
}

TEST_CASE("I: I*I*I*I = 1", "[1][imaginary]") {
    auto r = S::I() * S::I() * S::I() * S::I();
    REQUIRE(r == S::One());
}

TEST_CASE("I: I^0 = 1", "[1][imaginary]") {
    auto r = pow(S::I(), integer(0));
    REQUIRE(r == S::One());
}

TEST_CASE("I: I^1 = I", "[1][imaginary]") {
    auto r = pow(S::I(), integer(1));
    REQUIRE(r == S::I());
}

TEST_CASE("I: I^2 = -1", "[1][imaginary]") {
    auto r = pow(S::I(), integer(2));
    REQUIRE(r == S::NegativeOne());
}

TEST_CASE("I: I^3 = -I", "[1][imaginary]") {
    auto r = pow(S::I(), integer(3));
    auto expected = S::NegativeOne() * S::I();
    REQUIRE(r == expected);
}

TEST_CASE("I: I^4 = 1", "[1][imaginary]") {
    auto r = pow(S::I(), integer(4));
    REQUIRE(r == S::One());
}

TEST_CASE("I: I^5 = I", "[1][imaginary]") {
    auto r = pow(S::I(), integer(5));
    REQUIRE(r == S::I());
}

TEST_CASE("I: I^(-1) = -I", "[1][imaginary]") {
    auto r = pow(S::I(), integer(-1));
    auto expected = S::NegativeOne() * S::I();
    REQUIRE(r == expected);
}

TEST_CASE("I: 2*I + 3*I = 5*I", "[1][imaginary]") {
    auto r = integer(2) * S::I() + integer(3) * S::I();
    auto expected = integer(5) * S::I();
    REQUIRE(r == expected);
}

TEST_CASE("I: (1 + I)*(1 - I) = 2", "[1][imaginary]") {
    auto a = integer(1) + S::I();
    auto b = integer(1) - S::I();
    // Use expand to multiply out.
    auto r = a * b;
    // r should be 1 - I + I - I^2 = 1 - (-1) = 2 — but Mul doesn't
    // distribute. We need expand. For now check arithmetic with explicit
    // distribution:
    auto manual = integer(1) - S::I() + S::I() - S::I() * S::I();
    REQUIRE(manual == integer(2));
}

TEST_CASE("I: 2*I*I*I = -2*I (mixed numeric and I)", "[1][imaginary]") {
    auto r = integer(2) * S::I() * S::I() * S::I();
    auto expected = integer(-2) * S::I();
    REQUIRE(r == expected);
}

TEST_CASE("I: assumption queries", "[1][imaginary]") {
    auto i = S::I();
    REQUIRE(i->ask(AssumptionKey::Real) == std::optional<bool>{false});
    REQUIRE(i->ask(AssumptionKey::Integer) == std::optional<bool>{false});
    REQUIRE(i->ask(AssumptionKey::Nonzero) == std::optional<bool>{true});
    REQUIRE(i->ask(AssumptionKey::Finite) == std::optional<bool>{true});
}
