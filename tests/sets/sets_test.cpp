// Set algebra tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/sets/sets.hpp>

using namespace sympp;

// ----- EmptySet --------------------------------------------------------------

TEST_CASE("set: empty_set contains nothing", "[10][sets][empty]") {
    auto s = empty_set();
    REQUIRE(s->contains(integer(0)) == std::optional<bool>{false});
    REQUIRE(s->contains(integer(42)) == std::optional<bool>{false});
    REQUIRE(s->str() == "EmptySet");
}

// ----- Reals -----------------------------------------------------------------

TEST_CASE("set: reals contains integers and rationals",
          "[10][sets][reals]") {
    auto s = reals();
    REQUIRE(s->contains(integer(5)) == std::optional<bool>{true});
    REQUIRE(s->contains(rational(1, 3)) == std::optional<bool>{true});
}

// ----- Interval --------------------------------------------------------------

TEST_CASE("set: closed interval [0, 10] contains 5",
          "[10][sets][interval]") {
    auto s = interval(integer(0), integer(10));
    REQUIRE(s->contains(integer(5)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(0)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(10)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(11)) == std::optional<bool>{false});
    REQUIRE(s->contains(integer(-1)) == std::optional<bool>{false});
}

TEST_CASE("set: open interval (0, 10) excludes endpoints",
          "[10][sets][interval]") {
    auto s = interval(integer(0), integer(10), true, true);
    REQUIRE(s->contains(integer(0)) == std::optional<bool>{false});
    REQUIRE(s->contains(integer(10)) == std::optional<bool>{false});
    REQUIRE(s->contains(integer(5)) == std::optional<bool>{true});
}

TEST_CASE("set: half-open [0, 10)", "[10][sets][interval]") {
    auto s = interval(integer(0), integer(10), false, true);
    REQUIRE(s->contains(integer(0)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(10)) == std::optional<bool>{false});
}

TEST_CASE("set: interval prints correctly", "[10][sets][interval]") {
    auto a = interval(integer(0), integer(10));
    REQUIRE(a->str() == "[0, 10]");
    auto b = interval(integer(0), integer(10), true, true);
    REQUIRE(b->str() == "(0, 10)");
    auto c = interval(integer(0), integer(10), false, true);
    REQUIRE(c->str() == "[0, 10)");
}

// ----- FiniteSet -------------------------------------------------------------

TEST_CASE("set: finite_set contains its elements",
          "[10][sets][finite]") {
    auto s = finite_set({integer(1), integer(2), integer(3)});
    REQUIRE(s->contains(integer(1)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(2)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(3)) == std::optional<bool>{true});
}

TEST_CASE("set: empty finite_set collapses to EmptySet",
          "[10][sets][finite]") {
    auto s = finite_set({});
    REQUIRE(s->kind() == SetKind::Empty);
}

TEST_CASE("set: finite_set prints with curly braces",
          "[10][sets][finite]") {
    auto s = finite_set({integer(1), integer(2), integer(3)});
    REQUIRE(s->str() == "{1, 2, 3}");
}

// ----- Union -----------------------------------------------------------------

TEST_CASE("set: union contains members of either",
          "[10][sets][union]") {
    auto a = finite_set({integer(1), integer(2)});
    auto b = finite_set({integer(3), integer(4)});
    auto u = set_union(a, b);
    REQUIRE(u->contains(integer(2)) == std::optional<bool>{true});
    REQUIRE(u->contains(integer(3)) == std::optional<bool>{true});
}

TEST_CASE("set: union with empty returns the other",
          "[10][sets][union]") {
    auto a = finite_set({integer(1), integer(2)});
    auto u = set_union(a, empty_set());
    REQUIRE(u->kind() == SetKind::FiniteSet);
}

// ----- Intersection ----------------------------------------------------------

TEST_CASE("set: intersection requires both",
          "[10][sets][intersection]") {
    auto a = finite_set({integer(1), integer(2), integer(3)});
    auto b = interval(integer(2), integer(10));
    auto i = set_intersection(a, b);
    REQUIRE(i->contains(integer(2)) == std::optional<bool>{true});
    REQUIRE(i->contains(integer(3)) == std::optional<bool>{true});
    // 1 ∈ a but 1 ∉ [2, 10] → 1 not in intersection
    REQUIRE(i->contains(integer(1)) == std::optional<bool>{false});
}

TEST_CASE("set: intersection with empty is empty",
          "[10][sets][intersection]") {
    auto a = finite_set({integer(1)});
    auto i = set_intersection(a, empty_set());
    REQUIRE(i->kind() == SetKind::Empty);
}

// ----- Complement ------------------------------------------------------------

TEST_CASE("set: complement removes elements",
          "[10][sets][complement]") {
    auto u = finite_set({integer(1), integer(2), integer(3)});
    auto r = finite_set({integer(2)});
    auto c = set_complement(u, r);
    REQUIRE(c->contains(integer(1)) == std::optional<bool>{true});
    REQUIRE(c->contains(integer(2)) == std::optional<bool>{false});
    REQUIRE(c->contains(integer(3)) == std::optional<bool>{true});
}

TEST_CASE("set: complement with empty returns universal",
          "[10][sets][complement]") {
    auto u = finite_set({integer(1), integer(2)});
    auto c = set_complement(u, empty_set());
    REQUIRE(c->kind() == SetKind::FiniteSet);
}
