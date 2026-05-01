// Boolean / Relational / Piecewise tests.
//
// References:
//   sympy/core/tests/test_relational.py
//   sympy/functions/elementary/tests/test_piecewise.py

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/boolean.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/piecewise.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>

using namespace sympp;

// ----- Boolean atoms ---------------------------------------------------------

TEST_CASE("S::True / S::False are stable singletons", "[3k][bool]") {
    REQUIRE(S::True().get() == S::True().get());
    REQUIRE(S::False().get() == S::False().get());
    REQUIRE_FALSE(S::True() == S::False());
}

TEST_CASE("is_boolean_true / is_boolean_false", "[3k][bool]") {
    REQUIRE(is_boolean_true(S::True()));
    REQUIRE_FALSE(is_boolean_true(S::False()));
    REQUIRE(is_boolean_false(S::False()));
    REQUIRE_FALSE(is_boolean_false(S::True()));

    auto x = symbol("x");
    REQUIRE_FALSE(is_boolean_true(x));
    REQUIRE_FALSE(is_boolean_false(x));
}

// ----- Relational evaluation -------------------------------------------------

TEST_CASE("eq: numeric collapse", "[3k][rel]") {
    REQUIRE(eq(integer(3), integer(3)) == S::True());
    REQUIRE(eq(integer(3), integer(4)) == S::False());
    REQUIRE(eq(rational(1, 2), rational(2, 4)) == S::True());
}

TEST_CASE("lt/le/gt/ge: numeric collapse", "[3k][rel]") {
    REQUIRE(lt(integer(2), integer(5)) == S::True());
    REQUIRE(lt(integer(5), integer(2)) == S::False());
    REQUIRE(le(integer(5), integer(5)) == S::True());
    REQUIRE(gt(integer(7), integer(3)) == S::True());
    REQUIRE(ge(integer(3), integer(3)) == S::True());
}

TEST_CASE("eq: same expression collapses to True", "[3k][rel]") {
    auto x = symbol("x");
    REQUIRE(eq(x, x) == S::True());
    REQUIRE(ne(x, x) == S::False());
}

TEST_CASE("rel: stays unevaluated for symbolic args", "[3k][rel]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = lt(x, y);
    REQUIRE(e->type_id() == TypeId::Relational);
    REQUIRE(e->str() == "x < y");
}

// ----- Piecewise -------------------------------------------------------------

TEST_CASE("piecewise: first True branch wins", "[3k][piecewise]") {
    auto x = symbol("x");
    auto e = piecewise({{integer(1), S::True()}, {integer(2), S::True()}});
    // First branch's True condition wins immediately.
    REQUIRE(e == integer(1));
}

TEST_CASE("piecewise: skip False branches", "[3k][piecewise]") {
    auto e = piecewise({
        {integer(1), S::False()},
        {integer(2), S::False()},
        {integer(3), S::True()},
    });
    REQUIRE(e == integer(3));
}

TEST_CASE("piecewise: stays symbolic for unknown conditions", "[3k][piecewise]") {
    auto x = symbol("x");
    auto e = piecewise({
        {integer(1), lt(x, integer(0))},
        {integer(2), ge(x, integer(0))},
    });
    REQUIRE(e->type_id() == TypeId::Piecewise);
    REQUIRE(e->args().size() == 4);  // 2 branches × 2
}

TEST_CASE("piecewise: numeric Relational simplifies branches",
          "[3k][piecewise]") {
    // Build a Piecewise with the condition already-known-True after numeric
    // reduction.
    auto e = piecewise({
        {integer(10), lt(integer(5), integer(0))},   // False -> drop
        {integer(20), eq(integer(3), integer(3))},   // True  -> wins
        {integer(30), S::True()},                     // unreachable
    });
    REQUIRE(e == integer(20));
}

TEST_CASE("piecewise: subs drives numeric reduction",
          "[3k][piecewise][subs]") {
    auto x = symbol("x");
    // Build symbolic Piecewise, then subs.
    auto e = piecewise({
        {integer(-1), lt(x, integer(0))},
        {integer(0),  eq(x, integer(0))},
        {integer(1),  gt(x, integer(0))},
    });
    REQUIRE(subs(e, x, integer(-5)) == integer(-1));
    REQUIRE(subs(e, x, integer(0)) == integer(0));
    REQUIRE(subs(e, x, integer(7)) == integer(1));
}

TEST_CASE("piecewise: structural printing", "[3k][piecewise]") {
    auto x = symbol("x");
    auto e = piecewise({
        {integer(1), lt(x, integer(0))},
        {integer(2), S::True()},
    });
    REQUIRE(e->str() == "Piecewise((1, x < 0), (2, True))");
}
