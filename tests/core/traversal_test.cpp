// Traversal tests — subs, xreplace, has, free_symbols.
//
// References:
//   sympy/core/tests/test_subs.py::test_subs            (line 28)
//   sympy/core/tests/test_subs.py::test_dict_set        (xreplace map)
//   sympy/core/tests/test_basic.py::test_has, test_free_symbols
//   MATLAB Symbolic Math Toolbox §3 — `subs(expr, old, new)`

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/expr_collections.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- subs ------------------------------------------------------------------

TEST_CASE("subs: replaces atomic Symbol", "[1h][subs]") {
    // Reference: test_subs:  e = x;  e.subs(x, 3) == 3
    auto x = symbol("x");
    auto e = subs(x, x, integer(3));
    REQUIRE(e == integer(3));
}

TEST_CASE("subs: re-canonicalizes after substitution", "[1h][subs]") {
    // Reference: test_subs:  (2*x).subs(x, 3) == 6
    auto x = symbol("x");
    auto e = integer(2) * x;
    auto r = subs(e, x, integer(3));
    REQUIRE(r == integer(6));
}

TEST_CASE("subs: x in (x + y) -> a + y", "[1h][subs]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto a = symbol("a");
    auto e = x + y;
    auto r = subs(e, x, a);
    REQUIRE(r->str() == "a + y");
}

TEST_CASE("subs: substitutes inside Pow", "[1h][subs]") {
    auto x = symbol("x");
    auto r = subs(pow(x, integer(2)), x, integer(3));
    REQUIRE(r == integer(9));
}

TEST_CASE("subs: leaves expression unchanged when target absent", "[1h][subs]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    auto e = x + y;
    auto r = subs(e, z, integer(0));
    // No occurrences of z — result must be the same Expr (and same shared_ptr)
    REQUIRE(r == e);
    REQUIRE(r.get() == e.get());
}

TEST_CASE("subs: substitutes a non-atomic subexpression", "[1h][subs]") {
    // sin(x) doesn't exist yet; use a structural compound: (x+1) inside Pow.
    auto x = symbol("x");
    auto y = symbol("y");
    auto t = x + integer(1);  // canonical: x + 1
    auto e = pow(t, integer(2));
    auto r = subs(e, t, y);
    REQUIRE(r == pow(y, integer(2)));
    REQUIRE(r->str() == "y**2");
}

// ----- xreplace --------------------------------------------------------------

TEST_CASE("xreplace: simultaneous substitution", "[1h][xreplace]") {
    // Reference: test_dict_set — xreplace doesn't recurse into replaced subtrees,
    // so e.g. xreplace(x*y, {x: y, y: x}) gives y*x (not x*x via re-substitution).
    auto x = symbol("x");
    auto y = symbol("y");
    ExprMap<Expr> m;
    m.emplace(x, y);
    m.emplace(y, x);
    auto r = xreplace(x * y, m);
    // x*y -> y*x (same product after Mul canonicalization)
    REQUIRE(r->str() == "x*y");
}

TEST_CASE("xreplace: empty map is identity", "[1h][xreplace]") {
    auto x = symbol("x");
    auto e = x * x + integer(3);
    REQUIRE(xreplace(e, ExprMap<Expr>{}) == e);
}

TEST_CASE("xreplace: replaces multiple distinct symbols", "[1h][xreplace]") {
    auto x = symbol("x");
    auto y = symbol("y");
    ExprMap<Expr> m;
    m.emplace(x, integer(2));
    m.emplace(y, integer(3));
    auto r = xreplace(x * x + y, m);
    // 2*2 + 3 = 7
    REQUIRE(r == integer(7));
}

// ----- has -------------------------------------------------------------------

TEST_CASE("has: locates atomic targets", "[1h][has]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    auto e = x + y;

    REQUIRE(has(e, x));
    REQUIRE(has(e, y));
    REQUIRE_FALSE(has(e, z));
}

TEST_CASE("has: locates non-atomic subexpressions", "[1h][has]") {
    auto x = symbol("x");
    auto p = pow(x, integer(2));
    auto e = p + integer(1);

    REQUIRE(has(e, p));
    REQUIRE(has(e, x));        // x lives inside p
    REQUIRE(has(e, integer(2)));
    REQUIRE_FALSE(has(e, integer(3)));
}

TEST_CASE("has: empty when no occurrences", "[1h][has]") {
    auto x = symbol("x");
    REQUIRE_FALSE(has(integer(0), x));
}

// ----- free_symbols ----------------------------------------------------------

TEST_CASE("free_symbols: collects atomic Symbols", "[1h][free_symbols]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto fs = free_symbols(x + y);
    REQUIRE(fs.size() == 2);
    REQUIRE(fs.count(x) == 1);
    REQUIRE(fs.count(y) == 1);
}

TEST_CASE("free_symbols: pure-numeric expression has none", "[1h][free_symbols]") {
    auto e = integer(5) + rational(1, 2);
    auto fs = free_symbols(e);
    REQUIRE(fs.empty());
}

TEST_CASE("free_symbols: deduplicates", "[1h][free_symbols]") {
    auto x = symbol("x");
    auto fs = free_symbols(x * x + x);
    REQUIRE(fs.size() == 1);
    REQUIRE(fs.count(x) == 1);
}

TEST_CASE("free_symbols: descends into Pow and Mul", "[1h][free_symbols]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    auto e = pow(x, integer(2)) * (y + z);
    auto fs = free_symbols(e);
    REQUIRE(fs.size() == 3);
    REQUIRE(fs.count(x) == 1);
    REQUIRE(fs.count(y) == 1);
    REQUIRE(fs.count(z) == 1);
}

// ----- Oracle-validated tests ------------------------------------------------

TEST_CASE("subs: matches SymPy on polynomial substitution", "[1h][subs][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");

    // (x**2 + 3*x - 7).subs(x, 5) == 25 + 15 - 7 == 33
    auto e = pow(x, integer(2)) + integer(3) * x - integer(7);
    auto r = subs(e, x, integer(5));
    REQUIRE(r == integer(33));
    REQUIRE(oracle.equivalent(r->str(), "33"));
}

TEST_CASE("subs: matches SymPy on Rational substitution", "[1h][subs][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");

    // (x + x).subs(x, 1/2) -> 1
    auto r = subs(x + x, x, rational(1, 2));
    REQUIRE(r == integer(1));
    REQUIRE(oracle.equivalent(r->str(), "1"));
}

TEST_CASE("subs: leaves symbolic when target replaced by symbol",
          "[1h][subs][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");

    auto e = pow(x, integer(2)) + integer(2) * x + integer(1);
    auto r = subs(e, x, y);
    REQUIRE(oracle.equivalent(r->str(), "y**2 + 2*y + 1"));
}
