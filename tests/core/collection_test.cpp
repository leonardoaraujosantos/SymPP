// Term collection (Add) and base collection (Mul) tests.
// References:
//   sympy/core/tests/test_arit.py::test_arit0   (a + a == 2*a, etc.)
//   sympy/core/tests/test_arit.py::test_pow_arit
// These tests check the *canonical form* directly — that our str() output
// already shows the collection, not just that it's value-equivalent.

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

// ----- Add term collection ---------------------------------------------------

TEST_CASE("Add: x + x = 2*x", "[1f2][add][collection]") {
    auto x = symbol("x");
    auto e = x + x;
    REQUIRE(e->type_id() == TypeId::Mul);  // 2*x is a Mul
    REQUIRE(e->str() == "2*x");
}

TEST_CASE("Add: x - x = 0", "[1f2][add][collection]") {
    auto x = symbol("x");
    auto e = x - x;
    REQUIRE(e == S::Zero());
}

TEST_CASE("Add: 2*x + 3*x = 5*x", "[1f2][add][collection]") {
    auto x = symbol("x");
    auto e = integer(2) * x + integer(3) * x;
    REQUIRE(e->str() == "5*x");
}

TEST_CASE("Add: 2*x + 3*x - 5*x = 0", "[1f2][add][collection]") {
    auto x = symbol("x");
    auto e = integer(2) * x + integer(3) * x - integer(5) * x;
    REQUIRE(e == S::Zero());
}

TEST_CASE("Add: a + b + a + b - 2*a = 2*b", "[1f2][add][collection]") {
    auto a = symbol("a");
    auto b = symbol("b");
    auto e = a + b + a + b - integer(2) * a;
    REQUIRE(e->str() == "2*b");
}

TEST_CASE("Add: x*y + y*x = 2*x*y after Mul canonicalization",
          "[1f2][add][collection]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = mul(x, y) + mul(y, x);  // both Muls canonicalize to the same args
    REQUIRE(e->type_id() == TypeId::Mul);
    REQUIRE(e->str() == "2*x*y");
}

TEST_CASE("Add: x + x**2 + x stays separate (different terms)",
          "[1f2][add][collection]") {
    auto x = symbol("x");
    auto e = x + pow(x, integer(2)) + x;
    // Expect: 2*x + x**2  (two args after collection)
    REQUIRE(e->type_id() == TypeId::Add);
    REQUIRE(e->args().size() == 2);
}

// ----- Mul base collection ---------------------------------------------------

TEST_CASE("Mul: x * x = x**2", "[1f2][mul][collection]") {
    auto x = symbol("x");
    auto e = x * x;
    REQUIRE(e->type_id() == TypeId::Pow);
    REQUIRE(e->str() == "x**2");
}

TEST_CASE("Mul: x**2 * x**3 = x**5", "[1f2][mul][collection]") {
    auto x = symbol("x");
    auto e = pow(x, integer(2)) * pow(x, integer(3));
    REQUIRE(e->str() == "x**5");
}

TEST_CASE("Mul: x * x**(-1) = 1", "[1f2][mul][collection]") {
    auto x = symbol("x");
    auto e = x * pow(x, integer(-1));
    REQUIRE(e == S::One());
}

TEST_CASE("Mul: x / x = 1", "[1f2][mul][collection]") {
    auto x = symbol("x");
    auto e = x / x;
    REQUIRE(e == S::One());
}

TEST_CASE("Mul: x**2 / x = x", "[1f2][mul][collection]") {
    auto x = symbol("x");
    auto e = pow(x, integer(2)) / x;
    REQUIRE(e == x);
}

TEST_CASE("Mul: 2*x * 3*x = 6*x**2", "[1f2][mul][collection]") {
    auto x = symbol("x");
    auto e = (integer(2) * x) * (integer(3) * x);
    REQUIRE(e->str() == "6*x**2");
}

TEST_CASE("Mul: x*y*x*y = x**2 * y**2", "[1f2][mul][collection]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = x * y * x * y;
    // Expect: x**2 * y**2 — two Pow nodes after collection.
    REQUIRE(e->type_id() == TypeId::Mul);
    REQUIRE(e->args().size() == 2);
}

// ----- Combined Add + Mul collection ----------------------------------------

TEST_CASE("Combined: (x + 1) - 1 = x", "[1f2][collection]") {
    auto x = symbol("x");
    auto e = (x + integer(1)) - integer(1);
    REQUIRE(e == x);
}

TEST_CASE("Combined: x * 2 / 2 = x", "[1f2][collection]") {
    auto x = symbol("x");
    auto e = x * integer(2) / integer(2);
    REQUIRE(e == x);
}

TEST_CASE("Combined: x + x*y - x = x*y", "[1f2][collection]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = x + x * y - x;
    REQUIRE(e->str() == "x*y");
}

// ----- Oracle re-validation -------------------------------------------------
// Now our canonical form should match SymPy's str representation more closely.
// Test that our str() output AND a SymPy-equivalent reference both reduce to
// the same collected form.

TEST_CASE("Collection: oracle agrees on canonical forms", "[1f2][collection][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");

    // Reference: test_arit0
    REQUIRE(oracle.equivalent((x + x)->str(), "2*x"));
    REQUIRE(oracle.equivalent((x - x)->str(), "0"));
    REQUIRE(oracle.equivalent((integer(2) * x + integer(3) * x)->str(), "5*x"));
    REQUIRE(oracle.equivalent((x * x)->str(), "x**2"));
    REQUIRE(oracle.equivalent((x * y * x)->str(), "x**2 * y"));
    REQUIRE(oracle.equivalent((pow(x, integer(2)) * pow(x, integer(3)))->str(), "x**5"));
    REQUIRE(oracle.equivalent(((x + integer(1)) * (x + integer(2)))->str(),
                              "(x + 1)*(x + 2)"));
}

TEST_CASE("Collection: deeper algebraic identities", "[1f2][collection][oracle]") {
    auto& oracle = Oracle::instance();
    auto a = symbol("a");
    auto b = symbol("b");

    // a*b + b*a + a*b + 5*b*a should collect to 8*a*b
    auto e1 = mul(a, b) + mul(b, a) + mul(a, b) + integer(5) * mul(b, a);
    REQUIRE(oracle.equivalent(e1->str(), "8*a*b"));

    // (a + b)^2 stays unevaluated (no expand yet — Phase 1i)
    auto e2 = pow(a + b, integer(2));
    REQUIRE(oracle.equivalent(e2->str(), "(a + b)**2"));

    // Subtraction and re-collection
    auto e3 = integer(2) * a + integer(3) * b - a - b;
    REQUIRE(oracle.equivalent(e3->str(), "a + 2*b"));
}
