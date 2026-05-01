// Float / evalf tests — patterns drawn from
//   sympy/core/tests/test_numbers.py (Float construction, comparisons)
//   sympy/core/numbers.py docstring examples for Float and evalf
// Numeric assertions cross-checked against the SymPy oracle.

#include <stdexcept>
#include <string>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/symbol.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

TEST_CASE("Float: construction from double", "[1][float]") {
    auto f = float_value(1.5);
    REQUIRE(f->type_id() == TypeId::Float);

    auto z = float_value(0.0);
    auto& zr = static_cast<const Float&>(*z);
    REQUIRE(zr.is_zero());
}

TEST_CASE("Float: construction from string preserves precision", "[1][float]") {
    // Reference: sympy Float('1.23456789012345', 15) — exact at default precision
    auto f = float_value("1.23456789012345");
    REQUIRE(f->type_id() == TypeId::Float);
    auto& fr = static_cast<const Float&>(*f);
    REQUIRE(fr.precision_dps() == kDefaultDps);
}

TEST_CASE("Float: variable precision (dps)", "[1][float]") {
    // SymPy: Float('3.14', 30).evalf() -> 30 digits of pi
    auto f = float_value("3.14159265358979323846264338327950288", 50);
    auto& fr = static_cast<const Float&>(*f);
    REQUIRE(fr.precision_dps() >= 50);
}

TEST_CASE("Float: predicates", "[1][float]") {
    auto z = float_value(0.0);
    auto p = float_value(1.5);
    auto n = float_value(-2.7);

    auto& zr = static_cast<const Float&>(*z);
    auto& pr = static_cast<const Float&>(*p);
    auto& nr = static_cast<const Float&>(*n);

    REQUIRE(zr.is_zero());
    REQUIRE(pr.is_positive());
    REQUIRE(nr.is_negative());

    REQUIRE(zr.sign() == 0);
    REQUIRE(pr.sign() == 1);
    REQUIRE(nr.sign() == -1);

    // Float is real but not rational (per SymPy convention).
    REQUIRE(pr.is_real());
    REQUIRE_FALSE(pr.is_rational());
}

TEST_CASE("Float: rejects invalid string", "[1][float]") {
    REQUIRE_THROWS_AS(float_value("not a number"), std::invalid_argument);
}

TEST_CASE("Float: hash stable for equal values", "[1][float]") {
    auto a = float_value(2.5);
    auto b = float_value(2.5);
    REQUIRE(a == b);
    REQUIRE(a->hash() == b->hash());
}

// ----- evalf -------------------------------------------------------------------

TEST_CASE("evalf: Integer -> Float", "[1][evalf]") {
    auto z = integer(7);
    auto f = evalf(z, 15);
    REQUIRE(f->type_id() == TypeId::Float);
    REQUIRE(static_cast<const Float&>(*f).is_one() == false);
    REQUIRE(static_cast<const Float&>(*f).sign() == 1);
}

TEST_CASE("evalf: Rational -> Float", "[1][evalf]") {
    // Reference: SymPy Rational(1, 3).evalf() ≈ 0.333333333333333
    auto q = rational(1, 3);
    auto f = evalf(q, 15);
    REQUIRE(f->type_id() == TypeId::Float);
}

TEST_CASE("evalf: Float -> Float at new precision", "[1][evalf]") {
    auto f1 = float_value("3.14", 15);
    auto f2 = evalf(f1, 50);
    REQUIRE(f2->type_id() == TypeId::Float);
    REQUIRE(static_cast<const Float&>(*f2).precision_dps() >= 50);
}

TEST_CASE("evalf: Symbol passes through unchanged", "[1][evalf]") {
    // Phase 1f will add proper handling for compound expressions; for now
    // non-numeric atoms are returned as-is.
    auto x = make<Symbol>(std::string{"x"});
    Expr xe = x;
    auto out = evalf(xe, 15);
    REQUIRE(out == xe);
}

// ----- Oracle-validated tests --------------------------------------------------

TEST_CASE("Float: 1 over 3 at 30 digits matches SymPy", "[1][float][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto q = rational(1, 3);
    auto f = evalf(q, 30);
    // Build the SymPy expression "Rational(1,3).evalf(30)" via oracle helper.
    auto sympy_value = oracle.evalf("Rational(1, 3)", 30);
    INFO("SymPP:  " << f->str());
    INFO("SymPy:  " << sympy_value);
    // Compare numerically: equivalent simplifies a-b == 0 — works on Floats too.
    REQUIRE(oracle.equivalent(f->str(), sympy_value));
}

TEST_CASE("Float: 1/7 at 50 digits matches SymPy", "[1][float][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto q = rational(1, 7);
    auto f = evalf(q, 50);
    auto sympy_value = oracle.evalf("Rational(1, 7)", 50);
    INFO("SymPP:  " << f->str());
    INFO("SymPy:  " << sympy_value);
    REQUIRE(oracle.equivalent(f->str(), sympy_value));
}

TEST_CASE("Float: small integer evalf matches SymPy", "[1][float][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    for (int v : {0, 1, -1, 7, -42, 100}) {
        auto z = integer(v);
        auto f = evalf(z, 20);
        auto sympy_value = oracle.evalf(std::to_string(v), 20);
        INFO("v=" << v << "  SymPP=" << f->str() << "  SymPy=" << sympy_value);
        REQUIRE(oracle.equivalent(f->str(), sympy_value));
    }
}

TEST_CASE("Float: large integer evalf matches SymPy", "[1][float][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto z = integer("123456789012345678901234567890");
    auto f = evalf(z, 40);
    auto sympy_value = oracle.evalf("123456789012345678901234567890", 40);
    INFO("SymPP:  " << f->str());
    INFO("SymPy:  " << sympy_value);
    REQUIRE(oracle.equivalent(f->str(), sympy_value));
}
