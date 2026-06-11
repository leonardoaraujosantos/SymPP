// exp / log tests.
//
// Reference: sympy/functions/elementary/tests/test_exponential.py
//   * test_exp_values
//   * test_log_values
//   * test_log_assumptions

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/functions/exponential.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Auto-eval -------------------------------------------------------------

TEST_CASE("exp: exp(0) = 1", "[3c][exp]") {
    REQUIRE(exp(S::Zero()) == S::One());
}

TEST_CASE("exp: exp(1) = E", "[3c][exp]") {
    REQUIRE(exp(S::One()) == S::E());
}

// ----- Euler identity exp(r·I·π) (regression, EXP-1) -------------------------
// exp at imaginary multiples of π used to stay unevaluated. Now folds for a
// half-integer coefficient (q ∈ {1,2}) to ±1 / ±I, matching SymPy; π/3, π/4
// stay symbolic.
TEST_CASE("exp: Euler identity at imaginary π multiples", "[3c][exp][regression]") {
    auto pi = S::Pi();
    auto I = S::I();
    REQUIRE(exp(I * pi) == S::NegativeOne());            // e^{iπ} = −1
    REQUIRE(exp(integer(2) * I * pi) == S::One());       // e^{2iπ} = 1
    REQUIRE(exp(integer(3) * I * pi) == S::NegativeOne());
    REQUIRE(exp(I * pi / integer(2)) == I);              // e^{iπ/2} = I
    REQUIRE(exp(integer(5) * I * pi / integer(2)) == I); // periodic
}

TEST_CASE("exp: imaginary π/2 → −I via oracle", "[3c][exp][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto e = exp(mul(rational(-1, 2), mul(S::I(), S::Pi())));
    REQUIRE(oracle.equivalent(e->str(), "-I"));
}

TEST_CASE("exp: π/3 imaginary exponent stays symbolic", "[3c][exp][regression]") {
    // Denominator 3 → SymPy keeps exp(I*pi/3) unevaluated; so must SymPP.
    auto e = exp(S::I() * S::Pi() / integer(3));
    REQUIRE(e->type_id() == TypeId::Function);
}

TEST_CASE("log: log(1) = 0", "[3c][log]") {
    REQUIRE(log(S::One()) == S::Zero());
}

TEST_CASE("log: log(E) = 1", "[3c][log]") {
    REQUIRE(log(S::E()) == S::One());
}

// ----- log of negative / imaginary arguments (regression, LOG-1) -------------
// log(−x) = log(x) + Iπ and log(b·I) = log(|b|) + sign(b)·Iπ/2 used to stay
// unevaluated. The inverse counterpart to EXP-1; pairs with SQRT-3.
TEST_CASE("log: negative real → log|x| + Iπ", "[3c][log][oracle][regression]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(log(integer(-1))->str(), "I*pi"));
    REQUIRE(oracle.equivalent(log(integer(-2))->str(), "log(2) + I*pi"));
    REQUIRE(oracle.equivalent(log(mul(S::NegativeOne(), S::E()))->str(),
                              "1 + I*pi"));
}

TEST_CASE("log: imaginary axis → log|b| + sign·Iπ/2",
          "[3c][log][oracle][regression]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(log(S::I())->str(), "I*pi/2"));
    REQUIRE(oracle.equivalent(log(mul(S::NegativeOne(), S::I()))->str(),
                              "-I*pi/2"));
    REQUIRE(oracle.equivalent(log(mul(integer(2), S::I()))->str(),
                              "log(2) + I*pi/2"));
}

TEST_CASE("log: positive and symbolic args unaffected", "[3c][log][regression]") {
    auto x = symbol("x");
    REQUIRE(log(S::One()) == S::Zero());
    REQUIRE(log(integer(2))->type_id() == TypeId::Function);
    REQUIRE(log(x)->str() == "log(x)");
}

TEST_CASE("exp: stays unevaluated on a generic symbol", "[3c][exp]") {
    auto x = symbol("x");
    auto e = exp(x);
    REQUIRE(e->type_id() == TypeId::Function);
    REQUIRE(e->str() == "exp(x)");
}

TEST_CASE("log: stays unevaluated on a generic symbol", "[3c][log]") {
    auto x = symbol("x");
    auto e = log(x);
    REQUIRE(e->type_id() == TypeId::Function);
    REQUIRE(e->str() == "log(x)");
}

// ----- Inverse pair simplification ------------------------------------------

TEST_CASE("exp(log(x)) = x for positive x", "[3c][exp][log]") {
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    REQUIRE(exp(log(x)) == x);
}

TEST_CASE("exp(log(x)) stays unevaluated for unknown x", "[3c][exp][log]") {
    auto x = symbol("x");
    auto e = exp(log(x));
    REQUIRE(e->type_id() == TypeId::Function);
    REQUIRE(e->str() == "exp(log(x))");
}

TEST_CASE("log(exp(x)) = x for real x", "[3c][exp][log]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(log(exp(x)) == x);
}

// ----- Assumptions -----------------------------------------------------------

TEST_CASE("exp(real) is positive real", "[3c][exp][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true));
    REQUIRE(is_real(exp(x)) == true);
    REQUIRE(is_positive(exp(x)) == true);
    REQUIRE(is_nonzero(exp(x)) == true);
}

TEST_CASE("log(positive) is real", "[3c][log][assumptions]") {
    auto x = symbol("x", AssumptionMask{}.set_positive(true));
    REQUIRE(is_real(log(x)) == true);
}

// ----- Numeric evalf ---------------------------------------------------------

TEST_CASE("exp: numeric evalf via mpfr_exp", "[3c][exp][evalf]") {
    auto e = exp(float_value(1.0));
    REQUIRE(e->type_id() == TypeId::Float);
}

TEST_CASE("log: numeric evalf via mpfr_log", "[3c][log][evalf]") {
    auto e = log(float_value(2.0));
    REQUIRE(e->type_id() == TypeId::Float);
}

TEST_CASE("exp: high-precision matches SymPy", "[3c][exp][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = exp(float_value("1", 30));
    auto sympy_value = oracle.evalf("exp(Integer(1))", 30);
    INFO("SymPP: " << e->str());
    INFO("SymPy: " << sympy_value);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

TEST_CASE("log: log(2) at 50 digits matches SymPy", "[3c][log][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = log(float_value("2", 50));
    auto sympy_value = oracle.evalf("log(Integer(2))", 50);
    INFO("SymPP: " << e->str());
    INFO("SymPy: " << sympy_value);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

// ----- Substitution ----------------------------------------------------------

TEST_CASE("exp/log: subs replaces inside argument", "[3c][exp][log][subs]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = exp(x) + log(x);
    auto r = subs(e, x, y);
    REQUIRE(r == exp(y) + log(y));
}

TEST_CASE("exp(0) via subs from generic argument", "[3c][exp][subs]") {
    auto x = symbol("x");
    auto e = exp(x);
    REQUIRE(subs(e, x, S::Zero()) == S::One());
}

// ----- Oracle structural ----------------------------------------------------

TEST_CASE("exp/log: structural output matches SymPy",
          "[3c][exp][log][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(exp(x)->str(), "exp(x)"));
    REQUIRE(oracle.equivalent(log(x)->str(), "log(x)"));
    REQUIRE(oracle.equivalent((exp(x) + log(x))->str(), "exp(x) + log(x)"));
}
