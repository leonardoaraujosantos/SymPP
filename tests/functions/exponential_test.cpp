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
#include <sympp/core/operators.hpp>
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

TEST_CASE("log: log(1) = 0", "[3c][log]") {
    REQUIRE(log(S::One()) == S::Zero());
}

TEST_CASE("log: log(E) = 1", "[3c][log]") {
    REQUIRE(log(S::E()) == S::One());
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
