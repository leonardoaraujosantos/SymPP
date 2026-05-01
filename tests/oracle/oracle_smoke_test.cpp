// Smoke tests for the Oracle harness itself.
// If these fail, nothing else can be trusted — the validation rig is broken.

#include <catch2/catch_test_macros.hpp>

#include "oracle/oracle.hpp"

using sympp::testing::Oracle;

TEST_CASE("oracle: ping returns SymPy version", "[0][oracle]") {
    auto& oracle = Oracle::instance();
    auto resp = oracle.send({{"op", "ping"}});

    REQUIRE(resp.ok);
    REQUIRE(resp.result_str() == "pong");
    REQUIRE_FALSE(resp.raw.value("sympy_version", "").empty());
    INFO("SymPy version reported by oracle: " << resp.raw.value("sympy_version", "?"));
}

TEST_CASE("oracle: srepr matches SymPy", "[0][oracle]") {
    // Reference: sympy/core/tests/test_symbol.py::test_Symbol
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.srepr("x") == "Symbol('x')");
    REQUIRE(oracle.srepr("y") == "Symbol('y')");
}

TEST_CASE("oracle: equivalent recognizes canonical identities", "[0][oracle]") {
    // From SymPy's docs and tests: sin(x)**2 + cos(x)**2 == 1
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent("sin(x)**2 + cos(x)**2", "1"));
    REQUIRE(oracle.equivalent("(x + 1)**2", "x**2 + 2*x + 1"));
    REQUIRE_FALSE(oracle.equivalent("x", "y"));
}

TEST_CASE("oracle: error path returns ok=false", "[0][oracle]") {
    auto& oracle = Oracle::instance();
    auto resp = oracle.send({{"op", "no_such_op"}});
    REQUIRE_FALSE(resp.ok);
    REQUIRE(resp.error() == "UnknownOp");
}

TEST_CASE("oracle: handles many sequential requests", "[0][oracle]") {
    // Ensures the long-lived subprocess + line-buffered I/O actually persists.
    auto& oracle = Oracle::instance();
    for (int i = 0; i < 50; ++i) {
        auto resp = oracle.send({{"op", "ping"}});
        REQUIRE(resp.ok);
    }
}
