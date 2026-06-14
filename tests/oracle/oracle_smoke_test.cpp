// Smoke tests for the Oracle harness itself.
// If these fail, nothing else can be trusted — the validation rig is broken.

#include <string>
#include <utility>

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

// VERSION-GUARD-1 — pin the oracle to SymPy versions we have actually
// validated against. The whole suite trusts SymPy to adjudicate "equivalent",
// so a silent upgrade to an unvetted release could change that verdict (a new
// auto-simplification, a different canonical form) without any test noticing.
// Fail loudly here instead: an untested version forces a deliberate
// re-validation and an allowlist bump rather than letting drift slip through.
// CI installs 1.14; local development commonly runs 1.13.x — both are vetted.
namespace {
[[nodiscard]] std::pair<int, int> parse_major_minor(const std::string& v) {
    const auto d1 = v.find('.');
    if (d1 == std::string::npos) return {-1, -1};
    const auto d2 = v.find('.', d1 + 1);
    const int major = std::stoi(v.substr(0, d1));
    const int minor = std::stoi(v.substr(d1 + 1, d2 - d1 - 1));
    return {major, minor};
}
}  // namespace

TEST_CASE("oracle: SymPy version is one we validate against (VERSION-GUARD-1)",
          "[0][oracle]") {
    auto& oracle = Oracle::instance();
    const std::string version = oracle.sympy_version();
    INFO("oracle SymPy version: " << version);
    const auto [major, minor] = parse_major_minor(version);
    REQUIRE(major == 1);
    // Allowlist of validated minor releases. When CI/dev moves to a new SymPy,
    // re-run the suite against it and add its minor here in the same commit.
    const bool validated = (minor == 13 || minor == 14);
    if (!validated) {
        FAIL("Untested SymPy version "
             << version
             << " — re-validate the oracle suite against it, then add its "
                "major.minor to the allowlist in VERSION-GUARD-1.");
    }
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
