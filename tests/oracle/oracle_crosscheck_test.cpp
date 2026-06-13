// Independent-oracle cross-checks.
//
// Every other oracle-backed test compares SymPP output against a *human-written*
// expected literal (SymPy only adjudicates the equivalence). That catches
// different-but-equal forms, but a wrong literal that happens to match a wrong
// SymPP answer slips through. Here SymPy computes the reference answer ITSELF
// via oracle.diff / oracle.integrate / oracle.simplify, and we require SymPP's
// independently-computed result to be equivalent. Neither side is hand-authored,
// so the two engines genuinely check each other.

#include <string>
#include <string_view>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/integrals/integrate.hpp>
#include <sympp/parsing/parser.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::parsing::parse;
using sympp::testing::Oracle;

// ORACLE-XCHECK-1 — d/dx computed by both engines must agree.
TEST_CASE("oracle cross-check: differentiation (ORACLE-XCHECK-1)",
          "[0][oracle][crosscheck]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    const std::vector<std::string_view> inputs = {
        "x**5 - 3*x**2 + 7",
        "sin(x)*cos(x)",
        "exp(x**2)",
        "log(x**2 + 1)",
        "x*exp(x)",
        "tan(x)",
        "atan(x)",
        "sinh(x)*cosh(x)",
        "x/(x**2 + 1)",
        "sqrt(x**2 + 1)",
        "x**x",
        "asin(x)",
    };
    for (auto s : inputs) {
        CAPTURE(s);
        const std::string sympp = diff(parse(s), x)->str();
        const std::string sympy = oracle.diff(s, "x");  // SymPy computes it
        REQUIRE(oracle.equivalent(sympp, sympy));
    }
}

// ORACLE-XCHECK-2 — antiderivatives computed by both engines must agree.
TEST_CASE("oracle cross-check: integration (ORACLE-XCHECK-2)",
          "[0][oracle][crosscheck]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    const std::vector<std::string_view> inputs = {
        "1/(x**2 + 1)",
        "x/(x**2 + 1)",
        "log(x)",
        "x*exp(x)",
        "1/(x*log(x))",
        "x**3 + 2*x",
        "cos(x)",
        "exp(x)*sin(x)",
    };
    for (auto s : inputs) {
        CAPTURE(s);
        const std::string sympp = integrate(parse(s), x)->str();
        const std::string sympy = oracle.integrate(s, "x");  // SymPy computes it
        // Indefinite integrals from these rules share the +C = 0 convention, so
        // the antiderivatives are equal, not merely equal up to a constant.
        REQUIRE(oracle.equivalent(sympp, sympy));
    }
}

// ORACLE-XCHECK-3 — simplified forms must be equivalent. simplify is heuristic,
// so the two engines may land on different shapes; equivalence is the contract.
TEST_CASE("oracle cross-check: simplification (ORACLE-XCHECK-3)",
          "[0][oracle][crosscheck]") {
    auto& oracle = Oracle::instance();
    const std::vector<std::string_view> inputs = {
        "sin(x)**2 + cos(x)**2",
        "(x**2 - 1)/(x - 1)",
        "exp(x)*exp(y)",
        "(x + 1)**2 - (x - 1)**2",
        "cos(x)**2 - sin(x)**2",
        "log(x) + log(y)",
    };
    for (auto s : inputs) {
        CAPTURE(s);
        const std::string sympp = simplify(parse(s))->str();
        const std::string sympy = oracle.simplify(s);  // SymPy computes it
        REQUIRE(oracle.equivalent(sympp, sympy));
    }
}
