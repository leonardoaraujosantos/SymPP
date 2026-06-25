// Structured code-generation AST + emitter tests.

#include <string>

#include <catch2/catch_test_macros.hpp>

#include <sympp/codegen/codegen.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/printing/printing.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using namespace sympp::codegen;

namespace {

bool contains(const std::string& haystack, const std::string& needle) {
    return haystack.find(needle) != std::string::npos;
}

// f(x, y) = sin(x*y) + (x*y)^2 + (x*y) — x*y is the repeated subterm.
Expr sample_expr(const Expr& x, const Expr& y) {
    Expr xy = x * y;
    return sin(xy) + pow(xy, integer(2)) + xy;
}

}  // namespace

TEST_CASE("cse_codeblock: repeated x*y becomes one temporary", "[codegen][cse]") {
    auto x = symbol("x");
    auto y = symbol("y");

    CodeBlock block = cse_codeblock(sample_expr(x, y));

    // Exactly one shared subterm (x*y) → one temporary assignment.
    REQUIRE(block.body.size() == 1);
    REQUIRE(block.body[0].lhs->str() == "t0");
    REQUIRE(block.ret != nullptr);
}

TEST_CASE("emit_c: signature, CSE temp and return present", "[codegen][emit_c]") {
    auto x = symbol("x");
    auto y = symbol("y");

    FunctionDefinition fn{"f", {x, y}, cse_codeblock(sample_expr(x, y))};
    std::string code = emit_c(fn);

    REQUIRE(contains(code, "double f(double x, double y)"));
    REQUIRE(contains(code, "const double t0 ="));
    REQUIRE(contains(code, "return "));
    // The repeated factor is reified, so the rhs of t0 mentions both x and y.
    REQUIRE(contains(code, "x*y"));
}

TEST_CASE("emit_cxx: uses cxxcode for power", "[codegen][emit_cxx]") {
    auto x = symbol("x");
    auto y = symbol("y");

    FunctionDefinition fn{"g", {x, y}, cse_codeblock(sample_expr(x, y))};
    std::string code = emit_cxx(fn);

    REQUIRE(contains(code, "double g(double x, double y)"));
    REQUIRE(contains(code, "std::pow"));
    REQUIRE(contains(code, "return "));
}

TEST_CASE("emit_c: no-CSE expression has no temporaries", "[codegen][emit_c]") {
    auto x = symbol("x");

    FunctionDefinition fn{"h", {x}, cse_codeblock(x + integer(2))};
    std::string code = emit_c(fn);

    REQUIRE(fn.block.body.empty());
    REQUIRE(contains(code, "double h(double x)"));
    REQUIRE(contains(code, "return x + 2;"));
}

TEST_CASE("ccode of the t0 definition matches SymPy srepr-equivalent form",
          "[codegen][oracle]") {
    auto x = symbol("x");
    auto y = symbol("y");

    CodeBlock block = cse_codeblock(sample_expr(x, y));
    // The temporary's rhs is x*y; ccode renders it deterministically.
    REQUIRE(printing::ccode(block.body[0].rhs) == "x*y");

    // Cross-check the building block against SymPy via the oracle.
    auto& oracle = sympp::testing::Oracle::instance();
    REQUIRE(printing::srepr(x * y) == oracle.srepr("x*y"));
}
