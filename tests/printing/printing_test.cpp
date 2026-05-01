// Phase 13 — printing/codegen tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/printing/printing.hpp>

using namespace sympp;
using namespace sympp::printing;

// ----- ccode -----------------------------------------------------------------

TEST_CASE("ccode: x + 2", "[13][ccode]") {
    auto x = symbol("x");
    REQUIRE(ccode(x + integer(2)) == "x + 2");
}

TEST_CASE("ccode: x^2 → pow(x, 2)", "[13][ccode]") {
    auto x = symbol("x");
    REQUIRE(ccode(pow(x, integer(2))) == "pow(x, 2)");
}

TEST_CASE("ccode: sin(x) → sin(x)", "[13][ccode]") {
    auto x = symbol("x");
    REQUIRE(ccode(sin(x)) == "sin(x)");
}

TEST_CASE("ccode: pi → M_PI", "[13][ccode]") {
    REQUIRE(ccode(S::Pi()) == "M_PI");
}

TEST_CASE("ccode: rational forces double division",
          "[13][ccode]") {
    auto x = symbol("x");
    auto e = rational(1, 2) * x;
    auto out = ccode(e);
    // Some form of "1.0/2.0" division should appear.
    REQUIRE(out.find("1.0/2.0") != std::string::npos);
}

TEST_CASE("ccode: x/y emits division", "[13][ccode]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto out = ccode(x / y);
    REQUIRE(out.find("/") != std::string::npos);
    REQUIRE(out.find("x") != std::string::npos);
    REQUIRE(out.find("y") != std::string::npos);
}

// ----- cxxcode --------------------------------------------------------------

TEST_CASE("cxxcode: sin(x) → std::sin(x)", "[13][cxxcode]") {
    auto x = symbol("x");
    REQUIRE(cxxcode(sin(x)) == "std::sin(x)");
}

TEST_CASE("cxxcode: pow(x, 3) → std::pow", "[13][cxxcode]") {
    auto x = symbol("x");
    REQUIRE(cxxcode(pow(x, integer(3))) == "std::pow(x, 3)");
}

TEST_CASE("cxxcode: pi uses std::numbers::pi_v",
          "[13][cxxcode]") {
    auto out = cxxcode(S::Pi());
    REQUIRE(out.find("std::numbers::pi_v") != std::string::npos);
}

// ----- fcode ----------------------------------------------------------------

TEST_CASE("fcode: x^2 uses **", "[13][fcode]") {
    auto x = symbol("x");
    REQUIRE(fcode(pow(x, integer(2))) == "x**2.0d0");
}

TEST_CASE("fcode: integer literal has d0 suffix", "[13][fcode]") {
    REQUIRE(fcode(integer(42)) == "42.0d0");
}

// ----- latex ----------------------------------------------------------------

TEST_CASE("latex: x + 2 plain text", "[13][latex]") {
    auto x = symbol("x");
    REQUIRE(latex(x + integer(2)) == "x + 2");
}

TEST_CASE("latex: 1/2 → \\frac{1}{2}", "[13][latex]") {
    REQUIRE(latex(rational(1, 2)) == "\\frac{1}{2}");
}

TEST_CASE("latex: sin(x) → \\sin(x)", "[13][latex]") {
    auto x = symbol("x");
    REQUIRE(latex(sin(x)) == "\\sin(x)");
}

TEST_CASE("latex: sqrt(x) → \\sqrt{x}", "[13][latex]") {
    auto x = symbol("x");
    REQUIRE(latex(sqrt(x)) == "\\sqrt{x}");
}

TEST_CASE("latex: greek symbol alpha", "[13][latex]") {
    auto a = symbol("alpha");
    REQUIRE(latex(a) == "\\alpha");
}

TEST_CASE("latex: x^2", "[13][latex]") {
    auto x = symbol("x");
    REQUIRE(latex(pow(x, integer(2))) == "x^{2}");
}

TEST_CASE("latex: pi → \\pi", "[13][latex]") {
    REQUIRE(latex(S::Pi()) == "\\pi");
}

// ----- octave_code ----------------------------------------------------------

TEST_CASE("octave: x*y uses .*", "[13][octave]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto out = octave_code(x * y);
    REQUIRE(out.find(".*") != std::string::npos);
}

TEST_CASE("octave: x^3 uses .^", "[13][octave]") {
    auto x = symbol("x");
    REQUIRE(octave_code(pow(x, integer(3))) == "x.^3");
}

TEST_CASE("octave: sin(x)", "[13][octave]") {
    auto x = symbol("x");
    REQUIRE(octave_code(sin(x)) == "sin(x)");
}

// ----- pretty ---------------------------------------------------------------

TEST_CASE("pretty: x + 2", "[13][pretty]") {
    auto x = symbol("x");
    REQUIRE(pretty(x + integer(2)) == "x + 2");
}

TEST_CASE("pretty: x^2 uses **", "[13][pretty]") {
    auto x = symbol("x");
    REQUIRE(pretty(pow(x, integer(2))) == "x**2");
}

// ----- function emission ----------------------------------------------------

TEST_CASE("c_function: scalar polynomial", "[13][c_function]") {
    auto x = symbol("x");
    auto body = pow(x, integer(2)) + integer(3) * x + integer(1);
    auto out = c_function("f", body, {x});
    REQUIRE(out.find("double f(double x)") != std::string::npos);
    REQUIRE(out.find("return") != std::string::npos);
    REQUIRE(out.find("pow(x, 2)") != std::string::npos);
}

TEST_CASE("cxx_function: with std::sin", "[13][cxx_function]") {
    auto x = symbol("x");
    auto body = sin(x);
    auto out = cxx_function("g", body, {x});
    REQUIRE(out.find("double g(double x)") != std::string::npos);
    REQUIRE(out.find("std::sin(x)") != std::string::npos);
}

TEST_CASE("fortran_function: declares args + return type",
          "[13][fortran_function]") {
    auto x = symbol("x");
    auto body = pow(x, integer(2));
    auto out = fortran_function("h", body, {x});
    REQUIRE(out.find("function h(x)") != std::string::npos);
    REQUIRE(out.find("real(8), intent(in) :: x") != std::string::npos);
    REQUIRE(out.find("end function h") != std::string::npos);
}

TEST_CASE("octave_function: y = name(args)",
          "[13][octave_function]") {
    auto x = symbol("x");
    auto body = sin(x);
    auto out = octave_function("k", body, {x});
    REQUIRE(out.find("function y = k(x)") != std::string::npos);
    REQUIRE(out.find("y = sin(x);") != std::string::npos);
    REQUIRE(out.find("end") != std::string::npos);
}

TEST_CASE("c_function: multiple args", "[13][c_function]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto body = x * y + sin(x);
    auto out = c_function("foo", body, {x, y});
    REQUIRE(out.find("double foo(double x, double y)") != std::string::npos);
}
