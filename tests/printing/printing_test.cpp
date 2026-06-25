// Phase 13 — printing/codegen tests.

#include <string>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

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

#include "oracle/oracle.hpp"

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

// ----- rust_code ------------------------------------------------------------

TEST_CASE("rust: x^3 → x.powi(3)", "[13][rust]") {
    auto x = symbol("x");
    REQUIRE(rust_code(pow(x, integer(3))) == "x.powi(3)");
}

TEST_CASE("rust: x^y → x.powf(y)", "[13][rust]") {
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(rust_code(pow(x, y)) == "x.powf(y)");
}

TEST_CASE("rust: sqrt(x) → x.sqrt()", "[13][rust]") {
    auto x = symbol("x");
    REQUIRE(rust_code(sqrt(x)) == "x.sqrt()");
}

TEST_CASE("rust: integer literal is f64", "[13][rust]") {
    REQUIRE(rust_code(integer(42)) == "42.0");
}

TEST_CASE("rust: rational is float division", "[13][rust]") {
    REQUIRE(rust_code(rational(1, 2)) == "1.0/2.0");
}

TEST_CASE("rust: sin(x) → x.sin()", "[13][rust]") {
    auto x = symbol("x");
    REQUIRE(rust_code(sin(x)) == "x.sin()");
}

TEST_CASE("rust: exp/log are method calls", "[13][rust]") {
    auto x = symbol("x");
    REQUIRE(rust_code(exp(x)) == "x.exp()");
    REQUIRE(rust_code(log(x)) == "x.ln()");
}

TEST_CASE("rust: pi and e from f64::consts", "[13][rust]") {
    REQUIRE(rust_code(S::Pi()) == "f64::consts::PI");
    REQUIRE(rust_code(S::E()) == "f64::consts::E");
}

TEST_CASE("rust: 1/(1+x^2) uses division", "[13][rust]") {
    auto x = symbol("x");
    auto e = integer(1) / (integer(1) + pow(x, integer(2)));
    REQUIRE(rust_code(e) == "1.0/(x.powi(2) + 1.0)");
}

TEST_CASE("rust: x^3 + 2*x*y + 1", "[13][rust]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = pow(x, integer(3)) + integer(2) * x * y + integer(1);
    REQUIRE(rust_code(e) == "2.0*x*y + x.powi(3) + 1.0");
}

TEST_CASE("rust: exp(-x^2/2)*sin(2*pi*x)", "[13][rust]") {
    auto x = symbol("x");
    auto e = exp(rational(-1, 2) * pow(x, integer(2)))
             * sin(integer(2) * S::Pi() * x);
    REQUIRE(rust_code(e)
            == "(-1.0/2.0*x.powi(2)).exp()*(2.0*x*f64::consts::PI).sin()");
}

// ----- julia_code -----------------------------------------------------------

TEST_CASE("julia: x^3 uses ^", "[13][julia]") {
    auto x = symbol("x");
    REQUIRE(julia_code(pow(x, integer(3))) == "x^3");
}

TEST_CASE("julia: sqrt(x)", "[13][julia]") {
    auto x = symbol("x");
    REQUIRE(julia_code(sqrt(x)) == "sqrt(x)");
}

TEST_CASE("julia: sin(x) is a function call", "[13][julia]") {
    auto x = symbol("x");
    REQUIRE(julia_code(sin(x)) == "sin(x)");
}

TEST_CASE("julia: pi and MathConstants.e", "[13][julia]") {
    REQUIRE(julia_code(S::Pi()) == "pi");
    REQUIRE(julia_code(S::E()) == "MathConstants.e");
}

TEST_CASE("julia: 1/(1+x^2) uses division", "[13][julia]") {
    auto x = symbol("x");
    auto e = integer(1) / (integer(1) + pow(x, integer(2)));
    REQUIRE(julia_code(e) == "1/(x^2 + 1)");
}

TEST_CASE("julia: x^3 + 2*x*y + 1", "[13][julia]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = pow(x, integer(3)) + integer(2) * x * y + integer(1);
    REQUIRE(julia_code(e) == "2*x*y + x^3 + 1");
}

TEST_CASE("julia: exp(-x^2/2)*sin(2*pi*x)", "[13][julia]") {
    auto x = symbol("x");
    auto e = exp(rational(-1, 2) * pow(x, integer(2)))
             * sin(integer(2) * S::Pi() * x);
    REQUIRE(julia_code(e) == "exp(-1/2*x^2)*sin(2*x*pi)");
}

// ----- mathml ---------------------------------------------------------------

TEST_CASE("mathml: symbol", "[13][mathml]") {
    auto x = symbol("x");
    REQUIRE(mathml(x) == "<math><mi>x</mi></math>");
}

TEST_CASE("mathml: integer", "[13][mathml]") {
    REQUIRE(mathml(integer(7)) == "<math><mn>7</mn></math>");
}

TEST_CASE("mathml: rational → mfrac", "[13][mathml]") {
    REQUIRE(mathml(rational(1, 2))
            == "<math><mfrac><mn>1</mn><mn>2</mn></mfrac></math>");
}

TEST_CASE("mathml: x^2 → msup", "[13][mathml]") {
    auto x = symbol("x");
    REQUIRE(mathml(pow(x, integer(2)))
            == "<math><msup><mi>x</mi><mn>2</mn></msup></math>");
}

TEST_CASE("mathml: sqrt(x) → msqrt", "[13][mathml]") {
    auto x = symbol("x");
    REQUIRE(mathml(sqrt(x)) == "<math><msqrt><mi>x</mi></msqrt></math>");
}

TEST_CASE("mathml: sin(x)", "[13][mathml]") {
    auto x = symbol("x");
    REQUIRE(mathml(sin(x))
            == "<math><mrow><mi>sin</mi><mo>&#x2061;</mo>"
               "<mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>");
}

TEST_CASE("mathml: pi", "[13][mathml]") {
    REQUIRE(mathml(S::Pi()) == "<math><mi>&#x3C0;</mi></math>");
}

TEST_CASE("mathml: 1/(1+x^2) → mfrac", "[13][mathml]") {
    auto x = symbol("x");
    auto e = integer(1) / (integer(1) + pow(x, integer(2)));
    REQUIRE(mathml(e)
            == "<math><mfrac><mn>1</mn><mrow><mrow><msup><mi>x</mi>"
               "<mn>2</mn></msup><mo>+</mo><mn>1</mn></mrow></mrow>"
               "</mfrac></math>");
}

TEST_CASE("mathml: x^3 + 2*x*y + 1", "[13][mathml]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = pow(x, integer(3)) + integer(2) * x * y + integer(1);
    REQUIRE(mathml(e)
            == "<math><mrow><mrow><mn>2</mn><mo>&#x2062;</mo><mi>x</mi>"
               "<mo>&#x2062;</mo><mi>y</mi></mrow><mo>+</mo>"
               "<msup><mi>x</mi><mn>3</mn></msup><mo>+</mo><mn>1</mn>"
               "</mrow></math>");
}

// ----- pretty ---------------------------------------------------------------

TEST_CASE("pretty: x + 2", "[13][pretty]") {
    auto x = symbol("x");
    REQUIRE(pretty(x + integer(2)) == "x + 2");
}

TEST_CASE("pretty: x^2 renders as a 2D superscript", "[13][pretty]") {
    auto x = symbol("x");
    REQUIRE(pretty(pow(x, integer(2))) == " 2\nx ");
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

// ----- glsl_code ------------------------------------------------------------

TEST_CASE("glsl: x^3 → pow with float exponent", "[13][glsl]") {
    auto x = symbol("x");
    REQUIRE(glsl_code(pow(x, integer(3))) == "pow(x, 3.0)");
}

TEST_CASE("glsl: sqrt(x)", "[13][glsl]") {
    auto x = symbol("x");
    REQUIRE(glsl_code(sqrt(x)) == "sqrt(x)");
}

TEST_CASE("glsl: ^(-1/2) → inversesqrt", "[13][glsl]") {
    auto x = symbol("x");
    REQUIRE(glsl_code(pow(x, rational(-1, 2))) == "inversesqrt(x)");
}

TEST_CASE("glsl: sin / exp are function calls", "[13][glsl]") {
    auto x = symbol("x");
    REQUIRE(glsl_code(sin(x)) == "sin(x)");
    REQUIRE(glsl_code(exp(x)) == "exp(x)");
}

TEST_CASE("glsl: pi and e use named float constants", "[13][glsl]") {
    REQUIRE(glsl_code(S::Pi()) == "PI");
    REQUIRE(glsl_code(S::E()) == "E");
}

TEST_CASE("glsl: integer and rational literals are floats", "[13][glsl]") {
    REQUIRE(glsl_code(integer(2)) == "2.0");
    REQUIRE(glsl_code(rational(2, 3)) == "2.0/3.0");
}

TEST_CASE("glsl: 1/(1+x^2) uses division", "[13][glsl]") {
    auto x = symbol("x");
    auto e = integer(1) / (integer(1) + pow(x, integer(2)));
    REQUIRE(glsl_code(e) == "1.0/(pow(x, 2.0) + 1.0)");
}

TEST_CASE("glsl: exp(-x^2/2)*sin(2*pi*x)", "[13][glsl]") {
    auto x = symbol("x");
    auto e = exp(rational(-1, 2) * pow(x, integer(2)))
             * sin(integer(2) * S::Pi() * x);
    REQUIRE(glsl_code(e) == "exp(-1.0/2.0*pow(x, 2.0))*sin(2.0*x*PI)");
}

// ----- dot ------------------------------------------------------------------

TEST_CASE("dot: emits a digraph wrapper", "[13][dot]") {
    auto x = symbol("x");
    auto out = dot(x);
    REQUIRE(out.find("digraph {") != std::string::npos);
    REQUIRE(out.find("label=\"Symbol(x)\"") != std::string::npos);
}

TEST_CASE("dot: 2*x node labels and edge count", "[13][dot]") {
    auto x = symbol("x");
    auto out = dot(integer(2) * x);
    // Three nodes: Mul (root), Integer(2), Symbol(x).
    REQUIRE(out.find("label=\"Mul\"") != std::string::npos);
    REQUIRE(out.find("label=\"Integer(2)\"") != std::string::npos);
    REQUIRE(out.find("label=\"Symbol(x)\"") != std::string::npos);
    // Two directed edges (root -> each child).
    std::size_t edges = 0;
    for (std::size_t p = out.find(" -> "); p != std::string::npos;
         p = out.find(" -> ", p + 1)) {
        ++edges;
    }
    REQUIRE(edges == 2);
    // Deterministic: same expression prints identically.
    REQUIRE(dot(integer(2) * x) == out);
}

TEST_CASE("dot: x^2 + 1 tree", "[13][dot]") {
    auto x = symbol("x");
    auto out = dot(pow(x, integer(2)) + integer(1));
    REQUIRE(out.find("label=\"Add\"") != std::string::npos);
    REQUIRE(out.find("label=\"Pow\"") != std::string::npos);
    REQUIRE(out.find("label=\"Symbol(x)\"") != std::string::npos);
    REQUIRE(out.find("label=\"Integer(2)\"") != std::string::npos);
    REQUIRE(out.find("label=\"Integer(1)\"") != std::string::npos);
    // Edges: Add->Pow, Pow->x, Pow->2, Add->1  = 4.
    std::size_t edges = 0;
    for (std::size_t p = out.find(" -> "); p != std::string::npos;
         p = out.find(" -> ", p + 1)) {
        ++edges;
    }
    REQUIRE(edges == 4);
}

// ----- srepr ----------------------------------------------------------------

TEST_CASE("srepr: cross-checked against SymPy", "[13][srepr][oracle]") {
    auto& oracle = sympp::testing::Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");

    // (SymPP expr, SymPy sympify string) — chosen so arg ordering coincides.
    REQUIRE(srepr(x + integer(1)) == oracle.srepr("x+1"));
    REQUIRE(srepr(pow(x, integer(2))) == oracle.srepr("x**2"));
    REQUIRE(srepr(integer(2) * x * y) == oracle.srepr("2*x*y"));
    REQUIRE(srepr(rational(2, 3)) == oracle.srepr("Rational(2,3)"));
    REQUIRE(srepr(sin(x)) == oracle.srepr("sin(x)"));
}

TEST_CASE("srepr: structural forms", "[13][srepr]") {
    auto x = symbol("x");
    REQUIRE(srepr(x + integer(1)) == "Add(Symbol('x'), Integer(1))");
    REQUIRE(srepr(pow(x, integer(2))) == "Pow(Symbol('x'), Integer(2))");
    REQUIRE(srepr(rational(2, 3)) == "Rational(2, 3)");
    REQUIRE(srepr(sin(x)) == "sin(Symbol('x'))");
    REQUIRE(srepr(S::Pi()) == "pi");
    REQUIRE(srepr(S::E()) == "E");
    // 1/(1+x^2) → reciprocal power of a sum.
    REQUIRE(srepr(integer(1) / (integer(1) + pow(x, integer(2))))
            == "Pow(Add(Pow(Symbol('x'), Integer(2)), Integer(1)), Integer(-1))");
    // NOTE: SymPP's canonical Mul ordering places Symbol before the pi
    // NumberSymbol, whereas SymPy's srepr is Mul(Integer(2), pi, Symbol('x')).
    // The set of factors is identical; we assert SymPP's own stable form.
    REQUIRE(srepr(integer(2) * S::Pi() * x)
            == "Mul(Integer(2), Symbol('x'), pi)");
}
