// 05 — Code generation: emit C / C++ / Fortran / LaTeX from an Expr.

#include <sympp/sympp.hpp>
#include <sympp/parsing/parser.hpp>
#include <iostream>

int main() {
    using namespace sympp;
    using namespace sympp::printing;

    // Parse a math expression from a string and emit it in five
    // target languages.
    auto e = parsing::parse("exp(-x^2/2) * sin(2*pi*x)");

    std::cout << "input string:  exp(-x^2/2) * sin(2*pi*x)\n\n";
    std::cout << "ccode:        " << ccode(e) << "\n";
    std::cout << "cxxcode:      " << cxxcode(e) << "\n";
    std::cout << "fcode:        " << fcode(e) << "\n";
    std::cout << "octave_code:  " << octave_code(e) << "\n";
    std::cout << "latex:        " << latex(e) << "\n\n";

    // Wrap as a usable C function
    auto x = symbol("x");
    auto poly = integer(3) * pow(x, integer(2)) + integer(2) * x + integer(1);
    std::cout << "C function:\n" << c_function("f", poly, {x}) << "\n";
}
