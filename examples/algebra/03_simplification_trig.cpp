// Engineering Algebra — 03: Simplification, trig & rewriting.
//
// Tidying expressions and converting between equivalent forms — essential when
// matching a derived result to a textbook formula or a transfer-function table.

#include <sympp/sympp.hpp>

#include <iostream>
#include <sympp/core/rewrite.hpp>

int main() {
    using namespace sympp;
    auto x = symbol("x");
    auto n = symbol("n");

    std::cout << "=== Algebraic simplification ===\n";
    std::cout << "sin^2 + cos^2          = "
              << simplify(pow(sin(x), integer(2)) + pow(cos(x), integer(2)))
                     ->str()
              << "\n";
    std::cout << "(x^2-1)/(x-1)          = "
              << simplify((pow(x, integer(2)) - integer(1)) / (x - integer(1)))
                     ->str()
              << "\n";
    std::cout << "e^x * e^y              = "
              << simplify(mul(exp(x), exp(symbol("y"))))->str() << "\n";

    std::cout << "\n=== Trig identities ===\n";
    std::cout << "2 sin(x) cos(x)        = "
              << simplify(mul(integer(2), mul(sin(x), cos(x))))->str()
              << "  (double angle)\n";
    std::cout << "cos^2 - sin^2          = "
              << simplify(pow(cos(x), integer(2)) - pow(sin(x), integer(2)))
                     ->str()
              << "\n";
    std::cout << "tan(x)                 = " << simplify(tan(x))->str()
              << "  -> sin/cos: "
              << rewrite(tan(x), RewriteTarget::SinCos)->str() << "\n";

    std::cout << "\n=== Rewriting between function families ===\n";
    // Euler: trig in terms of complex exponentials (phasor analysis).
    std::cout << "sin(x) -> exp form: "
              << rewrite(sin(x), RewriteTarget::Exp)->str() << "\n";
    // Combinatorial: factorial <-> Gamma (special-function tables).
    std::cout << "n!     -> Gamma:    "
              << rewrite(factorial(n), RewriteTarget::Gamma)->str() << "\n";
    std::cout << "C(n,2) -> factorial:"
              << rewrite(binomial(n, integer(2)), RewriteTarget::Factorial)
                     ->str()
              << "\n";

    std::cout << "\n=== Expansion ===\n";
    std::cout << "expand (x+1)^2 (x-1) = "
              << expand(mul(pow(x + integer(1), integer(2)), x - integer(1)))
                     ->str()
              << "\n";
    return 0;
}
