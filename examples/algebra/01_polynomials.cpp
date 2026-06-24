// Engineering Algebra — 01: Polynomials.
//
// Polynomials model transfer functions, characteristic equations and curve
// fits. Covers expansion, factoring, roots, division, GCD and partial-fraction
// decomposition (the key step before an inverse Laplace transform).

#include <sympp/sympp.hpp>

#include <iostream>

int main() {
    using namespace sympp;
    auto x = symbol("x");

    std::cout << "=== Expand & factor ===\n";
    auto p = pow(x + integer(2), integer(3));
    std::cout << "(x+2)^3 expanded = " << expand(p)->str() << "\n";
    std::cout << "factor(x^4 - 1)  = "
              << factor(pow(x, integer(4)) - integer(1), x)->str() << "\n";
    std::cout << "factor(x^3-6x^2+11x-6) = "
              << factor(add({pow(x, integer(3)),
                             mul(integer(-6), pow(x, integer(2))),
                             mul(integer(11), x), integer(-6)}),
                        x)
                     ->str()
              << "\n";

    std::cout << "\n=== Roots (characteristic equation) ===\n";
    auto charpoly = pow(x, integer(2)) - integer(5) * x + integer(6);
    std::cout << "roots of x^2-5x+6:";
    for (const auto& r : solve(charpoly, x)) std::cout << " " << r->str();
    std::cout << "\n";
    // Complex roots of a stable 2nd-order system.
    std::cout << "roots of x^2+2x+5:";
    for (const auto& r :
         solve(add({pow(x, integer(2)), mul(integer(2), x), integer(5)}), x))
        std::cout << " " << r->str();
    std::cout << "\n";

    std::cout << "\n=== Division & GCD ===\n";
    auto num = pow(x, integer(3)) - integer(1);
    auto den = x - integer(1);
    std::cout << "(x^3-1)/(x-1) = " << quo(num, den)->str() << "  rem "
              << rem(num, den)->str() << "\n";
    std::cout << "gcd(x^2-1, x^2+2x+1) = "
              << gcd(pow(x, integer(2)) - integer(1),
                     add({pow(x, integer(2)), mul(integer(2), x), integer(1)}))
                     ->str()
              << "\n";

    std::cout << "\n=== Partial fractions (pre-inverse-Laplace) ===\n";
    // 1/((x-1)(x+2)) -> A/(x-1) + B/(x+2).
    auto rat = pow(mul(x - integer(1), x + integer(2)), integer(-1));
    std::cout << "1/((x-1)(x+2)) = " << apart(rat, x)->str() << "\n";
    // A repeated pole.
    auto rat2 = pow(x + integer(1), integer(-2))
                * (integer(2) * x + integer(3));
    std::cout << "(2x+3)/(x+1)^2 = " << apart(rat2, x)->str() << "\n";
    return 0;
}
