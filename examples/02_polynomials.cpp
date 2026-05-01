// 02 — Polynomials: factor, find roots, partial fractions, Gröbner.

#include <sympp/sympp.hpp>
#include <iostream>

int main() {
    using namespace sympp;

    auto x = symbol("x");
    auto y = symbol("y");

    // Factor over ℤ
    auto p = pow(x, integer(4)) - integer(1);
    std::cout << "factor(x⁴ - 1)  = " << factor(p, x)->str() << "\n";

    // Roots — Cardano + Ferrari for deg ≤ 4
    auto rs = solve(pow(x, integer(2)) - integer(5)*x + integer(6), x);
    std::cout << "roots of x² - 5x + 6: ";
    for (auto& r : rs) std::cout << r->str() << " ";
    std::cout << "\n";

    // Partial fractions
    auto rat = (integer(3) * x + integer(5))
               / ((x - integer(1)) * (x + integer(2)));
    std::cout << "apart           = " << apart(rat, x)->str() << "\n";

    // Multivariate Gröbner basis (Buchberger, lex order)
    auto basis = groebner({x*y - integer(1), x + y - integer(3)}, {x, y});
    std::cout << "groebner basis:\n";
    for (auto& g : basis) std::cout << "  " << g->str() << "\n";
}
