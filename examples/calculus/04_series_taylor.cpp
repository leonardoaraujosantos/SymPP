// Engineering Calculus — 04: Series & Taylor approximation.
//
// Polynomial approximation of nonlinear behaviour (small-signal models),
// truncation error, and closed-form sums used in signals and DSP.

#include <sympp/sympp.hpp>

#include <iostream>

int main() {
    using namespace sympp;
    auto x = symbol("x");
    auto n = symbol("n");
    auto oo = S::Infinity();

    std::cout << "=== Maclaurin series (expansion about 0) ===\n";
    std::cout << "e^x      = " << series(exp(x), x, integer(0), 6)->str()
              << "\n";
    std::cout << "sin(x)   = " << series(sin(x), x, integer(0), 8)->str()
              << "\n";
    std::cout << "cos(x)   = " << series(cos(x), x, integer(0), 7)->str()
              << "\n";
    std::cout << "1/(1-x)  = "
              << series(pow(integer(1) - x, integer(-1)), x, integer(0), 5)
                     ->str()
              << "  (geometric)\n";
    std::cout << "ln(1+x)  = "
              << series(log(integer(1) + x), x, integer(0), 5)->str() << "\n";

    std::cout << "\n=== Small-angle & engineering approximations ===\n";
    // sqrt(1+x) ~ 1 + x/2  (binomial / first-order).
    std::cout << "sqrt(1+x) = "
              << series(pow(integer(1) + x, rational(1, 2)), x, integer(0), 3)
                     ->str()
              << "\n";
    // tan(x) ~ x + x^3/3 (used in control / optics).
    std::cout << "tan(x)    = " << series(tan(x), x, integer(0), 6)->str()
              << "\n";

    std::cout << "\n=== Closed-form sums (DSP / convergence) ===\n";
    // Geometric series sum (settling of a feedback loop).
    std::cout << "SUM_{k=0}^oo (1/2)^k     = "
              << summation(pow(rational(1, 2), n), n, integer(0), oo)->str()
              << "\n";
    // Basel problem: SUM 1/n^2 = pi^2/6.
    std::cout << "SUM_{n=1}^oo 1/n^2       = "
              << summation(pow(n, integer(-2)), n, integer(1), oo)->str()
              << "\n";
    // Finite power sum (Faulhaber).
    std::cout << "SUM_{k=1}^n k^2          = "
              << summation(pow(symbol("k"), integer(2)), symbol("k"), integer(1),
                           n)
                     ->str()
              << "\n";
    // Exponential generating series.
    std::cout << "SUM_{k=0}^oo x^k/k!      = "
              << summation(mul(pow(x, symbol("k")),
                               pow(factorial(symbol("k")), integer(-1))),
                           symbol("k"), integer(0), oo)
                     ->str()
              << "  (= e^x)\n";
    return 0;
}
