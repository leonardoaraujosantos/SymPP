// Engineering Calculus — 01: Limits & continuity.
//
// Limits underpin every derivative, integral and convergence argument in
// engineering analysis. This example covers ordinary limits, indeterminate
// forms (the L'Hôpital cases), one-sided limits at a pole, and limits at
// infinity (steady-state / asymptotic behaviour).

#include <sympp/sympp.hpp>

#include <iostream>

int main() {
    using namespace sympp;
    auto x = symbol("x");
    auto oo = S::Infinity();

    std::cout << "=== Ordinary & indeterminate limits ===\n";
    // The defining limit of the derivative of sin at 0.
    std::cout << "lim_{x->0} sin(x)/x          = "
              << limit(sin(x) / x, x, integer(0))->str() << "\n";
    // 0/0 indeterminate form (L'Hopital).
    std::cout << "lim_{x->0} (1-cos x)/x^2     = "
              << limit((integer(1) - cos(x)) / pow(x, integer(2)), x, integer(0))
                     ->str()
              << "\n";
    // The number e as a 1^infinity limit.
    std::cout << "lim_{x->oo} (1+1/x)^x        = "
              << limit(pow(integer(1) + integer(1) / x, x), x, oo)->str()
              << "\n";

    std::cout << "\n=== One-sided limits at a pole ===\n";
    // 1/x is +inf from the right, -inf from the left — a genuine discontinuity.
    std::cout << "lim_{x->0+} 1/x              = "
              << limit(integer(1) / x, x, integer(0), +1)->str() << "\n";
    std::cout << "lim_{x->0-} 1/x              = "
              << limit(integer(1) / x, x, integer(0), -1)->str() << "\n";
    // |x|/x: the sign function's jump at 0.
    std::cout << "lim_{x->0+} |x|/x            = "
              << limit(abs(x) / x, x, integer(0), +1)->str() << "\n";
    std::cout << "lim_{x->0-} |x|/x            = "
              << limit(abs(x) / x, x, integer(0), -1)->str() << "\n";

    std::cout << "\n=== Limits at infinity (steady state) ===\n";
    // Rational function: ratio of leading coefficients.
    auto rat = (integer(3) * pow(x, integer(2)) + x)
               / (pow(x, integer(2)) - integer(5));
    std::cout << "lim_{x->oo} (3x^2+x)/(x^2-5) = " << limit(rat, x, oo)->str()
              << "\n";
    // A decaying transient → 0.
    std::cout << "lim_{x->oo} exp(-x)*sin(x)   = "
              << limit(exp(mul(integer(-1), x)) * sin(x), x, oo)->str() << "\n";

    std::cout << "\n=== Continuity check ===\n";
    // f(x) = (x^2 - 1)/(x - 1) has a removable discontinuity at x = 1:
    // the limit exists (=2) even though f(1) is 0/0.
    auto f = (pow(x, integer(2)) - integer(1)) / (x - integer(1));
    std::cout << "f(x) = (x^2-1)/(x-1),  lim_{x->1} f = "
              << limit(f, x, integer(1))->str()
              << "  (removable at x=1)\n";
    return 0;
}
