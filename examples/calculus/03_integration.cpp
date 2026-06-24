// Engineering Calculus — 03: Integration.
//
// Accumulation: areas, work, averages, energy. Covers indefinite integrals,
// the Fundamental Theorem (definite integrals), improper integrals to infinity,
// the average value of a signal, and arc length.

#include <sympp/sympp.hpp>

#include <iostream>

int main() {
    using namespace sympp;
    auto x = symbol("x");
    auto oo = S::Infinity();

    std::cout << "=== Indefinite integrals (antiderivatives) ===\n";
    std::cout << "INT x^3 dx       = " << integrate(pow(x, integer(3)), x)->str()
              << "\n";
    std::cout << "INT 1/(1+x^2) dx = "
              << integrate(pow(integer(1) + x * x, integer(-1)), x)->str()
              << "  (arctan)\n";
    std::cout << "INT ln(x) dx     = " << integrate(log(x), x)->str()
              << "  (by parts)\n";
    std::cout << "INT x e^x dx     = " << integrate(x * exp(x), x)->str()
              << "\n";

    std::cout << "\n=== Definite integrals (Fundamental Theorem) ===\n";
    // Area under one hump of a sine.
    std::cout << "INT_0^pi sin(x) dx        = "
              << integrate(sin(x), x, integer(0), S::Pi())->str() << "\n";
    // Work done by a linear spring F = kx from 0 to L (k=1): 1/2 L^2.
    auto L = symbol("L");
    std::cout << "INT_0^L x dx (spring work)= "
              << integrate(x, x, integer(0), L)->str() << "\n";

    std::cout << "\n=== Improper integrals (to infinity) ===\n";
    // Exponential decay: total charge / energy of a transient.
    std::cout << "INT_0^oo e^(-x) dx        = "
              << integrate(exp(mul(integer(-1), x)), x, integer(0), oo)->str()
              << "\n";
    // The Gaussian integral (probability / diffusion).
    std::cout << "INT_0^oo e^(-x^2) dx      = "
              << integrate(exp(mul(integer(-1), pow(x, integer(2)))), x,
                           integer(0), oo)
                     ->str()
              << "  (= sqrt(pi)/2)\n";

    std::cout << "\n=== Average value of a signal ===\n";
    // RMS-style mean of sin^2 over a period: (1/pi) INT_0^pi sin^2 = 1/2.
    auto mean = integrate(pow(sin(x), integer(2)), x, integer(0), S::Pi())
                / S::Pi();
    std::cout << "mean of sin^2 over [0,pi] = " << simplify(mean)->str() << "\n";

    std::cout << "\n=== Arc length ===\n";
    // Length of y = x^2 from 0 to 1:  INT_0^1 sqrt(1 + (y')^2) dx.
    auto yfun = pow(x, integer(2));
    auto integrand = pow(integer(1) + pow(diff(yfun, x), integer(2)),
                         rational(1, 2));
    std::cout << "arclen y=x^2 on [0,1] = "
              << simplify(integrate(integrand, x, integer(0), integer(1)))
                     ->str()
              << "\n";
    return 0;
}
