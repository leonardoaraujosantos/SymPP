// 01 — Calculus: differentiate, integrate, take a limit, expand a series.

#include <sympp/sympp.hpp>
#include <iostream>

int main() {
    using namespace sympp;

    auto x = symbol("x");

    // Derivative: chain + product rule
    auto f = sin(x) * exp(x);
    std::cout << "f(x)        = " << f->str() << "\n";
    std::cout << "f'(x)       = " << diff(f, x)->str() << "\n";

    // Indefinite integral via heurisch / table dispatch
    auto g = integer(2) * x * exp(pow(x, integer(2)));
    std::cout << "∫g dx       = " << integrate(g, x)->str() << "\n";

    // Definite integral via Newton-Leibniz
    auto h = integrate(sin(x), x, S::Zero(), S::Pi());
    std::cout << "∫₀^π sin(x) = " << h->str() << "\n";

    // L'Hôpital limit
    auto lim = limit(sin(x) / x, x, S::Zero());
    std::cout << "lim sin(x)/x at 0 = " << lim->str() << "\n";

    // Taylor series
    auto s = series(exp(x), x, S::Zero(), 5);
    std::cout << "exp(x) ≈ " << s->str() << "\n";
}
