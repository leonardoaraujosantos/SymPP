// 07 — Partial differential equations.
//
// SymPP solves several canonical first- and second-order PDEs in closed form:
//   * first-order linear transport via the method of characteristics,
//   * the heat equation by separation of variables,
//   * the wave equation by d'Alembert's formula.
// Arbitrary functions in a general solution are returned as F(...)/G(...)
// (unspecified `function_symbol`s), exactly as in a hand derivation.

#include <sympp/sympp.hpp>

#include <iostream>

int main() {
    using namespace sympp;

    auto x = symbol("x");
    auto y = symbol("y");
    auto t = symbol("t");

    std::cout << "=== First-order linear (transport) ===\n";

    // Constant-coefficient a·u_x + b·u_y = c, solved by characteristics:
    //   u(x,y) = (c/a)·x + F(b·x − a·y),  F arbitrary.
    // Homogeneous transport 2·u_x + 3·u_y = 0  →  u = F(3x − 2y).
    auto transport =
        pdsolve_first_order_linear(integer(2), integer(3), integer(0), x, y);
    std::cout << "2 u_x + 3 u_y = 0   ->  u = " << transport->str() << "\n";

    // Inhomogeneous u_x + u_y = 1  →  u = x + F(x − y).
    auto transport_c =
        pdsolve_first_order_linear(integer(1), integer(1), integer(1), x, y);
    std::cout << "u_x + u_y = 1       ->  u = " << transport_c->str() << "\n";

    std::cout << "\n=== Heat equation u_t = k u_xx ===\n";

    // Separation of variables gives the product mode
    //   u(x,t) = exp(-k·λ·t)·(A·sin(√λ·x) + B·cos(√λ·x))
    // for a separation constant λ.
    auto k = symbol("k");
    auto lambda = symbol("lambda");
    auto heat = pdsolve_heat(k, lambda, x, t);
    std::cout << "u_t = k u_xx        ->  u = " << heat->str() << "\n";

    std::cout << "\n=== Wave equation u_tt = c^2 u_xx ===\n";

    // d'Alembert's general solution: left- and right-moving arbitrary waves.
    auto c = symbol("c");
    auto wave = pdsolve_wave(c, x, t);
    std::cout << "u_tt = c^2 u_xx     ->  u = " << wave->str() << "\n";

    return 0;
}
