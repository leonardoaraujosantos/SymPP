// Engineering Calculus — 02: Differentiation.
//
// Rates of change, optimisation and sensitivity analysis. Covers ordinary and
// higher-order derivatives, the product/chain rules, partial derivatives and
// the gradient/Jacobian, plus locating stationary points (optimisation).

#include <sympp/sympp.hpp>
#include <sympp/vector/vector_calculus.hpp>

#include <iostream>

int main() {
    using namespace sympp;
    auto x = symbol("x");
    auto y = symbol("y");

    std::cout << "=== Single-variable derivatives ===\n";
    auto f = mul(exp(x), sin(x));  // product rule
    std::cout << "d/dx [e^x sin x]     = " << diff(f, x)->str() << "\n";
    std::cout << "d/dx [sin(x^2)]      = "
              << diff(sin(pow(x, integer(2))), x)->str() << "  (chain rule)\n";
    // Higher-order: acceleration is the 2nd derivative of position.
    auto pos = pow(x, integer(3)) - integer(2) * x;
    std::cout << "position s(x)        = " << pos->str() << "\n";
    std::cout << "velocity s'(x)       = " << diff(pos, x)->str() << "\n";
    std::cout << "accel.   s''(x)      = " << diff(pos, x, 2)->str() << "\n";

    std::cout << "\n=== Partial derivatives & gradient ===\n";
    auto g = pow(x, integer(2)) * y + sin(x * y);
    std::cout << "g(x,y) = x^2 y + sin(xy)\n";
    std::cout << "  dg/dx = " << diff(g, x)->str() << "\n";
    std::cout << "  dg/dy = " << diff(g, y)->str() << "\n";
    auto grad = sympp::vector::gradient(g, {x, y});
    std::cout << "  grad g = {" << grad[0]->str() << ",  " << grad[1]->str()
              << "}\n";
    // Mixed partials are equal (Clairaut): d2g/dxdy == d2g/dydx.
    std::cout << "  d2g/dxdy = " << diff(diff(g, x), y)->str()
              << "  ==  d2g/dydx = " << diff(diff(g, y), x)->str() << "\n";

    std::cout << "\n=== Optimisation: stationary points ===\n";
    // Minimise the parabola-like h(x) = x^3 - 3x  ->  h'(x) = 3x^2 - 3 = 0.
    auto h = pow(x, integer(3)) - integer(3) * x;
    std::cout << "h(x) = x^3 - 3x,  h'(x) = " << diff(h, x)->str() << "\n";
    auto crit = stationary_points(h, x);
    std::cout << "  stationary points:";
    for (const auto& c : crit) std::cout << " x=" << c->str();
    std::cout << "\n";
    // Second-derivative test at each stationary point.
    auto h2 = diff(h, x, 2);
    for (const auto& c : crit) {
        auto curv = simplify(subs(h2, x, c));
        std::cout << "  h''(" << c->str() << ") = " << curv->str()
                  << (curv->str()[0] == '-' ? "  -> local max\n"
                                            : "  -> local min\n");
    }

    std::cout << "\n=== Linearisation (tangent / small-signal model) ===\n";
    // First-order Taylor about x0 = 0: f(x) ~ f(0) + f'(0) x.
    auto fx = exp(x);
    auto f0 = simplify(subs(fx, x, integer(0)));
    auto fp0 = simplify(subs(diff(fx, x), x, integer(0)));
    std::cout << "e^x ~ " << f0->str() << " + " << fp0->str()
              << "*x   near x=0\n";
    return 0;
}
