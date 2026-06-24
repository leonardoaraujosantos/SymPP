// Engineering Algebra — 02: Equations & systems.
//
// Solving for operating points and intersections. Covers polynomial equations,
// solveset (full solution sets including complex / periodic), linear systems,
// and a small nonlinear system.

#include <sympp/sympp.hpp>

#include <iostream>

int main() {
    using namespace sympp;
    auto x = symbol("x");
    auto y = symbol("y");

    std::cout << "=== Polynomial equations ===\n";
    // Quadratic formula at work.
    std::cout << "solve x^2 = 2:";
    for (const auto& r : solve(pow(x, integer(2)) - integer(2), x))
        std::cout << " " << r->str();
    std::cout << "\n";
    // Cubic (Cardano).
    std::cout << "solve x^3 - 2x - 5 = 0:";
    for (const auto& r :
         solve(add({pow(x, integer(3)), mul(integer(-2), x), integer(-5)}), x))
        std::cout << " " << r->str();
    std::cout << "\n";

    std::cout << "\n=== Solution sets (solveset) ===\n";
    // Complex domain by default: x^2 + 1 = 0 -> {-I, I}.
    std::cout << "solveset(x^2+1)  = "
              << solveset(pow(x, integer(2)) + integer(1), x)->str() << "\n";
    // Periodic trig solution as an ImageSet.
    std::cout << "solveset(sin x)  = " << solveset(sin(x), x)->str() << "\n";
    // Restricted to the reals.
    std::cout << "solveset(x^2+1, Reals) = "
              << solveset(pow(x, integer(2)) + integer(1), x, reals())->str()
              << "\n";

    std::cout << "\n=== Linear system  A x = b ===\n";
    // 2x + y = 5,  x - y = 1   ->  x=2, y=1.
    Matrix A = {{integer(2), integer(1)}, {integer(1), integer(-1)}};
    Matrix b = {{integer(5)}, {integer(1)}};
    Matrix sol = linsolve(A, b);
    std::cout << "2x + y = 5,  x - y = 1  ->  x=" << sol.at(0, 0)->str()
              << ", y=" << sol.at(1, 0)->str() << "\n";

    std::cout << "\n=== Nonlinear system (intersection) ===\n";
    // Circle x^2 + y^2 = 1 meets line y = x.
    auto eqs = nonlinsolve(
        {pow(x, integer(2)) + pow(y, integer(2)) - integer(1), y - x}, {x, y});
    std::cout << "x^2+y^2=1 and y=x ->\n";
    for (const auto& s : eqs) {
        std::cout << "  (x,y) = (" << s[0]->str() << ", " << s[1]->str()
                  << ")\n";
    }
    return 0;
}
