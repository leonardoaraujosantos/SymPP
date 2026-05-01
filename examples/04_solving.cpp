// 04 — Equation solving: solveset, nsolve, ODE.

#include <sympp/sympp.hpp>
#include <iostream>

int main() {
    using namespace sympp;

    auto x = symbol("x"), y = symbol("y"), yp = symbol("yp");

    // Univariate polynomial: solve returns roots, solveset wraps as Set
    auto eq = pow(x, integer(2)) - integer(4);
    std::cout << "solveset(x² - 4) = " << solveset(eq, x)->str() << "\n";

    // Trig solveset emits ImageSet for the full periodic family
    std::cout << "solveset(sin(x)) = " << solveset(sin(x), x)->str() << "\n";

    // Numeric Newton in MPFR
    auto r = nsolve(cos(x) - x, x, integer(1), 20);
    std::cout << "cos(x) = x at      " << r->str() << "\n";

    // ODE: y' + y = e^x — linear first-order
    auto ode_eq = yp + y - exp(x);
    auto sol = dsolve_first_order(ode_eq, y, yp, x);
    std::cout << "y'+y=e^x → y = " << sol->str() << "\n";

    // Verify by substitution
    auto residual = checkodesol(ode_eq, sol, y, yp, x);
    std::cout << "residual: " << residual->str() << " (expect 0)\n";
}
