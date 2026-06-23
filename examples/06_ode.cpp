// 06 — Ordinary differential equations.
//
// Builds and solves a spread of ODEs with SymPP's `dsolve` family, then
// verifies each general solution by substituting it back into the equation
// (the residual must be 0) and shows how to pin the free constants with
// initial conditions.
//
// API convention: the unknown function and its derivatives are passed as
// ordinary symbols (y, yp, ypp), except for the unified `dsolve(eq, y(x), x)`
// entry point, which takes an applied FunctionSymbol y(x) and reads the order
// off the highest derivative present.

#include <sympp/sympp.hpp>

#include <iostream>

int main() {
    using namespace sympp;

    auto x = symbol("x");
    auto y = symbol("y");
    auto yp = symbol("yp");  // y'

    std::cout << "=== First-order ODEs ===\n";

    // Separable: y' = x·y  →  y = C·exp(x²/2).
    auto sep = yp - x * y;
    auto sep_sol = dsolve_separable(sep, y, yp, x);
    std::cout << "y' = x*y           ->  y = " << sep_sol->str() << "\n";

    // Linear first-order: y' + y = e^x  →  y = e^x/2 + C·e^(-x).
    auto lin = yp + y - exp(x);
    auto lin_sol = dsolve_first_order(lin, y, yp, x);
    std::cout << "y' + y = e^x       ->  y = " << lin_sol->str() << "\n";
    // Verify: substitute the solution back; residual should be 0.
    std::cout << "   check residual   =  "
              << checkodesol(lin, lin_sol, y, yp, x)->str() << "\n";

    // Bernoulli: y' + y = y²  →  y = 1/(1 + C·e^x).
    auto bern = yp + y - pow(y, integer(2));
    std::cout << "y' + y = y^2       ->  y = "
              << dsolve_bernoulli(bern, y, yp, x)->str() << "\n";

    std::cout << "\n=== Second-order linear, constant coefficients ===\n";

    // Homogeneous y'' − 3y' + 2y = 0. Coefficients are [a0, a1, a2] for
    // a2·y'' + a1·y' + a0·y = 0  →  y = C0·e^(2x) + C1·e^x.
    auto homog = dsolve_constant_coeff({integer(2), integer(-3), integer(1)}, x);
    std::cout << "y'' - 3y' + 2y = 0 ->  y = " << homog->str() << "\n";

    // Harmonic oscillator y'' + y = 0  →  y = C0·cos x + C1·sin x, via the
    // unified entry point with an applied function y(x).
    FunctionSymbol Y("y");
    auto sho = diff(Y(x), x, 2) + Y(x);
    auto sho_sol = dsolve(sho, Y(x), x);
    std::cout << "y'' + y = 0        ->  y = " << sho_sol->str() << "\n";

    // Nonhomogeneous y'' + y = x (variation of parameters)  →  particular x.
    auto forced =
        dsolve_constant_coeff_nonhomogeneous({integer(1), integer(0), integer(1)},
                                             x, x);
    std::cout << "y'' + y = x        ->  y = " << forced->str() << "\n";

    std::cout << "\n=== Cauchy-Euler (equidimensional) ===\n";

    // x²·y'' + x·y' − y = 0  →  y = C0·x + C1/x.
    auto ce = dsolve_cauchy_euler({integer(-1), integer(1), integer(1)}, x);
    std::cout << "x^2 y'' + x y' - y = 0 ->  y = " << ce->str() << "\n";

    std::cout << "\n=== Linear system y' = A*y ===\n";

    // Rotation generator: y1' = y2, y2' = -y1  →  sin/cos columns.
    Matrix A = {{integer(0), integer(1)}, {integer(-1), integer(0)}};
    auto sys = dsolve_system(A, x);
    std::cout << "A = [[0,1],[-1,0]]:\n";
    for (std::size_t i = 0; i < sys.rows(); ++i) {
        std::cout << "   y" << (i + 1) << " = " << sys.at(i, 0)->str() << "\n";
    }

    std::cout << "\n=== Initial value problem ===\n";

    // Pin the free constant of y' + y = e^x with y(0) = 1.
    // General y = e^x/2 + C0·e^(-x); y(0) = 1/2 + C0 = 1 ⇒ C0 = 1/2,
    // so y = (e^x + e^(-x))/2 = cosh x.
    auto ivp = apply_ivp(lin_sol, x, {{integer(0), integer(1)}});  // y(0) = 1
    std::cout << "y' + y = e^x, y(0)=1 ->  y = " << ivp->str() << "\n";

    return 0;
}
