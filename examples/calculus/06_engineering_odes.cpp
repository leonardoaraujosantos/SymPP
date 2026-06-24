// Engineering Calculus — 06: Differential equations in engineering.
//
// The two workhorse models of linear systems:
//   * first-order RC circuit (exponential charging),
//   * second-order mass-spring-damper (oscillation / decay),
// solved both directly (dsolve) and via the Laplace transform method used in
// controls and circuit analysis.

#include <sympp/sympp.hpp>

#include <iostream>

int main() {
    using namespace sympp;
    auto t = symbol("t");
    auto s = symbol("s");
    auto y = symbol("y");
    auto yp = symbol("yp");  // y'

    std::cout << "=== First-order: RC circuit charging ===\n";
    // RC v' + v = Vs  (here RC=1, Vs=1):  v' + v = 1.
    auto rc = yp + y - integer(1);
    auto vsol = dsolve_first_order(rc, y, yp, t);
    std::cout << "v' + v = 1  ->  v(t) = " << vsol->str() << "\n";
    // Pin the initial condition v(0) = 0 (capacitor initially uncharged).
    auto vivp = apply_ivp(vsol, t, {{integer(0), integer(0)}});
    std::cout << "with v(0)=0 ->  v(t) = " << vivp->str()
              << "   (-> 1 as t->oo)\n";

    std::cout << "\n=== Second-order: mass-spring-damper ===\n";
    // m y'' + c y' + k y = 0. Coefficients [k, c, m] for the solver.
    // Underdamped example: y'' + 2y' + 5y = 0  (roots -1 +/- 2i).
    auto under = dsolve_constant_coeff({integer(5), integer(2), integer(1)}, t);
    std::cout << "y'' + 2y' + 5y = 0 (underdamped) ->\n  y(t) = " << under->str()
              << "\n";
    // Critically damped: y'' + 2y' + y = 0  (double root -1).
    auto crit = dsolve_constant_coeff({integer(1), integer(2), integer(1)}, t);
    std::cout << "y'' + 2y' + y = 0  (critically damped) ->\n  y(t) = "
              << crit->str() << "\n";

    std::cout << "\n=== Forced response (variation of parameters) ===\n";
    // Undamped oscillator driven by a ramp:  y'' + y = t.
    auto forced =
        dsolve_constant_coeff_nonhomogeneous({integer(1), integer(0), integer(1)},
                                             t, t);
    std::cout << "y'' + y = t  ->  y(t) = " << forced->str() << "\n";

    std::cout << "\n=== Laplace transform method (controls) ===\n";
    // Transfer-function building blocks.
    std::cout << "L{1}        = "
              << laplace_transform(integer(1), t, s)->str() << "\n";
    std::cout << "L{e^(-t)}   = "
              << laplace_transform(exp(mul(integer(-1), t)), t, s)->str()
              << "\n";
    std::cout << "L{sin(t)}   = " << laplace_transform(sin(t), t, s)->str()
              << "\n";
    // Inverse transform of a damped-oscillator transfer function.
    auto F = pow(pow(s + integer(1), integer(2)) + integer(1), integer(-1));
    std::cout << "L^-1{1/((s+1)^2+1)} = "
              << inverse_laplace_transform(F, s, t)->str()
              << "   (damped sine)\n";
    return 0;
}
