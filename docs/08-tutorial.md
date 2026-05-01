# 08 — Tutorial

Hands-on walkthrough of SymPP's shipped features. Every snippet below
is a complete, compilable program — drop into a `.cpp` file, link
against `sympp` + `gmpxx` + `gmp` + `mpfr`, and run.

```cmake
# Minimal CMakeLists.txt for a SymPP consumer:
cmake_minimum_required(VERSION 3.25)
project(my_app CXX)
add_subdirectory(/path/to/SymPP build_sympp EXCLUDE_FROM_ALL)
add_executable(my_app my_app.cpp)
target_link_libraries(my_app PRIVATE SymPP::sympp)
```

Or compile directly (assuming SymPP built at `build/`):

```bash
c++ -std=c++20 -I path/to/SymPP/include -I /opt/homebrew/include \
    my_app.cpp \
    -L path/to/SymPP/build/src -L /opt/homebrew/lib \
    -lsympp -lgmpxx -lgmp -lmpfr -o my_app
```

The single header `<sympp/sympp.hpp>` pulls in the full public API.

---

## 1. Symbols and basic algebra

```cpp
#include <sympp/sympp.hpp>
#include <iostream>

int main() {
    using namespace sympp;

    auto x = symbol("x");
    auto y = symbol("y");

    auto e = (x + y) * (x - y);
    std::cout << e->str() << "\n";          // (x - y)*(x + y)
    std::cout << expand(e)->str() << "\n";  // -y**2 + x**2

    // Arbitrary-precision integers and rationals
    auto big = pow(integer(2), integer(64));
    std::cout << big->str() << "\n";        // 18446744073709551616
    std::cout << (rational(1, 3) + rational(1, 6))->str() << "\n";  // 1/2

    // Symbolic constants
    std::cout << evalf(S::Pi(), 30)->str() << "\n";
    // 3.14159265358979323846264338328
}
```

## 2. Substitution and assumptions

```cpp
auto x = symbol("x", AssumptionMask{}.set_positive(true));
std::cout << is_positive(pow(x, integer(2)) + integer(1)).value() << "\n";  // 1

auto e = sin(x) + cos(x);
auto e_at_pi = subs(e, x, S::Pi());
std::cout << simplify(e_at_pi)->str() << "\n";   // -1
```

## 3. Calculus — diff, integrate, series, limit

```cpp
auto x = symbol("x");

// Derivatives — full chain/product/power rules
auto f = sin(pow(x, integer(2)));
std::cout << diff(f, x)->str() << "\n";
// 2*x*cos(x**2)

// Antiderivatives — table + linear u-sub + trig + parts + heurisch
std::cout << integrate(integer(2)*x*exp(pow(x, integer(2))), x)->str() << "\n";
// exp(x**2)

// Definite integrals via Newton-Leibniz
auto I = integrate(sin(x), x, S::Zero(), S::Pi());
std::cout << I->str() << "\n";                   // 2

// Numeric quadrature when symbolic fails
auto Inum = vpaintegral(exp(pow(x, integer(2))), x,
                          S::Zero(), S::One(), 15);
std::cout << Inum->str() << "\n";                 // 1.46265174590718

// Taylor series
auto s = series(exp(x), x, S::Zero(), 5);
std::cout << s->str() << "\n";
// 1 + x + x**2/2 + x**3/6 + x**4/24

// Limits — L'Hôpital handles 0/0
std::cout << limit(sin(x)/x, x, S::Zero())->str() << "\n";  // 1
std::cout << limit((integer(1) - cos(x))/pow(x, integer(2)),
                   x, S::Zero())->str() << "\n";            // 1/2

// Summation
auto k = symbol("k"), n = symbol("n");
std::cout << summation(pow(k, integer(2)), k,
                       integer(1), n)->str() << "\n";
// n*(n+1)*(2*n+1)/6

// Padé approximants
auto p = pade(exp(x), x, 1, 1);
std::cout << p->str() << "\n";                    // (1 + x/2)/(1 - x/2)

// Optimization helpers
auto crits = stationary_points(pow(x, integer(3)) - integer(3) * x, x);
// crits = {1, -1}
```

## 4. Polynomials

```cpp
auto x = symbol("x");

// Factor over ℤ
std::cout << factor(pow(x, integer(4)) - integer(1), x)->str() << "\n";
// (x - 1)*(x + 1)*(x**2 + 1)

// Roots — closed-form for deg ≤ 4
auto roots4 = solve(
    pow(x, integer(4)) - integer(5)*pow(x, integer(2)) + integer(4), x);
for (auto& r : roots4) std::cout << r->str() << "\n";  // ±1, ±2

// Lazy form for deg ≥ 5 with no rational roots
auto quintic = pow(x, integer(5)) - x - integer(1);
auto laz = root_of(quintic, x, 0);
std::cout << evalf(laz, 20)->str() << "\n";
// 1.16730397826141868426

// Partial fractions
auto rat = (integer(3) * x + integer(5))
           / ((x - integer(1)) * (x + integer(2)));
std::cout << apart(rat, x)->str() << "\n";
// 8/(3*(x - 1)) + 1/(3*(x + 2))

// Resultant / discriminant
auto res = resultant(pow(x, integer(2)) - integer(2),
                     x - integer(1), x);
std::cout << res->str() << "\n";   // -1   (no shared root)

// Multivariate Gröbner basis
auto y = symbol("y");
auto basis = groebner({x*y - integer(1), x + y - integer(3)}, {x, y});
for (auto& g : basis) std::cout << g->str() << "\n";
// x + y - 3
// y**2 - 3*y + 1
```

## 5. Simplification

```cpp
auto x = symbol("x"), y = symbol("y");

// Pythagorean identity
std::cout << simplify(pow(sin(x), integer(2))
                       + pow(cos(x), integer(2)))->str() << "\n";
// 1

// Denesting nested radicals
auto inner = integer(3) + integer(2) * sqrt(integer(2));
std::cout << sqrtdenest(sqrt(inner))->str() << "\n";
// 1 + sqrt(2)

// Rationalize denominators
auto e = pow(integer(1) + sqrt(integer(2)), integer(-1));
std::cout << radsimp(e)->str() << "\n";    // sqrt(2) - 1

// Common subexpression elimination
auto big = pow(x + y, integer(2)) + pow(x + y, integer(3));
auto cse_result = cse(big);
for (auto& [tmp, def] : cse_result.substitutions) {
    std::cout << tmp->str() << " = " << def->str() << "\n";
}
std::cout << "= " << cse_result.expr->str() << "\n";
// _cse_0 = x + y
// = _cse_0**2 + _cse_0**3

// Numeric → exact recovery
std::cout << nsimplify(make<Float>("0.78539816339745", 15))->str() << "\n";
// pi/4
```

## 6. Linear algebra

```cpp
Matrix A = {{integer(2), integer(1)}, {integer(1), integer(3)}};

std::cout << "det = "  << A.det()->str()  << "\n";  // 5
std::cout << "trace = " << A.trace()->str() << "\n"; // 5

// Eigendecomposition + diagonalization
auto [P, D] = A.diagonalize();
// P · D · P^(-1) reconstructs A

// Decompositions
auto [L, U] = A.lu();        // L·U = A
auto [Q, R] = A.qr();        // Q·R = A, Q orthonormal columns
auto chol = A.cholesky();    // chol·cholᵀ = A (A symmetric SPD)

// Subspaces
auto m = Matrix({{integer(1), integer(2), integer(3)},
                 {integer(2), integer(4), integer(6)},
                 {integer(0), integer(0), integer(1)}});
std::cout << "rank = " << m.rank() << "\n";  // 2

// Differential matrices
auto x = symbol("x"), y = symbol("y");
auto J = jacobian({pow(x, integer(2)) + y, integer(2)*x*y}, {x, y});
// J = [[2x, 1], [2y, 2x]]

// Symbolic matrix expression tree (no entry-by-entry materialization)
auto Ms = matrix_symbol("M", 3, 4);
auto Ns = matrix_symbol("N", 4, 2);
auto product = matmul(Ms, Ns);   // shape 3×2, no actual computation
std::cout << product->str() << "\n";   // M*N
```

## 7. Equation solvers

```cpp
auto x = symbol("x"), y = symbol("y");

// Univariate
auto roots = solve(pow(x, integer(2)) - integer(5)*x + integer(6), x);
// roots = {2, 3}

// solveset returns Set objects
auto S = solveset(pow(x, integer(2)) - integer(4), x,
                   interval(integer(0), integer(10)));
// FiniteSet({2})  — domain filter rejects -2

// Trig solveset → ImageSet over Integers (full periodic family)
auto sols = solveset(sin(x), x);
std::cout << sols->str() << "\n";
// ImageSet(Lambda(__n, __n*pi), Integers)

// Numeric Newton in MPFR
auto cosx_eq_x = nsolve(cos(x) - x, x, integer(1), 30);
std::cout << cosx_eq_x->str() << "\n";
// 0.739085133215160641655312087674

// Inequalities
auto ineq = reduce_inequalities({gt(x, integer(0)), lt(x, integer(5))},
                                  x, true);  // AND
std::cout << ineq->str() << "\n";   // intersection of (0,∞) and (-∞,5)

// Linear systems via matrix inverse
Matrix Aleq = {{integer(1), integer(2)}, {integer(3), integer(4)}};
Matrix beq = {{integer(5)}, {integer(11)}};
auto x_sol = linsolve(Aleq, beq);   // x = 1, y = 2

// Recurrences — char poly + general solution
auto n = symbol("n");
auto y_n = rsolve({integer(-2), integer(1)}, n);    // y(n+1) = 2*y(n)
std::cout << y_n->str() << "\n";   // __C0*2**n

// Multivariate nonlinear via Gröbner
auto sols_nl = nonlinsolve_groebner(
    {x*y - integer(1), x + y - integer(3)}, {x, y});
// 2 solutions: ((3 ± √5)/2, (3 ∓ √5)/2)

// Diophantine
auto dio = linear_diophantine(integer(3), integer(5), integer(1));
// → (x0, y0) with 3*x0 + 5*y0 = 1
```

## 8. ODEs and PDEs

```cpp
auto x = symbol("x"), y = symbol("y"), yp = symbol("yp"), ypp = symbol("ypp");

// First-order — classifier dispatches automatically
auto sol1 = dsolve_first_order(yp + y - exp(x), y, yp, x);
// y' + y = e^x  →  e^x/2 + C·e^(-x)

// Constant-coefficient (any order)
//   y'' - 3y' + 2y = 0  with coefficients [a₀, a₁, a₂] = [2, -3, 1]
auto sol2 = dsolve_constant_coeff({integer(2), integer(-3), integer(1)}, x);
// C₀·e^x + C₁·e^(2x)

// Cauchy-Euler  x²y'' - xy' + y = 0
auto sol3 = dsolve_cauchy_euler({integer(1), integer(-1), integer(1)}, x);

// Riccati  y' = 1 + y²
auto sol4 = dsolve_riccati(integer(1), integer(0), integer(1), x);
// → tan(x + C)

// Verify any solution by substituting back
auto residual = checkodesol(yp + y - exp(x), sol1, y, yp, x);
// residual should be 0

// IVP — apply initial condition
auto sol_ivp = apply_ivp(sol1, x, {{S::Zero(), integer(1)}});
// y(0) = 1 fixes the integration constant

// Linear systems  y' = A·y
Matrix A = {{integer(2), integer(0)}, {integer(0), integer(3)}};
auto sys_sol = dsolve_system(A, x);
// 2-vector of independent exponentials

// PDE — heat equation u_t = k·u_xx via separation of variables
auto t = symbol("t"), k = symbol("k"), lambda = symbol("lambda");
auto u = pdsolve_heat(k, lambda, x, t);
// exp(-k·λ·t) · (A·sin(√λ·x) + B·cos(√λ·x))

// Wave equation u_tt = c²·u_xx via d'Alembert
auto c = symbol("c");
auto wave = pdsolve_wave(c, x, t);
// F(x - c·t) + G(x + c·t)
```

## 9. Transforms

```cpp
auto t = symbol("t"), s = symbol("s"), w = symbol("w"), n = symbol("n"), z = symbol("z");
auto a = symbol("a");

// Laplace pairs
std::cout << laplace_transform(exp(a*t), t, s)->str() << "\n";
// (s - a)**(-1)

// Inverse via apart + table lookup
auto F = pow(s - integer(2), integer(-3));
std::cout << inverse_laplace_transform(F, s, t)->str() << "\n";
// t**2 * exp(2*t) / 2

// Fourier (engineering convention F(ω) = ∫ f(t) e^(-iωt) dt)
std::cout << fourier_transform(dirac_delta(t), t, w)->str() << "\n";  // 1
std::cout << fourier_transform(integer(5), t, w)->str() << "\n";
// 10*pi*DiracDelta(w)

// Mellin
std::cout << mellin_transform(exp(mul(S::NegativeOne(), t)),
                               t, s)->str() << "\n";   // gamma(s)

// Z-transform — discrete-time
std::cout << z_transform(pow(a, n), n, z)->str() << "\n";
// z/(z - a)

// Half-line sine / cosine
std::cout << sine_transform(exp(integer(-2)*t), t, w)->str() << "\n";
// w/(4 + w**2)
```

## 10. Units and physical constants

```cpp
using namespace sympp::physics;

// Quantity arithmetic — dimensions auto-checked
Quantity m(integer(10), kilogram());
Quantity acc(rational(98, 10), meter() / second().pow(2));
auto F = m * acc;
std::cout << F.value()->str() << " in " << F.unit().str() << "\n";
// 98 in kg*m/s**2

// Verify dimensional consistency
F.unit().dimension() == newton().dimension();   // true

// Conversions
auto km = prefix_unit(kilo(), "k", meter());
auto a_mile = Quantity(integer(1), mile());
auto in_km = convert(a_mile, km);
std::cout << in_km.value()->str() << "\n";   // 1.609344

// Affine temperature (NOT multiplicative)
std::cout << celsius_to_kelvin(integer(100))->str() << "\n";   // 7463/20
std::cout << fahrenheit_to_celsius(integer(212))->str() << "\n"; // 100

// Physical constants — exact post-2019-redef values
auto c = speed_of_light();              // 299792458 m/s exact
auto h = planck();                      // 6.62607015e-34 J·s exact
auto Na = avogadro();                   // 6.02214076e+23 1/mol exact

// Mixed-unit add — auto-converts
Quantity dist(integer(5), meter());
Quantity extra(integer(50), centimeter());
auto total = dist + extra;
std::cout << total.value()->str() << " " << total.unit().str() << "\n";
// 11/2 m
```

## 11. Code generation

```cpp
using namespace sympp::printing;

auto x = symbol("x"), y = symbol("y");
auto e = pow(x, integer(2)) + sin(x) * y;

std::cout << ccode(e)        << "\n";  // pow(x, 2) + sin(x)*y
std::cout << cxxcode(e)      << "\n";  // std::pow(x, 2) + std::sin(x)*y
std::cout << fcode(e)        << "\n";  // x**2.0d0 + sin(x)*y
std::cout << latex(e)        << "\n";  // x^{2} + \sin(x) \cdot y
std::cout << octave_code(e)  << "\n";  // x.^2 + sin(x).*y
std::cout << pretty(e)       << "\n";  // x**2 + sin(x)*y

// Greek symbol auto-mapping in LaTeX
std::cout << latex(symbol("alpha")) << "\n";    // \alpha
std::cout << latex(rational(1, 2)) << "\n";     // \frac{1}{2}
std::cout << latex(sqrt(x)) << "\n";            // \sqrt{x}

// Wrap as a function
auto poly = pow(x, integer(2)) + integer(3)*x + integer(1);
std::cout << c_function("f", poly, {x}) << "\n";
// double f(double x) {
//     return pow(x, 2) + 3*x + 1;
// }

std::cout << fortran_function("h", pow(x, integer(2)), {x}) << "\n";
// function h(x)
//     implicit none
//     real(8) :: h
//     real(8), intent(in) :: x
//     h = x**2.0d0
// end function h
```

## 12. Cross-validation against SymPy

The `[oracle]`-tagged tests run every numeric and structural assertion
through a Python subprocess running SymPy 1.14. From a consumer's
point of view this is invisible — the SymPP API never depends on
Python at runtime. The oracle runs only in the test harness.

If you want the same level of paranoia in your own code:

```cpp
#include <sympp/sympp.hpp>
// (oracle headers live in tests/oracle/, not the public API — copy
//  oracle.hpp + oracle.cpp + oracle.py into your test setup if you
//  want SymPy cross-checks in your own integration tests).
```

---

## Where to go next

- [04-roadmap.md](04-roadmap.md) — what's deferred, why, and the path
  to full SymPy parity.
- [02-architecture.md](02-architecture.md) — internals: hash-cons,
  visitor patterns, dispatch.
- [05-validation-strategy.md](05-validation-strategy.md) — the
  oracle harness and how every algorithm is cross-checked.
- [03-feature-mapping.md](03-feature-mapping.md) — MATLAB function
  → SymPy reference table for porting MATLAB scripts.
