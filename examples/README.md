# SymPP examples

Standalone consumer programs. Build them from the SymPP tree with
`-DSYMPP_BUILD_EXAMPLES=ON`, or as a standalone project pointed at an installed
SymPP via `find_package(SymPP)`.

## Quick-start programs (repo root of `examples/`)

| File | Topic |
|------|-------|
| `01_calculus.cpp` | diff / integrate / limit / series |
| `02_polynomials.cpp` | factor / roots / apart / Gröbner |
| `03_linear_algebra.cpp` | det / eigenvalues / LU / Jacobian |
| `04_solving.cpp` | solveset / nsolve / first-order ODE |
| `05_codegen.cpp` | C / C++ / Fortran / Octave / LaTeX printers |
| `06_ode.cpp` | ODE families + IVP + Laplace |
| `07_pde.cpp` | transport / heat / wave PDEs |

## Engineering-course series

A subject-organised curriculum suitable for an undergraduate engineering maths
sequence. Each program prints worked results for a coherent topic.

### `calculus/`
1. `01_limits_continuity` — limits, indeterminate forms, one-sided limits, continuity
2. `02_differentiation` — derivatives, partials, gradient, optimisation, linearisation
3. `03_integration` — antiderivatives, definite/improper integrals, average value, arc length
4. `04_series_taylor` — Maclaurin/Taylor series, approximations, closed-form sums
5. `05_vector_calculus` — gradient/divergence/curl/Laplacian + vanishing identities
6. `06_engineering_odes` — RC circuit, mass-spring-damper, forced response, Laplace method

### `algebra/`
1. `01_polynomials` — expand/factor/roots/division/GCD/partial fractions
2. `02_equations_systems` — polynomial & set solving, linear and nonlinear systems
3. `03_simplification_trig` — simplification, trig identities, function-family rewriting
4. `04_complex_numbers` — arithmetic, conjugate/magnitude, Euler's formula, roots of unity

### `linear_algebra/`
1. `01_matrix_operations` — product/transpose/det/trace/inverse/powers
2. `02_linear_systems` — `A x = b`, rank, RREF, null space
3. `03_eigen_vibration` — eigenvalues/vectors, characteristic polynomial, modal analysis
4. `04_decompositions` — LU, QR, SVD (each reconstruction-verified)

Each executable is named `<subject>_<file>` in the build tree, e.g.
`calculus_03_integration`, `linear_algebra_04_decompositions`.
