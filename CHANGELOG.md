# Changelog

All notable changes to SymPP. Format: [Keep a Changelog](https://keepachangelog.com/),
versioning per [SemVer](https://semver.org/).

## [0.5.0] — Unreleased

First public-readiness release. 13 of 15 phases shipped, parser landed,
`find_package(SymPP)` works, examples directory and CI in place.

### Added

- **Phase 0 — Foundation**: oracle harness (long-lived Python subprocess),
  CMake build, GMP/MPFR/Catch2 wired.
- **Phase 1 — Core**: `Basic`/`Expr`, full `Number` tower
  (`Integer`/`Rational`/`Float`/`NumberSymbol`), `Symbol`,
  `ImaginaryUnit`, hash-cons cache, `Add`/`Mul`/`Pow` with canonical
  form, operator overloads, `subs`/`xreplace`/`has`/`expand`,
  `evalf` to arbitrary precision via MPFR.
- **Phase 2 — Assumptions**: `AssumptionMask`, `is_real`/`is_positive`/
  `is_integer`/... queries, propagation through Add/Mul/Pow, `refine`.
  (SAT-based reasoning deferred.)
- **Phase 3 — Functions**: 30+ named functions across trig, hyperbolic,
  exponential, complex, integer-rounding, combinatorial, and special
  families.
- **Phase 4 — Polynomials**: univariate `Poly`, GCD via Euclidean
  algorithm, square-free factorization (Yun), rational roots, Cardano
  cubic, Ferrari quartic, `factor` over ℤ via Kronecker, partial
  fractions (`apart`), `cancel`/`together`/`horner`, resultant +
  discriminant, lazy `RootOf` with MPFR bisection evalf.
- **Phase 5 — Simplification**: `simplify` orchestrator, `collect`,
  `cancel`, `powsimp`, `expand_power_*`, `trigsimp`, `radsimp`,
  `sqrtdenest`, `combsimp`, `gammasimp`, `cse`, `nsimplify`.
- **Phase 6 — Calculus**: `diff` (full chain rule), `series`, `limit`
  with L'Hôpital, `Order` (Big-O), `summation`, `pade`,
  `stationary_points`/`minimum`/`maximum`/`is_increasing`/
  `is_decreasing`, `euler_lagrange`, `vertical_asymptotes`,
  `inflection_points`.
- **Phase 7 — Integration**: linearity + table with affine arguments,
  trig reductions, rational integration via `apart`, integration by
  parts, `heurisch` chain-rule-reverse u-substitution,
  `manualintegrate` orchestrator, `vpaintegral` adaptive Simpson in
  MPFR.
- **Phase 8 — Transforms**: Laplace, Fourier, Mellin, Z, sine/cosine —
  all forward + inverse, table-driven.
- **Phase 9 — Linear algebra**: dense `Matrix`, det/inverse/transpose/
  trace/adjugate, conjugate transpose, Frobenius norm,
  `gradient`/`hessian`/`jacobian`/`wronskian`, `rref`/`rank`/
  `nullspace`/`columnspace`/`rowspace`, char poly + eigenvals +
  eigenvects + diagonalize, **LU + QR (Gram-Schmidt) + Cholesky**,
  Hadamard + Kronecker products, Hilbert / Vandermonde / companion /
  rotation constructors, `MatrixSymbol` symbolic expression tree.
- **Phase 10 — Solvers**: `solve`, `solveset`, `linsolve`, `nsolve`
  (MPFR Newton), `solve_univariate_inequality`, `reduce_inequalities`,
  `rsolve`, `nonlinsolve` (resultant + Gröbner). Set algebra:
  `Interval`, `FiniteSet`, `Union`, `Intersection`, `Complement`,
  `Reals`, `Integers`, `ConditionSet`, `ImageSet`. Trig solveset
  emits `ImageSet` for periodic solutions. Linear Diophantine,
  Pythagorean triples, Buchberger Gröbner basis.
- **Phase 11 — ODE/PDE**: `dsolve_first_order` classifier with
  separable / linear / Bernoulli / exact / Riccati / homogeneous /
  Lie autonomous strategies. Higher-order constant-coefficient and
  Cauchy-Euler. Linear systems via eigendecomposition. IVP with
  `apply_ivp`. Verification via `checkodesol`. PDE: first-order
  linear, heat equation, wave equation. DAE Jacobians + structural
  index estimation. Hypergeometric ODE recognition.
- **Phase 12 — Units**: `Dimension` algebra, `Unit`, `Quantity`, full
  SI base + 13 derived units, prefixes (yocto..yotta), CGS + US
  customary, 12 physical constants with exact post-2019-redef
  values, `convert`, affine temperature C/F/K conversions.
- **Phase 13 — Code generation**: `ccode` / `cxxcode` / `fcode` /
  `latex` / `octave_code` / `pretty` printers, `c_function` /
  `cxx_function` / `fortran_function` / `octave_function` emission.
- **Phase 15 — Parser + MATLAB facade**:
  * `parse(string_view)` returning `Expr` — full recursive-descent
    with proper precedence (unary minus weaker than power, **
    right-associative, `^` aliased to `**`), function-call
    dispatch, undefined-function preservation.
  * `sympp::matlab::*` MATLAB-named facade (header-only) —
    variadic `syms`, `sym`, `Int` (renamed from the `int` keyword),
    `diff(f, x, n)` for n-th derivative, `taylor` (= series),
    `vpa` (= evalf), all four transform pairs (`laplace`/
    `ilaplace`/`fourier`/`ifourier`/`ztrans`/`iztrans`), `solve`/
    `dsolve`, the algebra wrappers, and print/codegen wrappers
    (`pretty`, `latex`, `ccode`, `fortran`, `matlabFunction`).

### Build / packaging

- CMake `find_package(SymPP)` integration: install rules,
  `SymPPConfig.cmake.in`, version config. Downstream consumers can
  `find_package(SymPP REQUIRED)` after install.
- `examples/` directory: 5 standalone consumer programs (calculus,
  polynomials, linear algebra, solving, codegen) demonstrating the
  public API and dual-use as in-tree builds or standalone via the
  install.
- GitHub Actions CI: Linux + macOS matrix on gcc + clang. Each
  job builds in Release mode and runs the full test suite (with the
  oracle exercising SymPy 1.14).

### Known limitations / deferred work

See [docs/04-roadmap.md](docs/04-roadmap.md) for the per-phase
deferred-deep list. Highlights:

- Phase 14 (Plotting bridge) dropped from scope.
- Phase 16 (full v1.0 hardening) in progress: Doxygen API docs,
  performance benchmarks, vcpkg/Conan recipes, pre-built binaries
  pending.
- Within shipped phases: full Gruntz limits at infinity, full Risch
  integration, Meijer G method, Lie symmetry classifier, Pantelides
  DAE index reduction, F4/F5 Gröbner, lambdify (LLVM JIT), full
  symbolic SVD, Jordan canonical form, sparse matrices,
  multivariate factorization (Wang), Berlekamp-Zassenhaus,
  hyperexpand. Each is independent and tracked as **Category A**
  in the roadmap.

### Dependencies

- C++20 compiler (gcc 11+, clang 14+, AppleClang 14+)
- CMake ≥ 3.25
- GMP + GMPXX
- MPFR ≥ 4.0
- Python 3.10+ with SymPy 1.14+ (test oracle only — runtime is
  pure C++)

### Testing

```
880 tests / 1724 assertions  all passing
```

Every numeric and structural assertion is cross-checked against
SymPy 1.14 via the oracle harness on `[oracle]`-tagged tests.
