# Changelog

All notable changes to SymPP. Format: [Keep a Changelog](https://keepachangelog.com/),
versioning per [SemVer](https://semver.org/).

## [0.5.0] — Unreleased

First public-readiness release. 14 of 15 phases shipped (Phase 16
partial — `find_package` install + CI matrix done; v1.0 tag pending),
parser landed, MATLAB facade extended, examples directory and CI in
place. 962 tests / 1872 assertions, 307 oracle-validated against
SymPy 1.14, all passing.

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
  * MATLAB facade extension — split into `matlab/parsing.hpp`,
    `matlab/assumptions.hpp`, `matlab/solvers.hpp`,
    `matlab/ode.hpp`, `matlab/pde.hpp` (all re-exported from the
    `matlab/matlab.hpp` umbrella):
    - `str2sym(s)` and `sym(string_view)` — string inputs containing
      operators / parens / whitespace are routed through the parser
      (a bare identifier still becomes a single Symbol).
    - `assume` / `assumeAlso` / `assumptions` / `clearAssumptions`
      / `refresh` — process-wide assumption registry on top of
      SymPP's immutable Symbols. `assume(x, "positive")` registers
      under `x`'s printed name; `x = refresh(x)` returns a fresh
      Symbol carrying the registered mask.
    - `linsolve(A, b)`, `solve({eqs}, {vars})` (routes to
      `nonlinsolve` for 2×2, `nonlinsolve_groebner` for n ≥ 3),
      `nsolve(eq, var, x0, dps=15)`, `vpasolve(eq, var, x0, dps=32)`.
    - `dsolve(eq, y, yp, x)` (first-order classifier),
      `dsolve(eq, y, yp, ypp, x)` (second-order with auto-classify
      to constant-coefficient or Cauchy-Euler via partial diff,
      falls through to an unevaluated `Dsolve(...)` marker for
      non-linear inputs), `dsolve_ivp` companions, `dsolve(A, x)`
      for linear systems.
    - `pdsolve(a, b, c, x, y)`, `pdsolve_variable`, `pdsolve_heat`,
      `pdsolve_wave`.

### SymPy parity — Category A deltas

- **Phase 5 — Fu trig rule extension** (`simplify/simplify.hpp`):
  `trigsimp` now also collapses `cos²(x) − sin²(x) → cos(2x)`,
  `2·sin(x)·cos(x) → sin(2x)`, and the constant-mixed forms
  `1 − 2·sin²(x)` / `2·cos²(x) − 1 → cos(2x)`. Picks Pythagorean
  vs. double-angle by leaf count. New `expand_trig` (angle-addition
  expansion: `sin(a+b)`, `cos(a+b)`, `tan(a+b)`) and `fu` orchestrator
  (competes between identity, `trigsimp`, and `trigsimp ∘ expand_trig`,
  returning the form with the smallest leaf count).
- **Phase 10 — solveset `_invert` chain** (`solvers/solve.cpp`):
  `solveset` now peels recognized function heads from the LHS before
  falling through to the polynomial solver — covering `log`, `exp`,
  `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`, `abs`. Periodic-trig
  cases emit `ImageSet` over ℤ; multi-valued inverses (`acosh`, `abs`)
  emit `FiniteSet` of both branches. Removes the transcendental-equation
  cliff for the common single-head shapes.
- **Phase 11 — Variation of parameters** (`ode/dsolve.hpp`): new
  entry points `dsolve_constant_coeff_nonhomogeneous(coeffs, rhs, x)`
  and `dsolve_cauchy_euler_nonhomogeneous(coeffs, rhs, x)` for
  second-order linear ODEs with non-zero RHS. Uses Wronskian-based
  variation of parameters: y_p = -y₁·∫(y₂·g/(a₂·W)) + y₂·∫(y₁·g/(a₂·W)).
  Wired into `matlab::dsolve(eq, y, yp, ypp, x)` — non-zero RHS now
  routes automatically. Also adds an `exp(a)·exp(b) → exp(a+b)` /
  `exp(a)^n → exp(n·a)` folding helper (internal to `dsolve.cpp`)
  to keep the variation-of-parameters integrand in a closed-form-friendly
  shape.
- **Phase 9 — Jordan canonical form + matrix exponential**
  (`matrices/matrix.hpp`): new `Matrix::jordan_form()` returning
  `(P, J)` such that `A = P·J·P⁻¹`, and `Matrix::exp(t)` returning
  `exp(A·t)` via `P · exp(J·t) · P⁻¹` with closed-form Jordan-block
  exponential. Currently supports chains of length up to 2 (the
  textbook 2×2 defective shape and combinations thereof). Updated
  `dsolve_system(A, x)` to use `A.exp(x)` so defective A is handled
  automatically — the previous eigenvect-only path silently dropped
  `t·exp(λt)` components.
- **Phase 5 — Hypergeometric / Meijer-G infrastructure**
  (`functions/hypergeometric.hpp`, `simplify/hyperexpand.hpp`):
  * `Hyper` and `MeijerG` are now proper `Function` subclasses with
    stable `FunctionId` tags (`Hyper`, `MeijerG`), `rebuild`/`str`/
    `diff_arg` overrides, and Pochhammer-style structural accessors
    (`p`, `q`, `ap`, `bq`, `z`). Replaces the old opaque
    `function_symbol("hyper")` placeholder.
  * Variadic factory `hyper(ap, bq, z)` with auto-evaluation:
    `z == 0 → 1`, parameter cancellation, `₀F₀(z) → exp(z)`,
    `₁F₀(a; ; z) → (1 − z)^(−a)`. The 4-arg overload `hyper(a, b, c, z)`
    routes to the same factory for the ₂F₁ Gauss case.
  * `meijerg(an, ap_rest, bm, bq_rest, z)` factory builds the opaque
    `MeijerG` node. Full Slater-theorem expansion is deferred-deep.
  * `hyperexpand(e)` walks the tree and rewrites recognized closed
    forms: `₁F₁(1; 2; z) → (eᶻ − 1)/z`, `₂F₁(1, 1; 2; z) → −log(1−z)/z`.
    Cancellation-driven identities (e.g. `₂F₁(a, b; b; z) → (1−z)^(−a)`)
    fall out of the factory's auto-eval.
  * `simplify(e)` chain extended to call `hyperexpand` after `sqrtdenest`.

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
