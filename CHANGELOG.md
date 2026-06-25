# Changelog

All notable changes to SymPP. Format: [Keep a Changelog](https://keepachangelog.com/),
versioning per [SemVer](https://semver.org/).

## [1.1.0] — 2026-06-25

Parity push toward SymPy. 1817 tests / 7702 assertions, 764 oracle-validated,
all passing. A 100-case Tier 1 + Tier 2 acceptance suite (`just parity`) is at
SymPy parity. These are increments toward SymPy parity, not full module
completion.

### Added

- **Code generation (Phase 13 complete)**: `rust_code`, `julia_code`, `mathml`,
  **`glsl_code`** (GLSL shader syntax), **`dot`** (Graphviz expression-tree
  digraph), and **`srepr`** (SymPy constructor form, oracle-validated) printers;
  a structured **codegen AST** (`Assignment`/`CodeBlock`/`FunctionDefinition`
  with CSE-based `emit_c`/`emit_cxx`); and **`autowrap`** — compile generated C
  through the system toolchain into a `dlopen`'d native callable.
- **Tooling**: a `justfile` task runner (build, unit vs. integration test split
  via the `[oracle]` tag, examples, parity gate, bench, docs, ci); an
  engineering-course example series under `examples/{calculus,algebra,linear_algebra}`.
- **Sets / solvers**: `Naturals`/`Naturals0` domain (ℕ ⊂ ℤ ⊂ ℝ ⊂ ℂ); `solveset`
  defaults to the complex domain (SymPy parity); explicit `reals()` still drops
  complex roots.
- **Integration**: Beta-integral binomial sums `Σ C(n,k)/(k+m)`; definite
  erf/erfc integrals over [0,∞); Risch log-integrals with coefficiented power
  denominators; a broad rational/log/special-function integral table.
- **Polynomials**: multivariate Wang factorization up to 4 variables (content
  extraction, linear-form roots, composite difference-of-squares); a cyclotomic
  `factor(x^N−1) = ∏_{d|N} Φ_d(x)` fast path (x⁶⁰−1: >30 s → 6 ms).
- **Summation**: Gosper's algorithm and Zeilberger creative telescoping
  (`Σ k²C(n,k)`, `Σ C(2n,2k)`); symbolic `Product`.
- **Combinatorics / group theory**: Pólya cycle index + `necklaces`; Sylow
  p-subgroups; conjugacy classes, center, derived series, `is_solvable`/
  `is_nilpotent`/`is_simple`/`abelian_invariants`.
- **Physics**: single + multi-mode bosonic/fermionic Fock operators and
  **normal ordering** of operator words (second quantization).
- **hyperexpand**: squared-argument elementary forms and polylog `₃F₂`.
- **Limits**: gamma-ratio asymptotics, competing infinity-power products, and
  exp-of-exp tower values.
- **rewrite**: `RewriteTarget::Gamma`/`Factorial` (factorial ⇄ Γ, binomial,
  rising/falling factorial).

### Fixed

- Install the full public header set (18 module headers were missing from the
  CMake `FILE_SET`), so the `find_package` consumer build no longer fails on
  `sympp/core/rewrite.hpp` — fixes CI's consumer-build step.
- A SIGFPE in limits (GMP exact-power overflow), a `solveset` non-expand bug,
  and a Gosper hang on parametric summands.

## [1.0.0] — 2026-06-23

First tagged release. 1705 tests / 7032 assertions, 740 oracle-validated against
SymPy, all passing; 15 of 15 phases shipped. Packaged on three channels (CMake
`find_package`, vcpkg, Conan) with a benchmark harness, Doxygen target and CI.

### Added — release engineering & late features

- **Packaging**: Conan 2 recipe (`conanfile.py`) alongside the existing
  `vcpkg.json`; optional Doxygen API-doc target (`-DSYMPP_BUILD_DOCS=ON` →
  `sympp_docs`); dependency-free benchmark harness (`sympp_bench`).
- **Meijer-G engine**: generic Slater reduction, Mellin-transform `∫₀^∞`,
  function recognition with the η·xᶜ substitution (Gaussian/Dirichlet/Fresnel),
  confluent → modified Bessel K, and the two-G `∫G₁·G₂` Mellin–Parseval
  convolution.
- **Special functions**: Bessel `besselj/bessely/besseli/besselk`; Hankel
  transform; Bell and tribonacci numbers.
- **Number theory**: CRT, discrete log, linear Diophantine, Pell, sums of
  two/three/four squares, Legendre symbol, quadratic residues, Carmichael λ,
  `nextprime`/`prevprime`/`primorial`/`multiplicity`.
- **Algebra/calculus**: Berlekamp–Zassenhaus with Hensel lifting, general
  Jordan form, full SVD, DomainMatrix, non-commutative `Mul`, exponential Risch
  differential equation, Gruntz termination guarantee.
- **Physics**: arbitrary-spin operators, Wigner 3-j/6-j/9-j, Racah, Gaunt,
  Dirac matrices, hydrogen/QHO states, qubit gates, Jordan–Wigner fermions,
  multi-coordinate Lagrangian EOM.
- **Other modules**: ECC (ECDH/ECDSA), combinatorics (permutation groups,
  orbits, Pólya/Burnside), LaTeX parser, discrete transforms.

## [0.5.0] — development baseline (folded into 1.0.0)

Public-readiness milestone. 14 of 15 phases shipped (Phase 16
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
