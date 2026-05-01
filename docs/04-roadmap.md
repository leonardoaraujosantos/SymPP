# 04 — Roadmap

16 phases. Each phase has: deliverables, SymPy source to study, dependencies, and a ballpark effort estimate (solo developer-weeks). Phases marked **PARALLELIZABLE** can run concurrently with neighbours after their dependency is met.

Total estimate: **~14–20 months solo**, **~6–9 months with 2–3 engineers**.

The two long poles are **Phase 4 (polynomials)** and **Phase 7 (integration)**. Aggressive scope cuts there compress the timeline materially.

---

## Phase 0 — Foundation & oracle harness (2–3 weeks)

**Deliverables**
- Repo bootstrap: CMake, directory layout per [02-architecture](02-architecture.md), `.gitignore`, `.editorconfig`, license file (BSD-3).
- Dependencies wired: GMP, MPFR, Catch2, pybind11, fmtlib (or std::format), CLI11 (for example/repl).
- CI: GitHub Actions matrix — Linux/macOS/Windows × {gcc, clang, msvc} × {Debug, Release} × {ASan, UBSan, plain}.
- clang-tidy + clang-format configs.
- **SymPy oracle harness** — Python subprocess wrapper called from C++ test runner. See [05-validation-strategy](05-validation-strategy.md).
- Empty `Basic`/`Expr` skeleton compiles and links.
- "Hello World" test: construct `Symbol("x")`, print it, validate against SymPy.

**Exit criteria**: `make test` runs Catch2 + oracle harness end-to-end on a trivial expression.

---

## Phase 1 — Core expression tree (4–5 weeks)

**Deliverables**
- `Basic`, `Expr` alias, `Atomic`, `Symbol`, `Dummy`, `Wild`.
- `Number` tower: `Integer` (GMP), `Rational` (GMP), `Float` (MPFR), `Complex`. Coercion rules.
- Singletons: `Zero`, `One`, `NegativeOne`, `Half`, `Pi`, `E`, `EulerGamma`, `Catalan`, `GoldenRatio`, `Infinity`, `NegativeInfinity`, `ComplexInfinity`, `NaN`, `ImaginaryUnit`.
- `Add`, `Mul`, `Pow` with canonical-form auto-evaluation (flatten, sort, combine numerics, basic power rules).
- Hash-consing cache.
- Operator overloads (`+`, `-`, `*`, `/`, unary `-`, `pow`).
- `subs`, `xreplace`, `has`, `free_symbols`, `args`.
- `expand` (basic: distribute Mul over Add; integer-power expansion).
- Structural equality, hash, `cmp`-style stable ordering.
- Basic `str()` printer.
- `evalf(precision)` for numeric expressions (mpmath-equivalent via MPFR).

**SymPy reference**: `core/{basic,expr,symbol,sympify,numbers,add,mul,power,function,operations,evalf,exprtools}.py`.

**Tests vs oracle**: every operation cross-checked against `sympy.sympify(s).<op>()`.

**Exit criteria**: 200+ oracle tests passing covering arithmetic, substitution, basic expand, evalf to 50 digits.

---

## Phase 2 — Assumptions (3–4 weeks)

**Deliverables**
- Old assumptions system: `Symbol("x", real=true, positive=true)` → consistent fact propagation.
- The assumption ontology table (real, integer, rational, positive, negative, nonnegative, nonpositive, zero, nonzero, even, odd, prime, composite, finite, infinite, complex, imaginary, hermitian, antihermitian, irrational, algebraic, transcendental, commutative).
- `is_real`, `is_positive`, etc. queries on every `Expr`.
- Propagation: `Mul(positive, positive).is_positive == true`; `Add(real, real).is_real == true`.
- `assuming` context manager equivalent.
- `refine(expr)` to apply current assumptions.
- **Deferred**: SAT-based `ask` system (port post-v1).

**SymPy reference**: ★`core/assumptions.py`, `assumptions/{assume,facts,ask,refine}.py`. Skip `sathandlers.py`/`satask.py`/`cnf.py`/`lra_satask.py` for v1.

**Exit criteria**: oracle tests on 50 assumption queries; demonstrable propagation through compound expressions.

---

## Phase 3 — Elementary & special functions (5–6 weeks) — PARALLELIZABLE with Phase 4

**Deliverables**
- All elementary functions per [03-feature-mapping](03-feature-mapping.md) Chapter 6 / Elementary.
- All special functions per Chapter 6 / Special.
- All combinatorial functions.
- Each function: auto-evaluation rules, `evalf` via MPFR/mpmath, derivative, `_rewrite_as_*` for major equivalent forms, series expansion.
- `Derivative` and `Integral` as unevaluated forms.

**Approach**: highly mechanical. Build a code-generation script that scaffolds new function classes from a YAML spec (name, args, eval rules, derivative, rewrites). Manual implementation of evaluator bodies. ~80 functions × ~200 LoC each.

**SymPy reference**: `functions/elementary/*.py`, `functions/special/*.py`, `functions/combinatorial/*.py`.

**Exit criteria**: all functions in oracle test suite; numeric values match SymPy to 30 digits.

---

## Phase 4 — Polynomials (8–10 weeks) — LONG POLE

**Deliverables**
- `Poly` class (sparse, distributed monomial form).
- Domains: ℤ, ℚ, ℝ (MPFR), ℂ, ℤ_p, ℚ_alg.
- DUP/DMP arithmetic: add, sub, mul, div, quo, rem, gcd_ex.
- GCD: subresultant + heuristic + modular over ℤ.
- Square-free factorization.
- **Factorization over ℤ**: Berlekamp-Zassenhaus + Hensel lifting (the hard core; ~3 weeks alone).
- **Factorization over ℚ**: clearing denominators + ℤ.
- Multivariate factorization: Wang's algorithm.
- Resultants, discriminants, subresultants.
- Roots: rational roots, low-degree closed-form (≤4), numeric (mpsolve-style isolation), `RootOf` lazy form.
- `apart` (partial fractions).
- `cancel`, `together`.
- `horner`.
- **Deferred to Phase 4b**: Gröbner bases (`polys/groebnertools.py`, `fglmtools.py`) — port if Phase 10 polynomial systems need them.

**SymPy reference**: ★`polys/{polytools,densebasic,densearith,densetools,euclidtools,heuristicgcd,modulargcd,factortools,sqfreetools,partfrac,polyfuncs,polyroots,rootisolation,rootoftools,multivariate_resultants,subresultants_qq_zz}.py` and all of `polys/domains/`.

**Exit criteria**: oracle tests for: factor of polynomials of degree ≤20 with integer/rational coeffs, partial fractions, multivariate GCD on bivariate cases.

---

## Phase 5 — Simplification (5–6 weeks)

**Deliverables**
- `simplify(expr)` orchestrator with strategy switching.
- `collect`, `together`, `cancel` (uses Phase 4).
- `powsimp`, `expand_power_*`.
- `trigsimp` via Fu et al. table-driven rules + classical reduction.
- `radsimp`, `sqrtdenest`.
- `partfrac` (re-export from Phase 4).
- `combsimp`, `gammasimp`.
- `hyperexpand`.
- `cse` (common subexpression elimination).
- `nsimplify` (numeric → exact symbolic guess).

**SymPy reference**: ★`simplify/simplify.py`, plus `radsimp.py`, `trigsimp.py`, `fu.py`, `powsimp.py`, `sqrtdenest.py`, `combsimp.py`, `gammasimp.py`, `hyperexpand.py`, `cse_main.py`.

**Exit criteria**: pass SymPy's `simplify` test suite at ≥80% on shared cases; oracle tests on 100+ standard simplification patterns.

---

## Phase 6 — Calculus (5–7 weeks)

**Deliverables**
- `diff` (already partial from Phase 1; finalize chain rule, partials, `Derivative` propagation).
- `series` — Taylor & Laurent up to user-specified order; `O(x^n)` (Big-O).
- `limit` — port Gruntz algorithm. **Hard**; budget 2 weeks.
- `summation` — closed-form via Gosper, hypergeometric, telescoping.
- `pade` — Padé approximant.
- Functional derivatives — Euler-Lagrange.
- `minimum`, `maximum`, `stationary_points`, `is_increasing`, `is_decreasing`.
- Asymptote / inflection helpers.

**SymPy reference**: ★`series/{limits,gruntz,series,formal,approximants,order,residues}.py`, `concrete/{summations,gosper,expr_with_intlimits}.py`, `calculus/{euler,util}.py`.

**Exit criteria**: 200+ oracle tests including Gruntz limits (Wester test set), classical Taylor expansions, closed-form sums.

---

## Phase 7 — Integration (8–10 weeks) — LONG POLE

**Deliverables (in port order)**
1. Polynomial integration — trivial.
2. Rational function integration — Lazard-Rioboo-Trager. ~1 week.
3. **`manualintegrate`** — rule-based table lookup (powers, trig, exp, log substitution patterns). Covers ~70% of textbook integrals. ~2 weeks. **Port this first** for early user value.
4. Trigonometric reductions.
5. Risch heuristic (`heurisch`) — practical incomplete Risch. ~2 weeks.
6. Meijer G-function method (`meijerint`) — handles most special-function integrals; needed by Laplace transform. ~2 weeks.
7. Definite integration: Newton-Leibniz orchestration; use Meijer G when antiderivative not closed.
8. Numeric quadrature (`vpaintegral`) — wraps mpmath-equivalent.
9. **Deferred**: full Risch transcendental (`risch.py`/`prde.py`/`rde.py`) — research-level; port post-v1 if needed.

**SymPy reference**: ★`integrals/{integrals,manualintegrate,heurisch,meijerint,rationaltools,trigonometry,risch,prde,rde,quadrature}.py`.

**Exit criteria**: oracle tests on the standard CAS integration test set (Wester + SymPy's own); document which classes pass/fail.

---

## Phase 8 — Transforms (3–4 weeks)

**Deliverables**
- `laplace_transform`, `inverse_laplace_transform` (uses Phase 7 Meijer G).
- `fourier_transform`, `inverse_fourier_transform`.
- `mellin_transform`, `inverse_mellin_transform`.
- `z_transform`, `inverse_z_transform`.
- (Optional) sine/cosine/Hankel transforms.

**SymPy reference**: ★`integrals/laplace.py`, `integrals/transforms.py`.

**Exit criteria**: standard tables of transform pairs all match SymPy.

---

## Phase 9 — Linear algebra (3–4 weeks) — PARALLELIZABLE with Phase 8

**Deliverables**
- Dense `Matrix<Expr>` and sparse variant.
- `det`, `trace`, `inv`, `adjugate`, `transpose`, `conjugate_transpose`.
- `rref`, `nullspace`, `columnspace`, `rowspace`, `rank`.
- `eigenvals`, `eigenvects` (uses Phase 4 polys).
- `jordan_form`, `diagonalize`, `is_diagonalizable`.
- SVD, LU, QR, Cholesky.
- Norms (Frobenius, p-norm).
- `MatrixSymbol` and matrix expression tree (`MatMul`, `MatAdd`, `MatPow`, `Inverse`, `Transpose`, `Trace`, `Determinant`, `BlockMatrix`).
- `jacobian`, `hessian`, `wronskian`, `gradient`.
- Hadamard, Kronecker products.
- Rotation matrices.
- Hilbert / Vandermonde / companion matrix constructors.

**SymPy reference**: ★`matrices/{matrices,dense,sparse,immutable,matrixbase,determinant,inverse,reductions,subspaces,eigen,decompositions,normalforms}.py`, `matrices/expressions/`.

**Exit criteria**: oracle parity on the matrix examples in the MATLAB toolbox (Hilbert, SVD, Jordan, eigenvalues of Laplace operator example).

---

## Phase 10 — Equation solvers (5–6 weeks)

**Deliverables**
- `solveset(expr, var, domain)` — preferred path; returns Set objects.
- `linsolve` — using Phase 9.
- `nonlinsolve` — Gröbner-based for polynomial systems.
- `solve` — legacy heuristic API (thin wrapper for MATLAB compatibility).
- `nsolve` — mpmath-equivalent Newton (uses MPFR).
- Inequalities: `solve_univariate_inequality`, `reduce_inequalities`.
- Recurrences: `rsolve`.
- Diophantine (basic Pythagorean, linear) — optional.
- Set algebra (Interval, FiniteSet, Union, Intersection, Complement, ImageSet, ConditionSet) — built here as solver outputs need it.

**SymPy reference**: ★`solvers/{solveset,solvers,polysys,inequalities,simplex,recurr}.py`, `sets/`.

**Exit criteria**: oracle tests on MATLAB toolbox solver examples; the parametric `solve` with `ReturnConditions` mode equivalent.

---

## Phase 11 — ODE/PDE (6–8 weeks)

**Deliverables**
- `dsolve(eq, func)` for single-equation ODEs:
  - Classifier ports `single.py` patterns: separable, linear 1st-order, exact, Bernoulli, Riccati, homogeneous, Lie symmetry-reducible, hypergeometric, Cauchy-Euler, constant-coefficient linear of any order.
  - Per-class solvers.
- `dsolve` for ODE systems (linear constant-coefficient, autonomous).
- Initial conditions / IVP support.
- `pdsolve` — single-PDE classifier (separable, first-order linear, etc.); coverage matches MATLAB's actual reach.
- DAE manipulation: index reduction, Jacobian extraction. (Numeric DAE solving delegated to consumer.)
- `checksol`, `checkodesol`.

**SymPy reference**: ★`solvers/ode/{ode,single,systems,nonhomogeneous,riccati,hypergeometric,lie_group,subscheck}.py`, `solvers/pde.py`, `solvers/deutils.py`.

**Exit criteria**: oracle tests on the 10 ODE patterns and 5 PDE patterns from the MATLAB User's Guide chapter 3.

---

## Phase 12 — Units (1–2 weeks)

**Deliverables**
- `Quantity` type wrapping `(Expr, Dimension)`.
- `Dimension` algebra (length, mass, time, current, temperature, amount, luminosity).
- SI/CGS/US unit systems with conversion.
- Unit prefixes (kilo, milli, micro, …).
- Standard physical constants.
- Definitions table (length, mass, energy, force, pressure, …).

**SymPy reference**: `physics/units/{quantities,dimensions,prefixes,unitsystem,util}.py`, `physics/units/definitions/`, `physics/units/systems/`.

**Exit criteria**: parity with MATLAB unit examples (temperature conversion, SI/CGS/US conversion, unit-aware physics calculation).

---

## Phase 13 — Code generation (4–5 weeks)

**Deliverables**
- Visitor pattern: `CodePrinter` base, `precedence` table.
- Concrete printers: `CCodePrinter`, `CXXCodePrinter`, `FortranCodePrinter`, `LatexPrinter`, `PrettyPrinter` (ASCII), `OctavePrinter` (= MATLAB-compatible for `matlabFunction`).
- Higher-level codegen: function/module emission with arg lists, return types, includes.
- `CSE` (common subexpression elimination).
- `lambdify` — JIT compile expression to native function via LLVM ORC.
- `autowrap` equivalent — emit + invoke C compiler, dlopen.
- `jacobian`/`hessian` code emission for optimization.

**SymPy reference**: ★`printing/{codeprinter,c,cxx,fortran,latex,octave,str,precedence,pycode}.py`, `printing/pretty/`, `codegen/{algorithms,ast,cnodes,cxxnodes,fnodes,cutils,futils,rewriting}.py`, `utilities/{lambdify,autowrap}.py`, `simplify/cse_main.py`.

**Exit criteria**: round-trip — symbolic expression → C code → compile → execute → numeric match SymPy's `lambdify` output. Same for Fortran, LaTeX rendering, Octave (validated by running through MATLAB if available).

---

## Phase 14 — Plotting bridge (1 week)

**Deliverables**
- `SampledFunction` — sample `Expr` over a domain into a buffer.
- Adaptors: gnuplot stream, matplotlib-cpp, raw CSV/JSON export.
- 1D, 2D, 3D, contour, implicit, parametric variants.

**Exit criteria**: example consumer renders the MATLAB toolbox pendulum animation using SymPP + matplotlib-cpp.

---

## Phase 15 — Parser & MATLAB facade (2 weeks)

**Deliverables**
- `parse("x^2 + sin(x) == 0")` — recursive descent or PEG grammar; returns `Expr`.
- MATLAB-compatible top-level functions in `sympp::matlab::` namespace: `syms`, `sym`, `solve`, `dsolve`, `int`, `diff`, `simplify`, `subs`, `factor`, `expand`, `taylor`, `laplace`, `fourier`, `ztrans`, `vpa`, `eval`, etc.
- Each MATLAB function is a thin wrapper over the SymPP core API to ease porting MATLAB scripts.

**SymPy reference**: `parsing/sympy_parser.py` (algorithm only; can't reuse the AST/`eval` approach in C++).

**Exit criteria**: 20-example MATLAB regression suite (selected from the User's Guide) parses and runs.

---

## Phase 16 — Hardening & v1.0 release (4–6 weeks)

**Deliverables**
- Run SymPy's own test suite subset against SymPP via the oracle harness; document pass rates.
- Performance benchmarks vs SymEngine baseline.
- Doxygen / public-header docs complete.
- Tutorial: "Porting a MATLAB script to SymPP" — walk through 5 examples.
- vcpkg port, Conan recipe, CMake `find_package` example.
- Pre-built binaries for Linux/macOS/Windows.
- Tagged v1.0 release.

**Exit criteria**: a clean Linux/macOS/Windows install + sample app + green tests on a fresh machine.

---

## Critical path summary

```
0 → 1 → 2 ─┬─ 3 ─┐
           ├─ 4 ─┼─ 5 ─ 6 ─ 7 ─┬─ 8 ──┐
           │     │              ├─ 9 ──┼─ 10 ─ 11 ─ 12 ─ 13 ─ 14 ─ 15 ─ 16
           └─────┘
```

3 and 4 can run in parallel after 2. 8 and 9 can run in parallel after 7.

## Aggressive-cut v1.0

If solo developer + 9-month deadline, cut:
- Phase 4: skip multivariate factorization, Gröbner. (Saves 3 weeks.)
- Phase 7: ship only manualintegrate + rational. Skip Risch and Meijer G. (Saves 4 weeks; **but Phase 8 transforms become much weaker**.)
- Phase 11: ODE only, no PDE, no DAE. (Saves 2 weeks.)
- Phase 13: C and LaTeX printers only. (Saves 2 weeks.)
- Skip Phase 15 MATLAB facade — users can use the SymPP API directly. (Saves 2 weeks.)

This produces a "useful for engineers but not MATLAB-equivalent" v0.5 in ~8 months.
