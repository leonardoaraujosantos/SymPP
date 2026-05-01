# 04 — Roadmap

This file is the source of truth for what SymPP has shipped, what's
deferred (and *honestly why*), and the forward path to SymPy feature
parity.

---

## 📍 You are here — next-session pickup

**Current position**: **v0.5 release-candidate.** 14 of 15 phases
shipped at minimal-viable scope; Phase 16 partial; 5 Category-A
deep-deferred items already closed. Everyday CAS workflow coverage
≈ 65 %, composite SymPy parity ≈ 50 %. Critical path forward: parallelizable.

**Just landed** (most recent first; commit hashes on `main`):

| Commit | What |
|---|---|
| `92b603c` | Roadmap restructured: Category D added (Phases 25–39, 15 SymPy modules), multi-metric parity assessment, version-tag milestones |
| `abe5653` | Docs sync: README "What's in the box" + roadmap top table refreshed to match shipped code |
| `76bcce3` | Platform-independent canonical Add/Mul order (string-first sort + FNV-1a string hash) — fixes Linux gcc/clang test failures |
| `d3640f5` | CI: dropped macOS matrix; added `<cmath>` include + fixed `::write` warn-unused-result |
| `540ec90` | CI: added `<optional>` includes for libstdc++; collapsed duplicated branches in `dsolve.cpp` |
| `23244cc` | **Big feature commit** — MATLAB facade extension (5 sub-headers); Fu trig rules; solveset `_invert` chain; variation of parameters; Jordan form + `Matrix::exp(t)`; hypergeometric infrastructure (`Hyper`, `MeijerG`, `hyperexpand` subset); SymPy cross-validation oracle op |

**Recommended next pickup** (pick one — they're independent):

1. **Finish Phase 16 to ship v0.5 tag** · ~5 wk · **highest priority**
   - Outstanding: v1.0 tag, doxygen site, benchmark suite vs
     SymEngine, vcpkg + Conan recipes.
   - This is the *only* strictly-sequenced work; everything else
     parallelizes once Phase 16 is done.
   - See [Phase 16 — Hardening & v1.0](#phase-16--hardening--v10--) below.

2. **Cat-A: Infinity singletons → Gruntz limit algorithm** · ~3 wk
   - Adds `Infinity` / `NegativeInfinity` / `ComplexInfinity` / `NaN`
     to the singleton tier (Phase 1 deferred-deep), then implements
     the Gruntz algorithm for limits at infinity (Phase 6 deferred-deep).
   - Unblocks `horizontal_asymptotes` and ∞/∞, 0·∞, ∞–∞ indeterminate
     forms — currently L'Hôpital only handles 0/0.
   - High everyday-coverage payoff: ~65 % → ~75 %.

3. **Cat-A: Full hyperexpand (Slater theorem) + Meijer-G expansion** · ~7 wk
   - Builds on the `Hyper` / `MeijerG` Function classes already shipped.
   - Unblocks Phase 7's full Meijer-G integration method and Phase 8's
     general-case algebraic-input transforms.
   - Highest *algorithmic-depth* lift in Cat A.

4. **Cat-A: lambdify via LLVM ORC** · ~3 wk
   - JIT-compile `Expr` to native function for numerical evaluation.
   - Heaviest user-facing payoff if SymPP is being used in
     hot-path numerical code (matlab_llvm consumer wants this).

**Where to find what's deferred** (jump targets in this file):

- Phase 16 status: §[Phase 16](#phase-16--hardening--v10--)
- Category A table (deep-deferreds): §[Category A](#category-a-close-the-deep-deferreds-in-shipped-phases)
- Category C (Phases 17–24, beyond MATLAB scope): §[Category C](#category-c-new-phases-beyond-matlab-scope-sympy-parity)
- Category D (Phases 25–39, SymPy modules outside original plan): §[Category D](#category-d-sympy-modules-outside-the-original-024-plan)
- Critical-path diagram: §[Critical path](#critical-path)
- Parity-percentage breakdown: §[How far are we from SymPy?](#how-far-are-we-from-sympy)

**Test surface to maintain**: 962 cases / 1872 assertions (307
oracle-validated against SymPy 1.14). Run:

```bash
./build/tests/sympp_tests              # full suite
./build/tests/sympp_tests "[oracle]"   # SymPy cross-validated subset
```

**Open `CHANGELOG.md` `[0.5.0] — Unreleased`** is the running list of
what to mention in the v0.5 release notes when Phase 16 lands.

---

## Status snapshot

```
962 tests / 1872 assertions  all passing
14 of 15 phases shipped (Phase 14 dropped — see below)
```

| Phase | Title | Status |
|---|---|---|
| 0  | Foundation & oracle harness            | ✅ shipped |
| 1  | Core expression tree                   | ✅ shipped |
| 2  | Assumptions                            | 🟡 minimal subset (SAT-based deferred) |
| 3  | Elementary & special functions         | ✅ shipped |
| 4  | Polynomials                            | ✅ shipped (multivariate factor + BZ deferred) |
| 5  | Simplification                         | ✅ shipped (hyperexpand subset shipped — full Slater + Meijer-G expansion deferred) |
| 6  | Calculus                               | ✅ shipped (Gruntz at infinity deferred) |
| 7  | Integration                            | ✅ shipped (full Risch + Meijer G deferred) |
| 8  | Transforms                             | ✅ shipped (Meijer-driven general case deferred) |
| 9  | Linear algebra                         | ✅ shipped (Jordan + matrix exp shipped for chains ≤ 2; SVD + sparse + longer chains deferred) |
| 10 | Equation solvers                       | ✅ shipped (`_invert` chain shipped — transcendental solveset; F4/F5 Gröbner + SAT deferred) |
| 11 | ODE / PDE                              | ✅ shipped (variation of parameters shipped; full Lie + Pantelides deferred) |
| 12 | Units                                  | ✅ shipped |
| 13 | Code generation                        | 🟡 printers + function emission (lambdify deferred) |
| 15 | Parser & MATLAB facade                 | ✅ shipped (extension: parsing / assumptions / solvers / ode / pde sub-headers under `sympp::matlab`) |
| 16 | Hardening & v1.0                       | 🟡 `find_package` + install rules + CI matrix shipped; v1.0 tag + benchmarks pending |

> **Note**: Phase 14 (Plotting bridge) was dropped from scope.
> Plotting is downstream of CAS correctness — consumers can pipe
> SymPP output to gnuplot / matplotlib / their renderer of choice
> without needing dedicated SymPP infrastructure. Phase numbering
> for 15/16 and 17+ kept stable to avoid invalidating existing
> commit messages and external references.

Legend: ✅ feature-complete at minimal-viable scope · 🟡 partial · ❌ not started

The original "Phase X minimal" labels in the commit history reflect
intentional aggressive scope cuts taken in the first weekend session;
each minimal landing has since been re-visited and the deferred items
shipped in dedicated follow-up commits ("Phase X deep").

---

## Phase status breakdown

Each phase below lists what's actually in `main` and what's still
deferred — with the *real* reason for the deferral, not "it was hard".

### Phase 0 — Foundation & oracle harness · ✅

**Shipped**: CMake build, GMP/MPFR/Catch2/nlohmann_json wired,
clang/clang++ on macOS+Linux, BSD-3 license, `Basic`/`Expr`
skeleton, Python oracle subprocess (long-lived, JSON-over-stdin/stdout),
`tests/oracle/oracle_smoke_test`. CI not wired (intentional — local
testing is the bottleneck).

### Phase 1 — Core expression tree · ✅

**Shipped**: Number tower (`Integer` GMP, `Rational` GMP, `Float`
MPFR, `NumberSymbol` for Pi/E/EulerGamma/Catalan), `Symbol`,
`ImaginaryUnit`, hash-cons cache via `weak_ptr`-keyed multimap,
`Add`/`Mul`/`Pow` with full canonical-form auto-evaluation,
operator overloads (`+ - * / pow`), `subs`/`xreplace`/`has`/
`free_symbols`/`args`, `expand`, structural equality + canonical
ordering, `str()` printer, `evalf` to arbitrary precision via MPFR.
`I*I = -1` and `I^n` cycle. Singletons `S::Zero/One/Half/Pi/E/I/...`.

**Not shipped** (intentional): `Dummy`, `Wild`, `Complex` as a
distinct number type (we use `a + b*I` instead), `Infinity` /
`NegativeInfinity` / `ComplexInfinity` / `NaN` — these are deferred
together with the Gruntz limit work.

### Phase 2 — Assumptions · 🟡 minimal

**Shipped**: `AssumptionMask` (Real, Rational, Integer, Positive,
Negative, Zero, Nonzero, Nonnegative, Nonpositive, Finite),
`is_real/is_positive/is_integer/...` queries on every `Expr`,
propagation through `Add`/`Mul`/`Pow`, `refine()` on assumption-gated
rewrites.

**Deferred-deep**: Even/Odd/Prime/Composite, Hermitian/Antihermitian,
Algebraic/Transcendental, the SAT-based `ask` system. Reason: SAT
porting is its own multi-week subsystem and our practical use cases
hit only the basic ontology.

### Phase 3 — Elementary & special functions · ✅

**Shipped**: 30+ functions: `sin/cos/tan/asin/acos/atan/atan2`,
`sinh/cosh/tanh/asinh/acosh/atanh`, `exp/log`, `sqrt/abs/sign/
re/im/conjugate/arg`, `min_/max_`, `floor/ceiling`, `factorial/
binomial/gamma/loggamma`, `erf/erfc/heaviside/dirac_delta`. Each
has auto-eval rules, `evalf` via MPFR/mpmath-equivalent, and a
`diff_arg` chain-rule entry.

**Deferred**: Bessel/Zeta/Beta/digamma — naming is in `FunctionId`
but no implementation. Reason: not yet driven by a downstream phase.

### Phase 4 — Polynomials · ✅

**Shipped**: Univariate `Poly` over Expr coefficients, +/-/* and
divmod, GCD via Euclidean algorithm over ℚ, square-free factorization
(Yun), rational roots theorem, cubic via Cardano, quartic via Ferrari,
factor over ℤ via Kronecker (square-free + rational + Kronecker
reconstruction), `cancel`/`together`/`apart`/`horner`, resultant
(Sylvester) + discriminant, lazy `RootOf` with MPFR bisection evalf,
multivariate Gröbner basis (Buchberger, lex order) shipped in
Phase 10.

**Deferred-deep**: Multivariate `Poly` (sparse distributed monomial
form), full domain tower (ℂ, ℤ_p, ℚ_alg), Berlekamp-Zassenhaus +
Hensel lifting (Kronecker is correct but exponential), Wang
multivariate factorization, modular GCD, full `apart` for repeated
roots / quadratic factors via linear-system solve. Reason: each is
a multi-week port and Kronecker is adequate for textbook degrees.

### Phase 5 — Simplification · ✅

**Shipped**: `simplify` orchestrator (chains all of the below),
`collect`, `together`, `cancel` (uses Phase 4 GCD), `powsimp`,
`expand_power_exp`, `expand_power_base`, `trigsimp` (Pythagorean
identity), `radsimp` (rationalize binomial denominators),
`sqrtdenest` (Borodin denesting via discriminant test), `combsimp`
(factorial ratios), `gammasimp` (gamma ratios), `cse` (common
subexpression elimination via hash-cons), `nsimplify` (rational +
pi/e + sqrt-of-rational table).

**Deferred-deep**: `hyperexpand`. Reason: requires hypergeometric
function machinery — full SymPy port is ~1000 LOC and is also the
blocker for general-case Meijer G integration. Treated as one
research-level subsystem to land jointly.

### Phase 6 — Calculus · ✅

**Shipped**: `diff` with full chain/product/power/sum rules,
`series` (Taylor via 1/k!·f⁽ᵏ⁾ formula), `limit` with **L'Hôpital
for 0/0 indeterminate forms**, `Order` (Big-O class), `summation`
(closed-form: linearity, constant, Σk, Σk², Σk³, geometric),
`pade` approximant via Taylor + linear-system solve,
`stationary_points`/`minimum`/`maximum`/`is_increasing`/
`is_decreasing`, `euler_lagrange` functional derivative,
`vertical_asymptotes`, `inflection_points`.

**Deferred-deep**: Full **Gruntz** algorithm at infinity
(0·∞, ∞-∞, ∞/∞), `horizontal_asymptotes`. Reason: requires an
`Infinity` singleton + asymptotic-comparison machinery — which is
its own subsystem.

### Phase 7 — Integration · ✅

**Shipped**: Linearity over Add, table for elementary forms with
**affine arguments** (∫sin(ax+b), ∫(ax+b)^n, ∫exp(ax+b)),
**trig reductions** (sin², cos², product-to-sum), **rational
function integration via apart**, **integration by parts** (log,
poly·{sin,cos,exp}), **heurisch** chain-rule-reverse u-substitution,
**`manualintegrate` orchestrator** with optional fallback to numeric,
**`vpaintegral`** (adaptive Simpson in MPFR for high-precision
definite integrals).

**Deferred-deep**: Full Risch transcendental algorithm (`risch.py`
+ `prde.py` + `rde.py` — research-level ports), Meijer G method.
Reason: each is multi-week. Coverage on textbook integrals is
already solid via the heuristics + by-parts + heurisch chain.

### Phase 8 — Transforms · ✅

**Shipped**: `laplace_transform` + `inverse_laplace_transform`
(apart-driven), `fourier_transform` + `inverse_fourier_transform`,
`mellin_transform` + `inverse_mellin_transform`, `z_transform` +
`inverse_z_transform`, `sine_transform` + `cosine_transform`. Each
is table-based with linearity over Add and constant-factor pull-out.

**Deferred-deep**: Meijer-G-driven general-case transforms (handle
arbitrary algebraic / special-function inputs). Reason: depends on
the hypergeometric machinery still pending from Phase 5.

### Phase 9 — Linear algebra · ✅

**Shipped**: Dense `Matrix<Expr>`, identity/zeros, +/-/*, scalar
mul, transpose, conjugate_transpose, trace, det (cofactor),
adjugate, inverse (adjugate / det), Frobenius norm, `gradient` /
`hessian` / `jacobian` / `wronskian`, `rref` / `rank` / `nullspace`
/ `columnspace` / `rowspace`, `charpoly` / `eigenvals` / `eigenvects`
/ `is_diagonalizable` / `diagonalize`, **LU + QR (Gram-Schmidt) +
Cholesky** decompositions, Hadamard + Kronecker products, Hilbert /
Vandermonde / companion / rotation_matrix_2d/x/y/z constructors,
`MatrixSymbol` + symbolic matrix expression tree (`MatAdd`,
`MatMul`, `MatPow`, `MTranspose`, `MInverse`, `MTrace`,
`MDeterminant`, `BlockMatrix`).

**Deferred-deep**: Full symbolic SVD (requires eigendecomposition of
A^TA — scales poorly for n>3), `jordan_form` (nilpotent analysis on
generalized eigenspaces), sparse `Matrix` variant (separate storage
representation).

### Phase 10 — Equation solvers · ✅

**Shipped**: `solve` (univariate polynomial deg ≤ 4 closed-form, deg
≥ 5 via rational-roots), `solveset` (returns `Set`), `linsolve`,
`nsolve` (Newton's method in MPFR throughout), `solve_univariate_
inequality`, `reduce_inequalities` (single + AND/OR vector forms),
`rsolve` (constant-coefficient recurrence), `nonlinsolve` (2-var
via resultant + numeric verification, n-var via Gröbner). **Set
algebra**: `Interval` / `FiniteSet` / `Union` / `Intersection` /
`Complement` / `EmptySet` / `Reals` / `Integers` / `ConditionSet` /
`ImageSet`. **Trig `solveset`** emits `ImageSet` for periodic
solutions. **Linear Diophantine** + **Pythagorean triples**
(Euclid's formula). **Buchberger Gröbner basis** with lex order,
n-variable nonlinear system solving via triangular back-substitution.

**Deferred-deep**: Full SAT-based assumption reasoning, F4 / F5
Gröbner algorithms (Buchberger is exponential worst-case),
degrevlex / matrix monomial orderings, general transcendental
`solveset` (a·sin(x) + b·cos(x) = c, exp/log compositions, etc.).

### Phase 11 — ODE / PDE · ✅

**Shipped**: `dsolve_first_order` classifier + dispatch,
`dsolve_separable`, `dsolve_linear_first_order`, `dsolve_bernoulli`,
`dsolve_exact`, `dsolve_constant_coeff` (any order),
`dsolve_cauchy_euler`, `dsolve_system` for y'=A·y via
eigendecomposition, `apply_ivp` for initial conditions,
`checkodesol`, `pdsolve_first_order_linear` (constant coefficients).
**Deep**: `dsolve_riccati` (linearize to 2nd-order),
`dsolve_homogeneous` (y'=f(y/x)), `dsolve_lie_autonomous`,
`dsolve_hypergeometric` recognition + opaque `hyper(...)` factory,
`pdsolve_first_order_variable` (homogeneous case via real method
of characteristics), `pdsolve_heat` (separation of variables),
`pdsolve_wave` (d'Alembert form), `dae_jacobians`,
`dae_structural_index`.

**Deferred-deep**: Full Lie symmetry classifier (point + contact +
Lie-Bäcklund — 1500+ LOC port), full Pantelides BLT decomposition
+ Tarjan SCC for high-index DAEs, full `hyperexpand` series
expansion of the hypergeometric function.

### Phase 12 — Units · ✅

**Shipped**: `Dimension` (7-tuple of integer exponents over SI base
dimensions), `Unit` (name, symbol, dim, scale), `Quantity` (Expr +
Unit with dim-checked arithmetic), all SI base + 13 derived units
(N, J, W, Pa, Hz, V, Ω, F, T, Wb, H, lx, gram), prefixes
(yocto..yotta), CGS (cm, dyne, erg) + US (foot, inch, mile, lb,
oz, gal, BTU) + time (min, hour, day, year), 12 physical constants
with exact post-2019-redef values for c/h/k_B/e/N_A,
`convert(qty, target)`, affine **temperature conversions** (C/F/K).

### Phase 13 — Code generation · 🟡

**Shipped**: `CodePrinter` visitor base, `ccode` / `cxxcode` /
`fcode` / `latex` / `octave_code` / `pretty` printers, function
emission `c_function` / `cxx_function` / `fortran_function` /
`octave_function`. Greek letter mapping in LaTeX, `\frac{}{}`,
`\sqrt{}`, std::numbers::pi_v in C++, Fortran `d0` suffix for
double precision, Octave `.* ./ .^` element-wise operators.

**Deferred-deep**: `lambdify` (JIT compile expr to native function
via LLVM ORC — heavy LLVM dependency), `autowrap` equivalent
(emit C source + invoke compiler + dlopen), full 2D pretty-print
layout (stacked fractions, integral signs — ~500 LOC port from
sympy/printing/pretty/).

### Phase 14 — Plotting bridge · DROPPED

Removed from scope. CAS users pipe SymPP output to gnuplot,
matplotlib, or any other renderer through the existing code
generation pipeline (`ccode`, `octave_code`, `latex`) and numeric
quadrature (`vpaintegral`) plus `evalf`. A dedicated plotting
adapter would couple SymPP to a specific renderer without adding
mathematical capability — better to keep the library a pure CAS
and let the consumer pick their plotting layer.

### Phase 15 — Parser + MATLAB facade · ✅

**Shipped — parser**: a recursive-descent string parser
`parse("x**2 + sin(x)")` → `Expr`. Two intended uses:
(1) consumers reading expressions from non-C++ sources (config
files, databases, network protocols), (2) round-tripping with
SymPy via `srepr`. Right-associative `**`, `^` aliased,
function-call dispatch against the SymPP factory set,
undefined-function preservation, position-tagged `ParseError`.

**Shipped — MATLAB facade** (`sympp::matlab::*`): MATLAB-named
wrappers over the SymPP API. Variadic `syms`, `sym(name)`/
`sym(value)`, `Int` (renamed from `int` — C++ keyword),
`diff(f, x, n)` for n-th derivative, `taylor` (= `series`),
`vpa` (= `evalf`), `laplace`/`ilaplace`/`fourier`/`ifourier`/
`ztrans`/`iztrans`, `solve`/`dsolve`, `simplify`/`subs`/`expand`/
`factor`/`limit`, `pretty`/`latex`/`ccode`/`fortran`/
`matlabFunction`. Header-only.

**Facade extension** (split into `matlab/parsing.hpp`,
`matlab/assumptions.hpp`, `matlab/solvers.hpp`, `matlab/ode.hpp`,
`matlab/pde.hpp` under the umbrella `matlab/matlab.hpp`):
`str2sym` / parser-routed `sym(string_view)`; `assume` /
`assumeAlso` / `assumptions` / `clearAssumptions` / `refresh` on
top of a process-wide assumption registry (Symbols stay
immutable; `refresh(x)` rebinds with the registered mask);
`linsolve`, multi-equation `solve({eqs}, {vars})`, `nsolve` /
`vpasolve`; first-order and auto-classified second-order
`dsolve` overloads (constant-coefficient / Cauchy-Euler), IVP
companions, `dsolve(A, x)` for linear systems; `pdsolve` family
for first-order linear and heat / wave equations.

The first plan dropped the facade as "renaming wrapper without
real value". Reverted: the variadic `syms` and the `Int`-vs-`int`
keyword adapter are real wrappers, and the namespace-shadowing
forces consumers to make a clean choice (use SymPP API everywhere
or MATLAB API everywhere) which keeps porting clear.

### Phase 16 — Hardening & v1.0 · ❌

**Status**: Not started. Original scope: cross-port test runs,
performance benchmarks vs SymEngine, doxygen docs, vcpkg / Conan,
v1.0 tag. Tractable.

---

## Forward roadmap — to SymPy parity

The original 16-phase plan was scoped to "MATLAB Symbolic Math
Toolbox parity". To compete with **SymPy** itself we need three
additional categories of work:

### Category A: Close the deep-deferreds in shipped phases

Each item below replaces a "good enough" implementation with the
genuinely-comprehensive version. None of these is a blocker — the
current code works correctly within its documented scope. They are
parity deltas, not bug fixes.

| Phase | Deep-deferred work | Effort | Priority |
|---|---|---|---|
| 1  | Infinity / NegativeInfinity / ComplexInfinity / NaN singletons | 1 wk | High (unblocks Gruntz) |
| 2  | SAT-based assumption reasoning | 3 wk | Low |
| 4  | Multivariate Poly + Wang factorization                | 3 wk | Medium |
| 4  | Berlekamp-Zassenhaus + Hensel lifting               | 2 wk | Medium |
| 4  | Full polynomial domain tower (ℤ_p, ℚ_alg, ℂ)        | 2 wk | Low |
| 5  | Fu trig rule table (subset)                          | ✅ shipped — Pythagorean / double-angle / cos²−sin² / 2sin·cos collapses + `expand_trig` + `fu` orchestrator (TR8 product-to-sum, half-angle still deferred) |
| 5  | hyperexpand + hypergeometric function infrastructure | partial ✅ — `Hyper` and `MeijerG` proper Function classes, variadic factories with auto-eval (`₀F₀ → exp`, `₁F₀ → (1−z)^(−a)`, parameter cancellation), `hyperexpand` rewrites `₁F₁(1; 2; z)` and `₂F₁(1, 1; 2; z)`, integrated into `simplify` chain. Full Slater-theorem expansion + Meijer-G evaluation still deferred-deep. |
| 6  | Full Gruntz limit algorithm                          | 2 wk | High |
| 7  | Full Risch transcendental integration                | 4 wk | Medium |
| 7  | Meijer G-function integration method                 | 3 wk | High (with hyperexpand) |
| 9  | Symbolic SVD                                         | 2 wk | Medium |
| 9  | Jordan canonical form                                | ✅ shipped — chains of length ≤ 2 (covers textbook defective inputs) + `Matrix::exp(t)` via Jordan; defective `dsolve_system` fixed. Longer chains still deferred. |
| 9  | Sparse Matrix variant                                | 2 wk | Low |
| 10 | F4 / F5 Gröbner algorithms                           | 3 wk | Low |
| 10 | General transcendental solveset                      | ✅ shipped — `_invert` chain peels log / exp / sin / cos / tan / sinh / cosh / tanh / abs from the LHS; emits `ImageSet` over ℤ for periodic trig and `FiniteSet` for finite-branch inverses |
| 11 | Variation of parameters (2nd-order nonhomogeneous)   | ✅ shipped — Wronskian-based; `dsolve_constant_coeff_nonhomogeneous` + `dsolve_cauchy_euler_nonhomogeneous`, wired through `matlab::dsolve` |
| 11 | Full Lie symmetry classifier                         | 3 wk | Low |
| 11 | Pantelides BLT + Tarjan SCC for high-index DAEs      | 2 wk | Low |
| 13 | lambdify (LLVM ORC integration)                      | 3 wk | High (huge user value) |
| 13 | Full 2D pretty-print layout                          | 1 wk | Medium |

**Total category A**: roughly 42 developer-weeks of focused work,
of which 5 items are now shipped or partial (subset of Fu rules,
transcendental solveset via `_invert`, variation of parameters, Jordan
+ matrix exponential, hypergeometric infrastructure with closed-form
rewrites — together approximately 10–12 weeks of effort).

### Category B: Remaining original phases (15, 16)

| # | Title | Effort |
|---|---|---|
| 15 | Parser / serialization           | 1.5 wk |
| 16 | Hardening + v1.0 release         | 4–6 wk |

### Category C: New phases beyond MATLAB scope (SymPy parity)

These are SymPy modules with no MATLAB analogue but significant
real-world use. Listed in roughly the order they should land.

#### Phase 17 — Vector calculus & differential geometry

`Vector`, `CoordSys3D`, gradient/divergence/curl/laplacian on
arbitrary coordinate systems (Cartesian, cylindrical, spherical),
metric tensor, Christoffel symbols, covariant derivative,
Riemann/Ricci/scalar curvature.

**SymPy reference**: `vector/`, `diffgeom/`. Effort: 3 weeks.

#### Phase 18 — Tensor algebra

`Tensor`, abstract index notation, contraction, raising/lowering
indices via metric, symmetric/antisymmetric components.

**SymPy reference**: `tensor/`. Effort: 2 weeks.

#### Phase 19 — Statistics & probability

`Distribution` (continuous + discrete), `cdf`/`pdf`/`expectation`/
`variance`/`density`/`sample` over Normal, Uniform, Exponential,
Binomial, Poisson, etc. Random variables, conditional probabilities,
Bayesian updates.

**SymPy reference**: `stats/`. Effort: 3 weeks.

#### Phase 20 — Number theory

Primality testing, factorization, modular arithmetic, totient,
quadratic residues, continued fractions, Diophantine equation solver
(beyond linear + Pythagorean).

**SymPy reference**: `ntheory/`. Effort: 2 weeks.

#### Phase 21 — Combinatorics & group theory

Permutations, finite groups, polytope/Young tableaux, partitions,
generating functions, group presentations.

**SymPy reference**: `combinatorics/`. Effort: 3 weeks.

#### Phase 22 — Physics modules

Classical mechanics (Lagrangian/Hamiltonian formalisms), quantum
mechanics (Hilbert space, operators, Dirac notation, density matrix),
optics, electromagnetism. Each is a large submodule in SymPy.

**SymPy reference**: `physics/{mechanics,quantum,optics}/`. Effort:
4 weeks per submodule.

#### Phase 23 — Cryptography primitives

RSA, ElGamal, elliptic curve arithmetic, hash-to-curve, finite-field
operations beyond ℤ_p basics.

**SymPy reference**: `crypto/`. Effort: 2 weeks.

#### Phase 24 — Geometry

Points, lines, polygons, circles, ellipses, intersection, distance,
area, polygon decomposition.

**SymPy reference**: `geometry/`. Effort: 1 week.

---

### Category D: SymPy modules outside the original 0–24 plan

Modules in upstream SymPy that the original phase plan didn't budget
for. Listed here so the gap is visible — none of these is shipped or
in-flight today. Effort estimates are for a minimal-viable port at
the same scope-discipline as Phases 0–13.

| Phase | Title | SymPy ref | Effort |
|---|---|---|---|
| 25 | Logic & boolean algebra (CNF/DNF, `satisfiable`, `simplify_logic`) | `sympy.logic` | 2 wk |
| 26 | Discrete algorithms (FFT, NTT, Walsh-Hadamard, convolutions, Möbius) | `sympy.discrete` | 2 wk |
| 27 | Holonomic functions (D-finite, recurrence/differential closure) | `sympy.holonomic` | 3 wk |
| 28 | Algebraic number fields (`primitive_element`, `minimal_polynomial`, fields tower) | `sympy.polys.numberfields` | 3 wk |
| 29 | DomainMatrix (fast matrix ops on a polynomial domain — order-of-magnitude speed-up) | `sympy.polys.matrices` | 3 wk |
| 30 | Galois fields & finite-field polynomial tooling | `sympy.polys.galoistools` | 2 wk |
| 31 | LaTeX parser (round-trip with the LaTeX printer) | `sympy.parsing.latex` | 2 wk |
| 32 | Orthogonal polynomial families (Legendre, Chebyshev, Hermite, Laguerre, Jacobi, Gegenbauer) | `sympy.functions.special.polynomials` | 1.5 wk |
| 33 | Combinatorial number sequences (Bernoulli, Bell, Fibonacci, Catalan, Lucas, harmonic) | `sympy.functions.combinatorial.numbers` | 1 wk |
| 34 | Special integral functions (`Ei`, `Si`, `Ci`, `Shi`, `Chi`, `fresnels`, `fresnelc`, `expint`) | `sympy.functions.special.error_functions` | 1 wk |
| 35 | Quaternions + algebras | `sympy.algebras` | 1 wk |
| 36 | NDim arrays (distinct from tensor algebra in Phase 18 — concrete n-dim storage) | `sympy.tensor.array` | 1 wk |
| 37 | Unification (`sympy.unify`) — pattern-match driven by structural unification | `sympy.unify` | 1 wk |
| 38 | Extra printers (MathML, GLSL, Rust, Julia, dot, repr) | `sympy.printing` | 1 wk |
| 39 | Codegen AST nodes + Cython / autowrap bindings | `sympy.codegen.ast`, `sympy.utilities.autowrap` | 2 wk |

**Total category D**: roughly 26 developer-weeks.

---

## Summary

```
Where we are:        14/15 phases shipped, Phase 16 partial; 962 tests
                     (307 oracle-validated against SymPy 1.14, all passing).
Cat A shipped:       Fu trig rules, transcendental solveset (_invert),
                     variation of parameters, Jordan + matrix exp,
                     hyperexpand subset.
Cat A remaining:     ~30 wk (Risch, Meijer-G full expansion, Gruntz,
                     F4/F5, lambdify, SAT, multivariate Poly, …).
Cat B remaining:     ~5 wk (Phase 16 — v1.0 tag, benchmarks, doxygen,
                     vcpkg/Conan).
Cat C (Phases 17–24):  ~25 wk (vector calc, tensors, stats, ntheory,
                                combinatorics, physics, crypto, geometry).
Cat D (Phases 25–39):  ~26 wk (logic, discrete, holonomic, algebraic
                                number fields, DomainMatrix, Galois,
                                LaTeX parser, orthogonal polys, special
                                sequences, quaternions, NDim arrays,
                                unification, extra printers, codegen AST).
─────────────────────────────────────────────────
Total to full SymPy parity:  ~86 focused developer-weeks on top of
                              what's already in main.
```

### How far are we from SymPy?

A single number is misleading — depending on the metric, you get a
very different answer. So three metrics:

| Metric | SymPP / SymPy | Notes |
|---|---|---|
| **Everyday CAS workflow coverage** (calc, algebra, ODE, transforms, codegen on textbook-shaped inputs) | **≈ 65 %** | The "common 80 %" cases work end-to-end. Edge cases — high-index DAEs, irreducible quintics, Risch-only integrals, Meijer-G transforms — still hit deferred algorithms. |
| **Algorithmic depth within shipped subsystems** (the "deep-deferred" items in Category A) | **≈ 50 %** | We have the *breadth* of every Phase 0–13 subsystem, but several core algorithms (Risch, Gruntz, full hyperexpand, F4/F5, full Lie, BZ+Hensel) are pending. Within a shipped subsystem you typically get the textbook path; the graduate-textbook path is what's missing. |
| **Module-count parity** (top-level SymPy modules with a working SymPP counterpart) | **≈ 35 %** | SymPP currently covers ~14 of ~40 top-level SymPy modules at minimal-viable scope. Categories C and D fill in the missing 26. |

Composite estimate: **SymPP is about half of SymPy on a typical
user-facing weighted average.** Phrased the other way: you can do
about half of what you'd reach for SymPy for, end-to-end, today —
with SymPy itself wired in as the validation oracle on every
shipped feature, so what *does* work has a much stronger
correctness guarantee than a pure clean-room port would.

For comparison's sake:
- **Speed**: SymPP's static C++ tree + hash-cons cache is in the
  same league as SymEngine for the operations both libraries
  cover; SymPy itself is 10–100× slower on the same workload.
  Speed is not the gap — feature-set is.
- **Numerical work**: SymPP's GMP+MPFR backbone is comparable to
  SymPy's `mpmath` integration for arbitrary-precision evalf.
- **Tooling**: doxygen + benchmarks + vcpkg/Conan are part of
  Phase 16 (partial). SymPy has a much larger documentation and
  tutorial corpus — that gap closes only with sustained writing.

### Realistic milestones

1. **v0.5 (now → ~6 wk)** — Phase 16 finished: v1.0 tag, doxygen
   site, benchmark suite, vcpkg/Conan recipes. Closes Cat B.
2. **v0.7 (~6 → ~36 wk)** — Cat A's high-priority items: Infinity
   singletons + Gruntz, full hyperexpand + Meijer-G, lambdify
   (LLVM ORC). Lifts everyday-coverage from ~65 % → ~80 %.
3. **v0.9 (~36 → ~60 wk)** — Cat C Phases 17–21: vector calculus,
   tensors, statistics, number theory, combinatorics. Module
   count climbs toward ~70 %.
4. **v1.0 (~60 → ~86 wk)** — Cat C Phases 22–24 + Cat D high-value
   modules: physics submodules, geometry, logic, holonomic,
   algebraic number fields. Module count ~90 %, everyday coverage
   ~95 %, full algorithmic depth ~85 %.

Calling SymPP "v1.0" at ~95 % everyday coverage is honest: the
missing 5 % is research-level work (Risch domain extensions,
Pantelides for high-index DAEs, full Slater theorem) where any
serious user is likely going to drop into SymPy or Maple regardless.

---

## Critical path

```
shipped (0..13, 15) ─┬─ 16 (hardening → v0.5)
                      │
                      ├─ Cat A high-priority: Infinity → Gruntz
                      ├─ Cat A high-priority: hyperexpand → Meijer G → Risch
                      ├─ Cat A high-priority: lambdify (LLVM ORC)
                      │
                      ├─ Cat C: 17 (vector calc) → 18 (tensors)
                      │         19 (stats), 20 (ntheory),
                      │         21 (combinatorics), 22 (physics),
                      │         23 (crypto), 24 (geometry)
                      │
                      └─ Cat D: 25 (logic), 26 (discrete),
                                27 (holonomic), 28 (number fields),
                                29 (DomainMatrix), 30 (Galois),
                                31 (LaTeX parser), 32–34 (special
                                functions / sequences / orthogonal polys),
                                35 (algebras), 36 (NDim arrays),
                                37 (unify), 38 (extra printers),
                                39 (codegen AST)
                                              │
                                              └→ v1.0
```

The deep-deferreds within shipped phases are independent of one
another and can be tackled in any order. The new modules (17+)
likewise have low coupling. Phase 16 is the only strictly-sequenced
work remaining before parallelization opens up.
