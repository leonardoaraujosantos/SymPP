# 04 — Roadmap

This file is the source of truth for what SymPP has shipped, what's
deferred (and *honestly why*), and the forward path to SymPy feature
parity.

## Status snapshot

```
820 tests / 1644 assertions  all passing
53 commits on origin/main
13 of 15 phases shipped (Phase 14 dropped — see below)
```

| Phase | Title | Status |
|---|---|---|
| 0  | Foundation & oracle harness            | ✅ shipped |
| 1  | Core expression tree                   | ✅ shipped |
| 2  | Assumptions                            | 🟡 minimal subset (SAT-based deferred) |
| 3  | Elementary & special functions         | ✅ shipped |
| 4  | Polynomials                            | ✅ shipped (multivariate factor + BZ deferred) |
| 5  | Simplification                         | ✅ shipped (hyperexpand deferred) |
| 6  | Calculus                               | ✅ shipped (Gruntz at infinity deferred) |
| 7  | Integration                            | ✅ shipped (full Risch + Meijer G deferred) |
| 8  | Transforms                             | ✅ shipped (Meijer-driven general case deferred) |
| 9  | Linear algebra                         | ✅ shipped (SVD, Jordan, sparse deferred) |
| 10 | Equation solvers                       | ✅ shipped (F4/F5 Gröbner + SAT deferred) |
| 11 | ODE / PDE                              | ✅ shipped (full Lie + Pantelides deferred) |
| 12 | Units                                  | ✅ shipped |
| 13 | Code generation                        | 🟡 printers + function emission (lambdify deferred) |
| 15 | Parser & MATLAB facade                 | ❌ not started |
| 16 | Hardening & v1.0                       | ❌ not started |

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

### Phase 15 — Parser / serialization · ❌

**Scope**: a recursive-descent string parser for math expressions:
`parse("x**2 + sin(x)")` → `Expr`. Useful for (1) consumers reading
expressions from non-C++ sources (config files, databases, network
protocols), (2) round-tripping with SymPy via `srepr` for
serialization. ~1.5 weeks.

**Out of scope**: the MATLAB facade (`sympp::matlab::solve`, etc.)
that was in the original Phase 15 plan. The MATLAB-toolbox-parity
framing has served its purpose; SymPP is now reframed as a C++
alternative to SymPy, and a thin renaming-wrapper namespace adds
maintenance burden without real value (MATLAB users porting
scripts will rewrite anyway, and `int`/`diff` collide with C++
conventions).

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
| 5  | hyperexpand + hypergeometric function infrastructure | 4 wk | High (unblocks Phase 7+8) |
| 6  | Full Gruntz limit algorithm                          | 2 wk | High |
| 7  | Full Risch transcendental integration                | 4 wk | Medium |
| 7  | Meijer G-function integration method                 | 3 wk | High (with hyperexpand) |
| 9  | Symbolic SVD                                         | 2 wk | Medium |
| 9  | Jordan canonical form                                | 2 wk | Medium |
| 9  | Sparse Matrix variant                                | 2 wk | Low |
| 10 | F4 / F5 Gröbner algorithms                           | 3 wk | Low |
| 10 | General transcendental solveset                      | 2 wk | Medium |
| 11 | Full Lie symmetry classifier                         | 3 wk | Low |
| 11 | Pantelides BLT + Tarjan SCC for high-index DAEs      | 2 wk | Low |
| 13 | lambdify (LLVM ORC integration)                      | 3 wk | High (huge user value) |
| 13 | Full 2D pretty-print layout                          | 1 wk | Medium |

**Total category A**: roughly 42 developer-weeks of focused work.

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

## Summary

```
Where we are:        13/15 phases shipped, 820 oracle-validated tests.
What's deferred:     ~42 weeks of Category A (deep-deferreds within shipped phases).
What's next:         ~6–8 weeks for Category B (Phases 15, 16).
What's beyond:       ~25 weeks for Category C (Phases 17–24, true SymPy parity).
```

**Total to full SymPy parity**: roughly 75 focused developer-weeks
on top of what's already in `main`.

The first compelling milestone is **v0.5 = Phases 15 + 16 landing**:
that gives us a parser and a tagged release. From there, Category
A and C can be parallelized — the deep-deferreds within shipped
phases are independent of the new modules.

---

## Critical path

```
shipped (0..13) ─┬─ 15 (parser / serialization)
                  ├─ 16 (hardening → v0.5)
                  │
                  ├─ Cat A: hyperexpand → Meijer G → Risch
                  ├─ Cat A: Infinity → Gruntz
                  ├─ Cat A: lambdify (LLVM ORC)
                  │
                  └─ Cat C: 17 (vector calc) → 18 (tensors)
                            19 (stats), 20 (ntheory),
                            21 (combinatorics), 22 (physics),
                            23 (crypto), 24 (geometry)
                                        │
                                        └→ v1.0
```

The deep-deferreds within shipped phases are independent of one
another and can be tackled in any order. The new modules (17+)
likewise have low coupling. Phases 15 and 16 are the only
strictly-sequenced work remaining.
