# 03 — Feature Mapping

MATLAB Symbolic Math Toolbox features → SymPy reference source. Use this as the porting checklist.

Paths are relative to `sympy/sympy/` in the repository root. ★ marks the algorithmic hot-spot to study; sibling files are usually wiring/dispatch.

---

## Chapter 1–2 — Core, numbers, assumptions, units

| MATLAB feature | SymPy source | Phase |
|---|---|---|
| `sym`, `syms`, expression tree | ★`core/basic.py`, `core/expr.py`, `core/symbol.py`, `core/sympify.py` | 1 |
| `+ - * / ^` evaluation | ★`core/add.py`, `core/mul.py`, `core/power.py` | 1 |
| Number tower (exact) | ★`core/numbers.py` | 1 |
| `vpa` (variable-precision arithmetic) | `core/evalf.py` (uses `mpmath`) | 1 |
| Singletons (Pi, E, oo, NaN, …) | `core/singleton.py`, `core/numbers.py` | 1 |
| `subs`, `xreplace` | `core/basic.py` (`_subs`, `xreplace`) | 1 |
| Assumptions (`real`, `positive`, `integer`, …) — old system | ★`core/assumptions.py`, `assumptions/assume.py` | 2 |
| Assumptions — SAT-based `ask` | `assumptions/sathandlers.py`, `satask.py`, `cnf.py`, `lra_satask.py` | 2 (deferred sub-phase) |
| `refine` (apply assumptions to expr) | `assumptions/refine.py` | 2 |
| Units, SI/CGS/US conversion | `physics/units/quantities.py`, `dimensions.py`, `unitsystem.py`, `definitions/` | 12 |

## Chapter 3 — Mathematics: expression rearrangement

| MATLAB | SymPy source | Phase |
|---|---|---|
| `expand` | `core/function.py::expand`, `core/mul.py::_eval_expand_mul`, `core/power.py::_eval_expand_power_*` | 1 |
| `factor` (algebraic) | ★`polys/factortools.py`, `polys/polytools.py::factor` | 4 |
| `collect` | `simplify/radsimp.py::collect` | 5 |
| `combine` (logs/powers/sin·cos) | `simplify/powsimp.py`, `simplify/simplify.py::logcombine`, `simplify/trigsimp.py` | 5 |
| `rewrite(target)` | `core/basic.py::rewrite` + per-function `_eval_rewrite_as_*` methods on each `Function` subclass | 3 |
| `partfrac` | ★`polys/partfrac.py` | 4 |
| `numden`, `simplifyFraction` | `core/expr.py::as_numer_denom`, `simplify/simplify.py::nsimplify` | 1, 5 |
| `horner` | `polys/polyfuncs.py::horner` | 4 |
| `simplify` | ★`simplify/simplify.py`, plus all of `simplify/` | 5 |
| `trigsimp`, `expand_trig` | `simplify/trigsimp.py`, `simplify/fu.py` (Fu et al. table-driven) | 5 |
| `radsimp`, `rationalize` | `simplify/radsimp.py`, `simplify/sqrtdenest.py` | 5 |
| `powsimp`, `expand_power_*` | `simplify/powsimp.py` | 5 |
| `combsimp`, `gammasimp` | `simplify/combsimp.py`, `simplify/gammasimp.py` | 5 |
| `hyperexpand` | `simplify/hyperexpand.py` | 5 |

## Chapter 3 — Calculus

| MATLAB | SymPy source | Phase |
|---|---|---|
| `diff` | `core/function.py::Derivative`, `core/expr.py::diff` | 6 |
| `int` (indefinite, orchestrator) | ★`integrals/integrals.py` | 7 |
| `int` (rational fns) | `integrals/rationaltools.py` (Lazard-Rioboo-Trager) | 7 |
| `int` (rule-based / table) | ★`integrals/manualintegrate.py` | 7 (port first) |
| `int` (Risch heuristic) | `integrals/heurisch.py` | 7 |
| `int` (full Risch transcendental) | `integrals/risch.py`, `prde.py`, `rde.py` | 7 (deferred) |
| `int` (special-fn integrals) | `integrals/meijerint.py` | 7 |
| `int` (trig) | `integrals/trigonometry.py` | 7 |
| Definite `int` | `integrals/integrals.py::Integral.doit` (Newton-Leibniz + Meijer G) | 7 |
| `taylor` | `series/series.py`, `series/formal.py`, `core/expr.py::series` | 6 |
| `limit` (one-sided, two-sided) | ★`series/limits.py` + `series/gruntz.py` (Gruntz algorithm) | 6 |
| `symsum` | `concrete/summations.py`, `concrete/gosper.py`, `concrete/expr_with_intlimits.py` | 6 |
| `pade` | `series/approximants.py` | 6 |
| `asymptotes`, critical/inflection | compose via `limit` + `solve(diff(...))`; `calculus/util.py` helpers | 6 |
| Functional derivatives (Euler-Lagrange) | `calculus/euler.py` | 6 |
| `feval`, multivar minima/maxima | `calculus/util.py::minimum/maximum/stationary_points` | 6 |
| `vpaintegral` (numeric quadrature) | `integrals/quadrature.py`, mpmath via `evalf` | 7 |

## Chapter 3 — Transforms

| MATLAB | SymPy source | Phase |
|---|---|---|
| `fourier`, `ifourier` | `integrals/transforms.py` (Mellin-based) | 8 |
| `laplace`, `ilaplace` | ★`integrals/laplace.py` (newer rewrite) | 8 |
| `ztrans`, `iztrans` | `integrals/transforms.py` (Z-transform helpers) | 8 |
| `mellin`, `imellin` | `integrals/transforms.py` | 8 |
| `hankel`, `sine_transform`, `cosine_transform` | `integrals/transforms.py` | 8 (optional) |

## Chapter 3 — Equation & ODE/PDE solving

| MATLAB | SymPy source | Phase |
|---|---|---|
| `solve` (algebraic) | ★`solvers/solveset.py::solveset` (preferred), `solvers/solvers.py::solve` (legacy heuristics) | 10 |
| `vpasolve` (numeric) | `solvers/solvers.py::nsolve` (mpmath Newton) | 10 |
| `linsolve` | `solvers/solveset.py::linsolve`, `polys/solvers.py` | 10 |
| Polynomial systems | `solvers/polysys.py` (Gröbner-based, depends on Phase 4 Gröbner) | 10 |
| Inequalities (`isolve`) | `solvers/inequalities.py`, `solvers/simplex.py` | 10 |
| Diophantine | `solvers/diophantine/` | 10 (optional) |
| `dsolve` (ODE single) | ★`solvers/ode/single.py` (classifier), `solvers/ode/ode.py` (dispatcher) | 11 |
| `dsolve` (ODE systems) | `solvers/ode/systems.py`, `nonhomogeneous.py` | 11 |
| Riccati, hypergeometric, Lie | `solvers/ode/riccati.py`, `hypergeometric.py`, `lie_group.py` | 11 |
| `pdsolve` | `solvers/pde.py` | 11 |
| Recurrences (Z-transform domain) | `solvers/recurr.py` | 11 |

## Chapter 3 — Linear algebra

| MATLAB | SymPy source | Phase |
|---|---|---|
| Matrix construction (dense, sparse, immutable) | `matrices/dense.py`, `sparse.py`, `immutable.py`, `matrices.py`, `matrixbase.py` | 9 |
| `det`, `inv`, `rref`, `null`, `rank` | `matrices/determinant.py`, `inverse.py`, `reductions.py`, `subspaces.py` | 9 |
| `eig`, eigenvectors | ★`matrices/eigen.py` (uses Phase 4 polys for charpoly roots) | 9 |
| `jordan` | `matrices/eigen.py::jordan_form` | 9 |
| SVD | `matrices/decompositions.py::singular_value_decomposition` | 9 |
| LU, QR, Cholesky | `matrices/decompositions.py` | 9 |
| Norms (Frobenius, p-norms) | `matrices/matrixbase.py::norm` | 9 |
| Rotation matrices | `matrices/dense.py::rot_axis1/2/3`, `rot_ccw_axis*` | 9 |
| Matrix expressions / `MatrixSymbol` | `matrices/expressions/` | 9 |
| Hadamard, Kronecker products | `matrices/expressions/hadamard.py`, `kronecker.py` | 9 |
| `jacobian`, `hessian` | `matrices/matrixbase.py::jacobian`, `matrices/dense.py::hessian` | 9, 13 |

## Chapter 3 — Polynomials & number theory

| MATLAB | SymPy source | Phase |
|---|---|---|
| `Poly` constructor & operations | ★`polys/polytools.py` (user API), `densebasic.py`, `densearith.py`, `densetools.py` | 4 |
| GCD (single & multivariate) | `polys/euclidtools.py`, `heuristicgcd.py`, `modulargcd.py` | 4 |
| Factorization over ℤ/ℚ | ★`polys/factortools.py` (Berlekamp-Zassenhaus, Wang) | 4 |
| Roots (closed-form & RootOf) | `polys/polyroots.py`, `rootoftools.py`, `rootisolation.py` | 4 |
| Domains (ℤ, ℚ, ℝ, ℂ, ℤ_p, algebraic) | `polys/domains/` | 4 |
| Rings, ideals, Gröbner | `polys/rings.py`, `groebnertools.py`, `fglmtools.py` | 4 (Gröbner deferred) |
| Resultants, subresultants | `polys/multivariate_resultants.py`, `subresultants_qq_zz.py` | 4 |
| Square-free factorization | `polys/sqfreetools.py` | 4 |
| `factor` (integers), `isprime`, primes | `ntheory/factor_.py`, `generate.py`, `primetest.py`, `residue_ntheory.py` | 4 |

## Chapter 4 — Plotting

| MATLAB | SymPy source | SymPP approach |
|---|---|---|
| `fplot`, `fplot3` | `plotting/plot.py` | Bridge: emit data via `lambdify` → consumer plots |
| `fsurf`, `fmesh`, `fcontour`, `fimplicit` | `plotting/plot.py`, `plot_implicit.py` | Same |
| `animatedline`, animations | (matplotlib in SymPy) | Out of scope; consumer responsibility |

We ship `sympp::plotting::SampledFunction` — a sampled `Expr` over a domain — and provide adaptors for gnuplot and matplotlib-cpp.

## Chapter 5 — Code generation

| MATLAB | SymPy source | Phase |
|---|---|---|
| `ccode`, `fortran` | ★`printing/c.py`, `printing/fortran.py`, `printing/codeprinter.py`, `printing/precedence.py` | 13 |
| `latex` | `printing/latex.py` | 13 |
| `pretty`, ASCII pretty-print | `printing/pretty/` | 13 |
| `matlabFunction` (Octave-syntax) | `printing/octave.py` | 13 |
| C++ printer | `printing/cxx.py` | 13 |
| `lambdify` | `utilities/lambdify.py` | 13 (LLVM ORC backend) |
| CSE (common subexpression elimination) | ★`simplify/cse_main.py` | 13 |
| Higher-level codegen (functions, modules, files) | `codegen/algorithms.py`, `ast.py`, `cnodes.py`, `fnodes.py`, `cutils.py`, `futils.py` | 13 |
| Numeric autowrap (compile + load) | `utilities/autowrap.py` | 13 |
| `jacobian` for optimization | `matrices/matrixbase.py::jacobian` | 9 |

## Chapter 6 — Function reference (built-ins to port)

### Elementary
`functions/elementary/`:
- `trigonometric.py` — sin, cos, tan, cot, sec, csc, asin, acos, atan, atan2, acot, asec, acsc, sinc
- `exponential.py` — exp, log (any base), LambertW, exp_polar
- `hyperbolic.py` — sinh, cosh, tanh, coth, sech, csch, asinh, acosh, atanh, …
- `complexes.py` — re, im, abs, arg, conjugate, sign, polar_lift, periodic_argument
- `integers.py` — floor, ceiling, frac, round
- `miscellaneous.py` — sqrt, root, Min, Max, real_root, cbrt
- `piecewise.py` — Piecewise

### Special
`functions/special/`:
- `gamma_functions.py` — gamma, loggamma, digamma, polygamma, lowergamma, uppergamma
- `beta_functions.py` — beta, betainc, betainc_regularized
- `error_functions.py` — erf, erfc, erfi, erf2, erfinv, fresnels, fresnelc, expint, Ei, Si, Ci, Shi, Chi, li, Li
- `bessel.py` — besselj, bessely, besseli, besselk, hankel1, hankel2, jn, yn, hn1, hn2, airyai, airybi, airyaiprime, airybiprime
- `zeta_functions.py` — zeta, dirichlet_eta, polylog, lerchphi, stieltjes
- `hyper.py` — hyper (generalized hypergeometric), meijerg
- `polynomials.py` — Legendre, Chebyshev T/U, Hermite, Laguerre, Jacobi, Gegenbauer, assoc_legendre, assoc_laguerre
- `mathieu_functions.py` — Mathieu C, S, Cprime, Sprime
- `elliptic_integrals.py` — elliptic_k, elliptic_e, elliptic_pi, elliptic_f
- `delta_functions.py` — DiracDelta, Heaviside
- `singularity_functions.py` — SingularityFunction
- `bsplines.py` — bspline_basis, bspline_basis_set
- `spherical_harmonics.py` — Ynm, Znm
- `tensor_functions.py` — KroneckerDelta, LeviCivita

### Combinatorial
`functions/combinatorial/`:
- `factorials.py` — factorial, factorial2, RisingFactorial, FallingFactorial, binomial, subfactorial
- `numbers.py` — Fibonacci, Lucas, Bell, Bernoulli, Euler, Catalan, Genocchi, Stirling, harmonic, motzkin, partition

These are mostly **mechanical to port**: each function is a class with `_eval_evalf`, `_eval_rewrite_as_*` methods, derivative rules, and series expansion rules. Pattern is highly repetitive.

---

## Mapping to SymPP phases

See [04-roadmap.md](04-roadmap.md) for the phase-by-phase plan. Each phase number above corresponds to that document.
