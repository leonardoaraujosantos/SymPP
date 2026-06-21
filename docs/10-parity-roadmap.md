# Remaining-Work Roadmap (SymPy parity)

A focused, execution-ordered view of what is left to reach full SymPy parity,
refreshed after the assumptions-parity push and the `lambdify` + number-theory
deliverables. Complements the per-phase history in
[04-roadmap.md](04-roadmap.md); this file is the *forward* plan.

Effort figures are focused developer-weeks at the project's scope discipline
(clean-room, oracle-validated, no regressions). Priority blends real-world
demand with cost.

## Recently shipped (this line of work)

- **Assumptions** — completed the predicate ontology (28 predicates incl.
  extended-signed, hermitian/antihermitian, commutative), boolean/SAT-style
  `ask`, and an oracle parity guard (`tests/core/parity_probe.cpp`).
- **`lambdify`** — portable closure-based numeric compiler (interpreter
  backend; the LLVM-JIT variant remains below).
- **Number theory** — `factorint`, `divisors`, `igcdex`, `jacobi_symbol`
  (GMP-backed, oracle-validated).

## How to read effort vs. session-size

Most items below are **multi-week research-grade algorithms** — they are *not*
completable in a single working session and must be landed as their own
scoped changes (ideally one OpenSpec proposal each). The handful of
**session-sized** items are flagged ⚡; tackle those first for quick parity
gains.

---

## Category A — depth inside shipped subsystems

| Item | SymPy ref | Effort | Priority | Notes |
|---|---|---|---|---|
| `lambdify` LLVM-JIT backend | `utilities.lambdify` | 3 wk | High | Interpreter backend shipped; JIT is the perf tier |
| Meijer-G integration | `integrals/meijerint.py` | 3 wk | High | Pairs with full hyperexpand; unlocks special-fn integrals |
| Full Risch transcendental integration | `integrals/risch.py` | 4 wk | Medium | Closes the remaining unsolved elementary integrals |
| Full Slater/Meijer-G `hyperexpand` | `simplify/hyperexpand.py` | 2 wk | Medium | General hypergeometric closed forms |
| Multivariate `Poly` + Wang factorization | `polys/` | 3 wk | Medium | Multivariate factoring |
| Berlekamp–Zassenhaus + Hensel lifting | `polys/factortools.py` | 2 wk | Medium | Robust univariate factoring over ℤ |
| Symbolic SVD | `matrices/` | 2 wk | Medium | Needs eigendecomposition of AᵀA |
| General Jordan form (chains > 2) | `matrices/eigen.py` | 1 wk | Medium | Filtration algorithm; eigenvalue-limited |
| Full 2D pretty-print layout | `printing/pretty` | 1 wk | Medium | ⚡ partial `pretty()` exists; extend incrementally |
| Last Gruntz mrv-set rewrite | `series/gruntz.py` | 1 wk | High | Most stages shipped; remaining `0·∞` divergent-exp subclass |
| F4/F5 Gröbner, sparse matrix, full polynomial domain tower, full Lie classifier, Pantelides DAE, full SAT-`ask` | various | ~15 wk | Low | Performance / edge-case depth |
| Non-commutative *algebra* (`Mul` ordering) | `core/mul` | 2–3 wk | Low–Med | The `commutative` *predicate* already ships |

## Category B — release engineering (Phase 16 → v1.0)

| Item | Effort | Priority | Notes |
|---|---|---|---|
| Benchmarks, doxygen API docs, vcpkg/Conan packaging, v1.0 tag | 4–6 wk | High | Mostly process/CI, not algorithmic; adoption-readiness |

## Category C — new modules (no MATLAB analogue)

| Module | SymPy ref | Effort | Priority | Notes |
|---|---|---|---|---|
| Number theory (continued, primitive_root, sqrt_mod, Diophantine) | `ntheory/` | 1.5 wk | High | ⚡ core primitives shipped; extend incrementally |
| Statistics & probability | `stats/` | 3 wk | High | Distributions, pdf/cdf/E/Var |
| Geometry | `geometry/` | 1 wk | Medium | ⚡ self-contained |
| Vector calculus & differential geometry | `vector/`, `diffgeom/` | 3 wk | Medium | grad/div/curl, curvature |
| Tensor algebra | `tensor/` | 2 wk | Medium | Index notation, contraction |
| Combinatorics & group theory | `combinatorics/` | 3 wk | Medium | Permutations, finite groups |
| Cryptography | `crypto/` | 2 wk | Low–Med | RSA, ECC, finite fields |
| Physics (mechanics, quantum, optics) | `physics/*` | 4 wk each | Medium | Large, niche submodules |

## Category D — modules outside the original 0–24 plan

| Module | SymPy ref | Effort | Priority | Notes |
|---|---|---|---|---|
| Special integral functions (`Ei,Si,Ci,Shi,Chi,fresnel,expint`) | `special.error_functions` | 1 wk | Medium | ⚡ enums exist; add eval/series |
| Orthogonal polynomials (Legendre, Chebyshev, …) | `special.polynomials` | 1.5 wk | Medium | ⚡ |
| `rewrite(target)` cross-cutting API (exp↔trig…) | core | ~1 wk | Medium | ⚡ notable missing convenience |
| Logic & boolean algebra (`satisfiable`, `simplify_logic`) | `logic` | 2 wk | Medium | Also feeds fuller SAT-`ask` |
| LaTeX parser (round-trip) | `parsing.latex` | 2 wk | Medium | Input ergonomics |
| Discrete (FFT/NTT/convolution/Möbius) | `discrete` | 2 wk | Low–Med | |
| DomainMatrix (fast poly-domain matrices) | `polys.matrices` | 3 wk | Medium | Order-of-magnitude speedups |
| Holonomic fns, algebraic number fields, Galois tools, quaternions, NDim arrays, unification, extra printers (MathML/Rust/Julia), codegen AST/autowrap | various | ~13 wk | Low | |

---

## Recommended execution order

Optimizing impact-per-week, starting with the ⚡ session-sized wins:

1. **Cheap, high-value, session-sized** (⚡): extend number theory
   (continued fractions, `sqrt_mod`, `primitive_root`), special integral
   functions eval, orthogonal polynomials, `rewrite(target)`, geometry. Each
   lands as its own oracle-validated change.
2. **Flagship algorithm tier**: Meijer-G integration + full `hyperexpand`
   (do together), then the last Gruntz rewrite, then `lambdify` LLVM-JIT.
3. **Polynomial depth**: Berlekamp–Zassenhaus, then multivariate Wang.
4. **Linear algebra depth**: general Jordan form, then symbolic SVD.
5. **New modules**: statistics → vector calculus → tensors → combinatorics.
6. **Release**: Category B (Phase 16 → v1.0) once the above stabilizes.

**Total remaining to full module-for-module parity: ~80 developer-weeks.**
Everyday-CAS-workflow parity is already ≈85% (see
[04-roadmap.md](04-roadmap.md#how-far-are-we-from-sympy)); the list above is
what separates that from *complete* SymPy coverage.
