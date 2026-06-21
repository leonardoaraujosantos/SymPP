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
- **Number theory** — `factorint`, `divisors`, `igcdex`, `jacobi_symbol`,
  `continued_fraction`, `n_order`, `primitive_root`, `sqrt_mod`
  (GMP-backed, oracle-validated).
- **Orthogonal polynomials** — `legendre`, `chebyshevt`, `chebyshevu`,
  `hermite`, `laguerre` (recurrence-based, oracle-validated to degree 8).
- **`rewrite(target)`** — trig/hyperbolic ↔ `exp` and tan/cot/sec/csc → sin/cos,
  value-preserving across the tree (numerically certified vs SymPy).
- **Logic & boolean algebra** — And/Or/Not (+xor/implies/equivalent),
  `satisfiable`, `to_cnf`/`to_dnf`, `simplify_logic` (Quine–McCluskey).
- **Full 2D pretty-print** — block-layout `pretty()` (stacked fractions,
  superscript powers, radicals) matching SymPy's ASCII output.
- **`lambdify` LLVM ORC-JIT backend** — `lambdify_jit` compiles to native code
  (optional, auto-on when LLVM is found).
- **Geometry** — Point/Line/Polygon: distance, area, slope, intersection, …
- **Statistics** — Normal/Uniform/Exponential + Bernoulli/Binomial/Poisson with
  mean/variance/pdf/cdf.
- **Symbolic singular values** — `Matrix::singular_values()` (full SVD
  factorization still pending).
- **Special integral functions** — Fresnel integrals + generalized `expint`
  (Ei/Si/Ci/Shi/Chi already shipped).
- **Vector calculus & differential geometry** — grad/div/curl/laplacian +
  Christoffel/Ricci.
- **Tensor algebra** — dense tensors: product, contraction, raise/lower.
- **Cryptography** — RSA, Diffie–Hellman, ElGamal.

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
| ✅ `lambdify` LLVM-JIT backend | `utilities.lambdify` | — | shipped (`core/lambdify_llvm.hpp`, optional/auto-on) |
| Meijer-G integration | `integrals/meijerint.py` | 3 wk | High | Pairs with full hyperexpand; unlocks special-fn integrals |
| Full Risch transcendental integration | `integrals/risch.py` | 4 wk | Medium | Closes the remaining unsolved elementary integrals |
| Full Slater/Meijer-G `hyperexpand` | `simplify/hyperexpand.py` | 2 wk | Medium | General hypergeometric closed forms |
| Multivariate `Poly` + Wang factorization | `polys/` | 3 wk | Medium | Multivariate factoring |
| Berlekamp–Zassenhaus + Hensel lifting | `polys/factortools.py` | 2 wk | Medium | Robust univariate factoring over ℤ |
| Symbolic SVD | `matrices/` | 1 wk | Medium | `singular_values()` ✅ shipped; full U·Σ·Vᴴ factorization remains |
| General Jordan form (chains > 2) | `matrices/eigen.py` | 1 wk | Medium | Filtration algorithm; eigenvalue-limited |
| ✅ Full 2D pretty-print layout | `printing/pretty` | — | shipped (block-layout `pretty()`) |
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
| Number theory (Diophantine, CRT, discrete_log, partitions) | `ntheory/` | 1 wk | Medium | core + factorint/primitive_root/sqrt_mod/continued_fraction shipped |
| ✅ Statistics & probability (core distributions) | `stats/` | — | shipped (`stats/stats.hpp`); extend w/ more distributions, sampling |
| ✅ Geometry (Point/Line/Polygon) | `geometry/` | — | shipped (`geometry/geometry.hpp`) |
| ✅ Vector calculus & differential geometry | `vector/`, `diffgeom/` | — | shipped (`vector/vector_calculus.hpp`; grad/div/curl/laplacian + Christoffel/Ricci) |
| ✅ Tensor algebra (dense) | `tensor/` | — | shipped (`tensor/tensor.hpp`; product/contraction/raise/lower) |
| Combinatorics & group theory | `combinatorics/` | 3 wk | Medium | Permutations, finite groups |
| ✅ Cryptography (RSA/DH/ElGamal) | `crypto/` | — | shipped (`crypto/crypto.hpp`); ECC remains |
| Physics (mechanics, quantum, optics) | `physics/*` | 4 wk each | Medium | Large, niche submodules |

## Category D — modules outside the original 0–24 plan

| Module | SymPy ref | Effort | Priority | Notes |
|---|---|---|---|---|
| ✅ Special integral functions (`Ei,Si,Ci,Shi,Chi,fresnels,fresnelc,expint`) | `special.error_functions` | — | shipped (derivatives + special values) |
| Orthogonal polynomials (Legendre, Chebyshev, …) | `special.polynomials` | — | — | ✅ shipped (`functions/orthopolys.hpp`) |
| `rewrite(target)` cross-cutting API (exp↔trig…) | core | — | — | ✅ shipped (`core/rewrite.hpp`); extend with more targets |
| ✅ Logic & boolean algebra (`satisfiable`, `simplify_logic`) | `logic` | — | shipped (`logic/logic.hpp`) |
| LaTeX parser (round-trip) | `parsing.latex` | 2 wk | Medium | Input ergonomics |
| Discrete (FFT/NTT/convolution/Möbius) | `discrete` | 2 wk | Low–Med | |
| DomainMatrix (fast poly-domain matrices) | `polys.matrices` | 3 wk | Medium | Order-of-magnitude speedups |
| Holonomic fns, algebraic number fields, Galois tools, quaternions, NDim arrays, unification, extra printers (MathML/Rust/Julia), codegen AST/autowrap | various | ~13 wk | Low | |

---

## Recommended execution order

Optimizing impact-per-week, starting with the ⚡ session-sized wins:

1. **Cheap, high-value, session-sized** (⚡): *number-theory extensions,
   orthogonal polynomials and `rewrite(target)` are done.* Remaining quick
   wins: special integral functions eval (`Ei/Si/Ci`), geometry core
   (Point/Line/Segment/Polygon), and more `rewrite` targets. Each lands as its
   own oracle-validated change.
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
