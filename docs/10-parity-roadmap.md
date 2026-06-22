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
- **Discrete transforms** — fft/ifft, ntt/intt, convolution, Möbius.
- **General Jordan form** — chains of any length (reconstruction-verified).
- **DomainMatrix** — fraction-free ℤ/ℚ matrices (Bareiss det, rank, rref).
- **Berlekamp–Zassenhaus** — `factor_zassenhaus` univariate ℤ factoring.
- **Physics quantum/atomic** — arbitrary-spin operators, Wigner 3-j/6-j/9-j,
  Racah W, Gaunt + Clebsch–Gordan, Dirac γ-matrices, hydrogen & QHO
  energies/wavefunctions, Hamiltonian, Gaussian-beam optics.
- **ECC** — elliptic curves over 𝔽ₚ (group law, scalar mul, ECDH, ECDSA).
- **Combinatorics & group theory** — permutations, permutation groups, standard
  groups, integer partitions.
- **Number-theory extensions** — CRT, discrete log, linear Diophantine.
- **Hankel transform** — `hankel_transform` (order-ν, table-driven; Bessel-enabled).
- **Meijer-G product convolution** — `meijerg_convolution`, ∫₀^∞ G₁·G₂ via Mellin–Parseval, wired into `integrate`.
- **Bessel functions** — besselj/bessely/besseli/besselk (special values, ±1/2
  elementary closed forms, derivative recurrences, parser registration).
- **Meijer-G engine** — generic-case Slater reduction (Phase 1), Mellin-transform
  definite integration (Phase 3), function→Meijer-G recognition (Phase 4), and
  the **confluent G^{2,0}_{0,2} → modified Bessel K** closed form (Phase 2);
  OpenSpec `add-meijerg-slater-engine`. Remaining: the two-G `∫G₁·G₂` convolution.
- **Full SVD** — `Matrix::svd()` (U·Σ·Vᵀ), reconstruction-verified.
- **LaTeX parser** — `parse_latex`, round-trips with the LaTeX printer.
- **Physics core** — quantum (commutators/Pauli/ladder), ABCD optics,
  conjugate momentum.

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
| 🟡 Meijer-G integration (general method) | `integrals/meijerint.py` | 3 wk | High | Pairs with full hyperexpand. **Mellin master formula shipped** (OpenSpec `add-meijerg-slater-engine`, Phase 3+4): `meijerg_mellin_transform`, `meijerg_integrate_0_inf`, function→Meijer-G recognition (`to_meijerg`), so `∫₀^∞ xᵃe⁻ˣ = Γ(a+1)` routes through Meijer-G. Classic reducible integrals already evaluated via the heuristic engine. Remaining: grow the recognition table (trig/Bessel/erf), product `∫G₁·G₂` convolution, and general `integrate(…,0,∞)` dispatch |
| Full Risch transcendental integration | `integrals/risch.py` | 4 wk | Medium | Closes the remaining unsolved elementary integrals |
| 🟡 Full Slater/Meijer-G `hyperexpand` | `simplify/hyperexpand.py` | 2 wk | Medium | General hypergeometric closed forms. **Phase 1 shipped** (OpenSpec `add-meijerg-slater-engine`): generic-case Slater reduction `meijerg → Σ z^{b_k}·pFq`, wired into `hyperexpand` (G^{1,1}_{1,1}→1/(z+1), G^{2,0}_{0,2}→√π(cosh−sinh)(2√z), …). Remaining: confluent/log case, Mellin–Barnes definite integration, function→Meijer-G recognition |
| Multivariate `Poly` + Wang factorization | `polys/` | 3 wk | Medium | Multivariate factoring |
| ✅ Berlekamp–Zassenhaus | `polys/factortools.py` | — | shipped `factor_zassenhaus` — Berlekamp mod a small prime, **multifactor Hensel lifting** to a prime power above the Landau–Mignotte bound, then recombination |
| ✅ Symbolic SVD | `matrices/` | — | shipped — `singular_values()` and full `svd()` (U·Σ·Vᵀ) |
| ✅ General Jordan form (chains > 2) | `matrices/eigen.py` | — | shipped (filtration algorithm; reconstruction-verified) |
| ✅ Full 2D pretty-print layout | `printing/pretty` | — | shipped (block-layout `pretty()`) |
| Last Gruntz mrv-set rewrite | `series/gruntz.py` | 1 wk | Medium | Most stages shipped; the previously-flagged `0·∞` divergent-exp case (eˣ(e^{1/x}−1)→∞) now evaluates, as do x^{1/x}→1, (1+1/x)ˣ→e. **Practical MRV-set values are covered** — the dominant-summand split resolves mixed polynomial/exp/log-rate sums exactly (log(x²+eˣ)/x→1, log(x²+x)/log x→2, nested-log sums…; see LIMIT-MRV-VALUES-1 / LIMIT-LOGSUMDOM-1). **Termination is guaranteed** — the search is memoized and bounded by a work budget, so deeply nested exp-of-exp towers return an honest nan instead of hanging. Only the *exact value* of those pathological towers needs the full mrv-set comparability algorithm |
| F4/F5 Gröbner, sparse matrix, full polynomial domain tower, full Lie classifier, Pantelides DAE, full SAT-`ask` | various | ~15 wk | Low | Performance / edge-case depth |
| ✅ Non-commutative *algebra* (`Mul` ordering) | `core/mul` | — | shipped: `mul` preserves non-commutative operator order (A·B ≠ B·A), folds commutative factors into a leading coefficient, merges adjacent equal bases (A·A → A²); commutative fast path untouched |

## Category B — release engineering (Phase 16 → v1.0)

| Item | Effort | Priority | Notes |
|---|---|---|---|
| Benchmarks, doxygen API docs, vcpkg/Conan packaging, v1.0 tag | 3–5 wk | High | Mostly process/CI, not algorithmic; adoption-readiness. **Shipped:** `find_package`/install/export, a dependency-free benchmark harness (`SYMPP_BUILD_BENCHMARKS` → `sympp_bench`), and a `vcpkg.json` manifest. **Remaining:** Doxygen API-doc generation, Conan recipe, and the v1.0 tag |

## Category C — new modules (no MATLAB analogue)

| Module | SymPy ref | Effort | Priority | Notes |
|---|---|---|---|---|
| 🟡 Number theory (CRT, discrete_log, linear Diophantine, partitions) | `ntheory/` | — | shipped: factorint/divisors/igcdex/jacobi/continued_fraction/n_order/primitive_root/sqrt_mod, **CRT** (`crt`), **discrete log** (`discrete_log`), **linear Diophantine** (`diop_linear`), and integer partitions (in `combinatorics`). Remaining: general (non-linear) Diophantine — Pell, sum-of-squares, ternary quadratics |
| ✅ Statistics & probability (core distributions) | `stats/` | — | shipped (`stats/stats.hpp`); extend w/ more distributions, sampling |
| ✅ Geometry (Point/Line/Polygon) | `geometry/` | — | shipped (`geometry/geometry.hpp`) |
| ✅ Vector calculus & differential geometry | `vector/`, `diffgeom/` | — | shipped (`vector/vector_calculus.hpp`; grad/div/curl/laplacian + Christoffel/Ricci) |
| ✅ Tensor algebra (dense) | `tensor/` | — | shipped (`tensor/tensor.hpp`; product/contraction/raise/lower) |
| 🟡 Combinatorics & group theory | `combinatorics/` | — | shipped (`combinatorics/combinatorics.hpp`): permutations (compose/inverse/sign/order/cyclic form), permutation groups (closure/order/membership/abelian), S/C/D/A standard groups, integer partitions. Remaining: Sylow/Schreier–Sims fast group algorithms, Polýa enumeration |
| ✅ Cryptography (RSA/DH/ElGamal/ECC) | `crypto/` | — | shipped (`crypto/crypto.hpp`): RSA, Diffie–Hellman, ElGamal, and **elliptic curves over 𝔽ₚ** (group law, scalar mul, ECDH, ECDSA) |
| 🟡 Physics (mechanics, quantum, optics) | `physics/*` | 2 wk each | Medium | shipped (`physics/physics.hpp`): commutators/Pauli/ladder, **arbitrary-spin operators (Jx/Jy/Jz/J±/J²)**, **Wigner 3-j/6-j/9-j, Racah W, Gaunt + Clebsch–Gordan**, **Dirac γ-matrices**, **hydrogen E/R_nl**, **QHO E/ψ_n**, ABCD + Gaussian-beam optics, conjugate momentum + Hamiltonian. Remaining: second quantization, full continuum/relativistic mechanics, unit-bearing quantum states |

## Category D — modules outside the original 0–24 plan

| Module | SymPy ref | Effort | Priority | Notes |
|---|---|---|---|---|
| ✅ Special integral functions (`Ei,Si,Ci,Shi,Chi,fresnels,fresnelc,expint`) | `special.error_functions` | — | shipped (derivatives + special values) |
| Orthogonal polynomials (Legendre, Chebyshev, …) | `special.polynomials` | — | — | ✅ shipped (`functions/orthopolys.hpp`) |
| `rewrite(target)` cross-cutting API (exp↔trig…) | core | — | — | ✅ shipped (`core/rewrite.hpp`); extend with more targets |
| ✅ Logic & boolean algebra (`satisfiable`, `simplify_logic`) | `logic` | — | shipped (`logic/logic.hpp`) |
| ✅ LaTeX parser (round-trip) | `parsing.latex` | — | shipped (`parsing/latex_parser.hpp`) |
| ✅ Discrete (FFT/NTT/convolution/Möbius) | `discrete` | — | shipped (`discrete/discrete.hpp`) |
| ✅ DomainMatrix (fraction-free ℤ/ℚ matrices) | `polys.matrices` | — | shipped (`polys/domain_matrix.hpp`; Bareiss det, rank, rref) |
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
