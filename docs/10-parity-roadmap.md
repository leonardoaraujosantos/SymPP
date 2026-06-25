# Remaining-Work Roadmap (SymPy parity)

A focused, execution-ordered view of what is left to reach full SymPy parity,
refreshed after the assumptions-parity push and the `lambdify` + number-theory
deliverables. Complements the per-phase history in
[04-roadmap.md](04-roadmap.md); this file is the *forward* plan.

Effort figures are focused developer-weeks at the project's scope discipline
(clean-room, oracle-validated, no regressions). Priority blends real-world
demand with cost.

## Recently shipped (this line of work)

- **Extra code printers** — `rust_code`, `julia_code`, `glsl_code`, `dot`,
  `srepr`, and Presentation `mathml`
  join the existing C/C++/Fortran/LaTeX/Octave set (`.powi/.powf/.sqrt`
  method-call form + `f64::consts` for Rust, `^`/`MathConstants.e` for Julia,
  `<msup>/<mfrac>/<msqrt>` MathML), each with exact-string regression tests.
- **solveset domains** — added a `Naturals`/`Naturals0` set domain wired into the
  named-domain lattice (ℕ ⊂ ℤ ⊂ ℝ ⊂ ℂ) and used to restrict polynomial roots
  (`solveset(…, naturals())`); `solveset` with no domain now defaults to the
  **complex domain** (SymPy parity: `solveset(x²+1) = {−I, I}`), while an explicit
  `reals()` domain still drops complex roots (SOLVESET-NAT-1, SOLVESET-COMPLEX-1).
- **Definite erf/erfc integrals** — `∫₀^∞ erfc(ax)² = (2−√2)/(√π·a)` and
  `∫₀^∞ x·e^{−cx²}·erf(bx) = b/(2c√(b²+c))`, recognized ahead of the Meijer-G
  fallback and verified numerically (INT-ERF-1).
- **Cubic+ bivariate factorization** — `factor_multivariate` now deflates a
  closed-form monomial root `±(p/q)·yᵏ` and recurses, so the quadratic
  discriminant path handles cubic/quartic inputs (x³−y³, x⁴−y⁴, …); every result
  is expand-verified, no-root inputs are rejected (WANG-2).
- **Pólya enumeration** — cycle-index polynomial `Z(G)` of a permutation group
  plus `necklaces(n, k)` ((1/n)Σ_{d|n} φ(d)·k^{n/d}), cross-checked against
  `colorings_count` and the oracle.
- **Second-quantization Fock operators** — single-mode bosonic `FockState`
  (a/a†/N, `[a,a†]=1`) and fermionic `FermionState` with Pauli exclusion
  (a/a†/N, `{a,a†}=1`), alongside the existing Jordan–Wigner fermions.
- **Normal ordering** — `normal_order(word, Statistics)` rewrites a product of
  (multi-mode) ladder operators into the equivalent linear combination of
  normal-ordered words, moving creations to the left via `a_i a_j† = δ_ij +
  a_j† a_i` (boson) / `δ_ij − a_j† a_i` (fermion) and accumulating contractions
  (fermion words with a repeated factor vanish).
- **Sylow p-subgroups** — `sylow_order(G, p)` and `sylow_subgroup(G, p)` (BFS
  p-closure to the exact prime-power order), validated on S₃/S₄/A₄/D₄.
- **Simplicity & abelianization** — `normal_closure(G, gens)` (BFS conjugation +
  subgroup closure), `is_simple(G)` (every nontrivial element's normal closure
  is all of G; A₅ simple, A₄/S₃/S₄/D₄/C₆ not), and `abelian_invariants(G)` (the
  prime-power invariant factors of G/G': C₆→[2,3], D₄→[2,2], S₃/S₄→[2], A₄→[3],
  A₅→[] perfect), cross-checked against SymPy.
- **hyperexpand squared-argument forms** — ₀F₁(;3/2;−z²/4)=sin(z)/z,
  ₀F₁(;1/2;−z²/4)=cos(z), ₂F₁(½,½;3/2;z²)=asin(z)/z, giving radical-free results
  where the generic √z forms would leave `sqrt(z²)` (HYPER-ELEM-1).
- **Gruntz competing infinity-powers** — products of two diverging powers sharing
  an exponent kernel (`xˣ/(x+1)ˣ → e⁻¹`, `(x+1)ˣ/xˣ → e`) are pre-normalized into
  the power-form machinery instead of hanging the L'Hôpital search
  (LIMIT-MRV-TOWER-1).
- **Risch coefficiented-power-denominator log integrals** — `(c·g)ⁿ`
  denominators (`log(x)/(2x²)`, `log(x)²/(3x³)`, …) are distributed to `cⁿ·gⁿ`
  before the rational-times-log dispatch, closing integrals that bare `gⁿ`
  matchers missed (INT-RLOG-1).
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
- **Symbolic Product (∏)** — `product(expr, var, lo, hi)` closed forms: factorial
  `∏k=n!`, constant-power, geometric-exponent `∏2^k=2^{n(n+1)/2}`, and
  polynomial-in-k via Gamma ratios (`∏(k²−1)`, telescoping `∏(k+1)/(k+2)→2/(n+2)`);
  the concrete counterpart to `summation` (PROD-1/PROD-2).
- **Gosper's algorithm** — finite hypergeometric summation (`Σk·k!=(n+1)!−1`,
  hockey-stick `ΣC(k,3)=C(n+1,4)`, `Σk·2^k`) via term-ratio → Gosper–Petkovšek
  normal form → certificate solve, gated by the rational Gosper identity and a
  pole guard; non-summable terms (`k!`, `k²·k!`) stay unevaluated as in SymPy
  (SUM-GOSPER-1).
- **Zeilberger creative telescoping** — parametric binomial sums whose value is
  first-order P-recursive: `Σ k²·C(n,k)=n(n+1)2^{n−2}`, `Σ k³·C(n,k)`,
  `Σ C(2n,2k)=2^{2n−1}`, found by exact recurrence discovery + `product()` solve,
  verified on held-out points; higher-order sums (Franel `Σ C(n,k)³`) stay
  unevaluated as in SymPy (SUM-ZEILBERGER-1).
- **Beta-integral binomial sums** — `Σ_{k=0}^n C(n,k)/(k+m)` for a positive
  integer `m`, whose value is order-2 P-recursive (mixes `2ⁿ` and a constant) and
  so escapes the first-order Zeilberger solve: `Σ C(n,k)/(k+1)=(2^{n+1}−1)/(n+1)`,
  `Σ C(n,k)/(k+2)=(2^{n+2}−1)/(n+2)−(2^{n+1}−1)/(n+1)`, from
  `∫₀¹ x^{m−1}(1+x)ⁿ dx`, verified numerically (SUM-BINOM-RECIP-1).
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
  OpenSpec `add-meijerg-slater-engine`, plus the two-G `∫G₁·G₂` Mellin–Parseval convolution. The Meijer-G engine is now feature-complete for the common cases.
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
| 🟡 Meijer-G integration (general method) | `integrals/meijerint.py` | 3 wk | High | Pairs with full hyperexpand. **Mellin master formula shipped** (OpenSpec `add-meijerg-slater-engine`, Phase 3+4): `meijerg_mellin_transform`, `meijerg_integrate_0_inf`, function→Meijer-G recognition (`to_meijerg`), so `∫₀^∞ xᵃe⁻ˣ = Γ(a+1)` routes through Meijer-G. Recognition table (exp/sin/cos via η·xᶜ → Gaussian/Dirichlet/Fresnel), product `∫G₁·G₂` convolution, and general `integrate(…,0,∞)` dispatch all shipped. **Definite erf/erfc integrals over [0,∞) shipped** — `∫₀^∞ erfc(ax)²=(2−√2)/(√π·a)`, `∫₀^∞ x·e^{−cx²}·erf(bx)=b/(2c√(b²+c))`, recognized ahead of the Meijer-G fallback (INT-ERF-1). Remaining: further erf-core table growth |
| 🟡 Full Risch transcendental integration | `integrals/risch.py` | 3 wk | Medium | **Exponential Risch differential equation shipped** — ∫R(x)·e^{ax} solved via Q′+aQ=R (undetermined coefficients, Q=U/gcd(V,V′)), closing rational·exp cases like ∫x·eˣ/(x+1)²=eˣ/(x+1). **Logarithmic-extension integration complete** — rational functions of log(x) reduce via u=log x to full rational integration (atan(log x), 1/log x, nested log towers; non-elementary cases stay unevaluated, matching SymPy; RISCH-LOG-1); coefficiented power denominators `(c·g)ⁿ` (log(x)/(2x²), log(x)²/(3x³)) are distributed to `cⁿ·gⁿ` before dispatch so the rational-times-log matchers fire (INT-RLOG-1). **Generalized-exponential tower shipped** — powers a(x)^{b(x)} with the variable in the exponent are treated as true exponentials e^{b·log a}, closing self-derivative cases ∫(1+log x)·x^x=x^x, ∫(2x·log x+x)·x^{x²}=x^{x²} (verify-guarded; RISCH-POWEXP-1). Mixed exp/log integrands are now at SymPy parity — closed forms where they exist, unevaluated where genuinely non-elementary. Remaining: deeper multi-level tower recursion and formal non-elementarity proofs (note: SymPy also returns these as special functions or unevaluated rather than proving non-elementarity) |
| 🟡 Full Slater/Meijer-G `hyperexpand` | `simplify/hyperexpand.py` | 2 wk | Medium | General hypergeometric closed forms. **Phase 1 shipped** (OpenSpec `add-meijerg-slater-engine`): generic-case Slater reduction `meijerg → Σ z^{b_k}·pFq`, wired into `hyperexpand` (G^{1,1}_{1,1}→1/(z+1), G^{2,0}_{0,2}→√π(cosh−sinh)(2√z), …). confluent→Bessel-K, Mellin–Barnes definite integration, and function→Meijer-G recognition all shipped. **Squared-argument elementary forms shipped** — ₀F₁(;3/2;−z²/4)=sin(z)/z, ₀F₁(;1/2;−z²/4)=cos(z), ₂F₁(½,½;3/2;z²)=asin(z)/z give radical-free results where the generic √z forms would leave `sqrt(z²)` (HYPER-ELEM-1) |
| 🟡 Multivariate `Poly` + Wang factorization | `polys/` | 2 wk | Medium | **Bivariate Wang shipped, now cubic-and-higher** — `factor()` handles genuinely-multivariate (symbolic-coefficient) inputs via content extraction over the second variable plus a quadratic-in-var solve through a perfect-square discriminant (WANG-1), extended with a monomial-root deflation branch that peels a closed-form root `±(p/q)·yᵏ` and recurses, so cubic/quartic inputs reduce into the discriminant path (x³−y³→(x−y)(x²+xy+y²), x⁴−y⁴, x³+yx²−x−y→(x−1)(x+1)(x+y), …; WANG-2). Every result is expand-verified, no-root inputs (x³+x+y) rejected. Remaining: general multivariate Hensel lifting and >2 variables |
| ✅ Berlekamp–Zassenhaus | `polys/factortools.py` | — | shipped `factor_zassenhaus` — Berlekamp mod a small prime, **multifactor Hensel lifting** to a prime power above the Landau–Mignotte bound, then recombination |
| ✅ Symbolic SVD | `matrices/` | — | shipped — `singular_values()` and full `svd()` (U·Σ·Vᵀ) |
| ✅ General Jordan form (chains > 2) | `matrices/eigen.py` | — | shipped (filtration algorithm; reconstruction-verified) |
| ✅ Full 2D pretty-print layout | `printing/pretty` | — | shipped (block-layout `pretty()`) |
| 🟡 Last Gruntz mrv-set rewrite | `series/gruntz.py` | 1 wk | Medium | Most stages shipped; the previously-flagged `0·∞` divergent-exp case (eˣ(e^{1/x}−1)→∞) now evaluates, as do x^{1/x}→1, (1+1/x)ˣ→e. **Practical MRV-set values are covered** — the dominant-summand split resolves mixed polynomial/exp/log-rate sums exactly (log(x²+eˣ)/x→1, log(x²+x)/log x→2, nested-log sums…; see LIMIT-MRV-VALUES-1 / LIMIT-LOGSUMDOM-1). **Competing infinity-powers shipped** — products of two diverging powers sharing an exponent kernel (xˣ/(x+1)ˣ→e⁻¹, (x+1)ˣ/xˣ→e) are merged into the power-form path instead of hanging the L'Hôpital search (LIMIT-MRV-TOWER-1). **Termination is guaranteed** — the search is memoized and bounded by a work budget, so deeply nested exp-of-exp towers return an honest nan instead of hanging. Only the *exact value* of those pathological towers needs the full mrv-set comparability algorithm |
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
| 🟡 Combinatorics & group theory | `combinatorics/` | — | shipped (`combinatorics/combinatorics.hpp`): permutations (compose/inverse/sign/order/cyclic form), permutation groups (closure/order/membership/abelian), S/C/D/A standard groups, integer partitions, **Schreier–Sims BSGS** (order/membership/transitivity), **Sylow p-subgroups** (`sylow_order`/`sylow_subgroup`), **conjugacy/center/derived series** (`conjugacy_classes`/`group_center`/`derived_subgroup`/`is_solvable`/`is_nilpotent`), **simplicity & abelianization** (`normal_closure`/`is_simple`/`abelian_invariants`), and **Pólya enumeration** (cycle-index polynomial `Z(G)` + `necklaces(n,k)`). Remaining: full Schreier–Sims base optimization, general Burnside/Pólya weight inventories |
| ✅ Cryptography (RSA/DH/ElGamal/ECC) | `crypto/` | — | shipped (`crypto/crypto.hpp`): RSA, Diffie–Hellman, ElGamal, and **elliptic curves over 𝔽ₚ** (group law, scalar mul, ECDH, ECDSA) |
| 🟡 Physics (mechanics, quantum, optics) | `physics/*` | 2 wk each | Medium | shipped (`physics/physics.hpp`): commutators/Pauli/ladder, **arbitrary-spin operators (Jx/Jy/Jz/J±/J²)**, **Wigner 3-j/6-j/9-j, Racah W, Gaunt + Clebsch–Gordan**, **Dirac γ-matrices**, **hydrogen E/R_nl**, **QHO E/ψ_n**, ABCD + Gaussian-beam optics, conjugate momentum + Hamiltonian, **second quantization (Jordan–Wigner fermions + single-mode bosonic `FockState` and fermionic `FermionState` ladder/number operators with the canonical [a,a†]/{a,a†} relations)**. Remaining: full continuum/relativistic mechanics, unit-bearing quantum states |

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
| ✅ Extra printers (MathML / Rust / Julia / GLSL / dot / srepr) | `sympy.printing` | — | shipped: `rust_code`, `julia_code`, `glsl_code`, `dot`, `srepr`, Presentation `mathml` join the C/C++/Fortran/LaTeX/Octave set (Phase 13) |
| Holonomic fns, algebraic number fields, Galois tools, quaternions, NDim arrays, unification, codegen AST/autowrap | various | ~12 wk | Low | |

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
