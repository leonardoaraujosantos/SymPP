# Tasks

Multi-session initiative. Each phase is an independent green, oracle-validated
increment on `main`. Only Phase 1 is implemented in the first session.

## Phase 1 — generic Slater reduction (this session)

- [x] 1.1 Extract `(an, ap_rest, bm, bq_rest, z)` from a MeijerG node
- [x] 1.2 Generic-case guard: bail (leave opaque) when two of `b₁…b_m` differ by
  an integer
- [x] 1.3 Build the Slater sum `Σ A_k z^{b_k} pF_{q−1}(…; (−1)^{p−m−n} z)` with the
  Gamma coefficients `A_k`
- [x] 1.4 Wire into `hyperexpand` (reduce the Meijer-G, then expand each hyper)
- [x] 1.5 Tests vs `sympy.hyperexpand`: exp, power·exp, rational, algebraic,
  two-pole sinh/cosh, and the confluent-case-left-opaque guard
- [x] 1.6 Docs: header comment, README, parity roadmap

## Phase 2 — confluent / logarithmic case

Reassessed: the confluent (integer-spaced lower parameter) Meijer-G functions are
**special-function-valued** (modified Bessel `K`, Bessel `Y`, …), not elementary.
SymPy's own `hyperexpand` likewise returns these Meijer-G nodes **unchanged**, so
SymPP is already at parity by detecting the confluent case and leaving it opaque
(Phase 1 guard, tested). A true closed form requires first introducing the
Bessel-family special functions as represented objects.

- [x] 2.1 Detect integer-spaced lower parameters
- [x] 2.2 Bessel-family special functions (besselj/y/i/k) — shipped
- [x] 2.3 Confluent G^{2,0}_{0,2} → 2·z^{(b1+b2)/2}·K_{b1−b2}(2√z) closed form

## Phase 4 — function → Meijer-G recognition

- [x] 4.1 `to_meijerg(f, var)` starter table (eᵍ, varᵃ·e^{−var}, 1/(1+var)),
  round-tripping through `hyperexpand`
- [x] 4.2 `meijerg_integrate_0_inf_of(f, var)` — ∫₀^∞ via recognition + Mellin
  (e.g. ∫₀^∞ xᵃe^{−x} = Γ(a+1))
- [x] 4.3a Dispatch the general `integrate(..., 0, ∞)` through Meijer-G (done in 3.4)
- [x] 4.3b `η·xᶜ` substitution + exp/sin/cos table → Gaussian / Dirichlet / Fresnel
  integrals, with per-core convergence guards (oscillatory strip)
- [x] 4.3c Bessel cores — confluent G^{2,0}_{0,2} → modified Bessel K (Phase 2,
  via the shipped Bessel family); erf-core ∫₀^∞ already covered (ERFC-INT-1)

## Phase 3 — Mellin–Barnes definite integration

- [x] 3.1 Mellin transform `M[G(η·x)](s)` master formula (`meijerg_mellin_transform`)
- [x] 3.2 `∫₀^∞ G(η·x) dx` = transform at s=1, with convergence guard
  (`meijerg_integrate_0_inf`)
- [x] 3.3 Product `∫₀^∞ G₁·G₂` convolution (Mellin–Parseval, `meijerg_convolution`)
- [x] 3.4 Wire into the general `integrate(expr, var, 0, ∞)` (Meijer-G rule, with an
  η>0 convergence guard so generic-scale integrals stay unevaluated)

## Phase 4 — function → Meijer-G recognition (done; see the section above)

- [x] 4.1 Rewrite elementary functions into Meijer-G form (`to_meijerg`,
  `recognize_form`)
- [x] 4.2 Integrator dispatch through the Meijer-G representation
  (`meijerg_integrate_0_inf_of`, wired into `integrate(·,0,∞)`)
