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

## Phase 2 — confluent / logarithmic case (later session)

- [ ] 2.1 Detect integer-spaced lower parameters; group into pole clusters
- [ ] 2.2 Derivative-of-Gamma (digamma) expansion producing `log`-carrying terms
- [ ] 2.3 Tests vs SymPy for Bessel-Y / log-bearing Meijer-G

## Phase 3 — Mellin–Barnes definite integration (later session)

- [ ] 3.1 `∫₀^∞ G(ax) dx` master formula
- [ ] 3.2 Product `∫₀^∞ G₁·G₂` convolution
- [ ] 3.3 Wire into `integrate` for the `meijerint` family

## Phase 4 — function → Meijer-G recognition (later session)

- [ ] 4.1 Rewrite elementary/special functions into Meijer-G form
- [ ] 4.2 Integrator dispatch through the Meijer-G representation
