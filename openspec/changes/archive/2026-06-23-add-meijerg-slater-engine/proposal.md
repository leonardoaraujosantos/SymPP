## Why

The Meijer-G function is SymPy's universal special-function representation and
the engine behind its general definite-integration (`meijerint`) and
`hyperexpand`. SymPP currently keeps `meijerg(...)` **opaque** — it has the data
structure but no evaluation. Building a real Meijer-G engine is the single
largest remaining parity gap (it unlocks the long tail of definite integrals and
arbitrary hypergeometric closed forms), and it is genuinely multi-week research,
so it is run as a **multi-session initiative** with independently shippable
phases rather than one monolithic change.

## What Changes

This change tracks the whole initiative; phases land incrementally and each is a
green, oracle-validated step. The capability spec grows with each phase.

- **Phase 1 (this session) — generic Slater reduction.** `meijerg → Σ z^{b_k}
  ·pFq` via Slater's theorem in the generic case (no two lower `b₁…b_m` differ by
  an integer), wired into `hyperexpand` so Meijer-G that reduces to elementary
  closed forms now evaluates (e.g. `G^{1,1}_{1,1}`, `G^{1,0}_{0,1}`,
  `G^{2,0}_{0,2}` → exp / powers / sinh-cosh).
- **Phase 2 (later) — confluent/logarithmic case.** Lower parameters differing
  by integers (repeated poles) producing `log`-carrying expansions.
- **Phase 3 (later) — Mellin–Barnes definite integration.** `∫₀^∞ f(x) dx` via
  the Meijer-G master formula, the `meijerint` core.
- **Phase 4 (later) — function → Meijer-G recognition.** Rewrite elementary /
  special functions into Meijer-G form so the integrator can reach them.

## Capabilities

### New Capabilities
- `meijer-g`: evaluate the Meijer-G function by reducing it, via Slater's
  theorem, to a sum of generalized hypergeometric functions, then to closed
  form through `hyperexpand`. (Phase 1 covers the generic lower-parameter case.)

## Non-goals (initiative-wide, revisited per phase)

- Numeric Mellin–Barnes contour evaluation; the engine is symbolic.
- The fully general confluent case with arbitrary integer-spaced parameters in
  one step — staged across phases 1–2.
