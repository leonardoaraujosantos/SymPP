## Context

The limit engine was a layered set of heuristics (growth-rank table, L'Hôpital,
reciprocal substitution). Converting it to a Gruntz-style engine in one step would
have been a high-risk rewrite, so this change adds **composable leading-term
stages** that are each gated to the nan-branch (only fire when the existing
machinery returns `nan`), keeping every increment zero-regression. Several stages
compose: the log-exp reduction emits a sum that the dominant-term rule and power
continuity then resolve.

## Key decisions

- **Soundness via numeric verification, not symbolic series.** SymPP has no full
  asymptotic-series engine, so stages that drop higher-order terms (small-angle
  substitution, harmonic asymptotics) verify the candidate by sampling the
  original at large `|x|`. The sample **precision scales with the point** so a
  difference like `sin(a+h) − sin(a)` with `h ∼ e⁻ˣ` is resolved rather than lost
  to catastrophic cancellation; the sample points stay moderate (≤ 600) so a big
  exponential does not overflow.
- **Gamma multiplication routes through the existing asymptotic.** Legendre/Gauss
  multiplication rewrites `Γ(kx+b)` into k slope-1 gammas plus a constant-base
  exponential, then resolves *directly* through `gamma_ratio_asymptotic` (never a
  recursive `limit`) — this is what stops the previous `Γ(2x)/Γ(x)` ~80 s case and
  the `Γ(2n)·n⁻ⁿ` loop.
- **Assuming context injected at one hook.** Every query funnels through
  `ask(e, k)`; the context is consulted in its `direct()` helper, so a scoped fact
  composes with the full implication chain. Compound nodes (`Mul`/`Add`/`Pow`)
  route their child queries through an exported `direct_ask` so the context
  propagates — and `direct_ask ≡ e->ask(k)` when no scope is active, making the
  change zero-regression by construction.
- **Incomplete gamma mirrors the existing Beta/PolyGamma two-arg pattern**; the
  half-integer forms reuse the existing `erf`/`erfc` via the recurrence
  `V(s+1) = s·V(s) ± xˢ·e⁻ˣ`.

## Risks / trade-offs

- The leading-term substitution can leave a slow-converging case (e.g.
  `(x^(1/x) − 1)·x/log x`) as honest `nan` when the numeric check cannot certify
  it within tolerance; correctness is preserved (no wrong answers), parity is not
  fully reached on those.
- `∂/∂s` of the incomplete gammas is non-elementary and returns `0` (the project's
  polygamma-order convention), not the Meijer-G form SymPy gives.
