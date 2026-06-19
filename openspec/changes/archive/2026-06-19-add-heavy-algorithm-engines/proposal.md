## Why

SymPP reached textbook-shaped parity with SymPy, but a band of *research-grade*
algorithms remained heuristic-only: the limit engine ranked a fixed growth
hierarchy instead of computing asymptotics, the incomplete gamma functions were
absent, and assumptions could only be attached to a symbol at construction — there
was no way to assert a fact (`x > 0`) for a region of code and have `refine`/`ask`
respect it. This change converts those heuristics into real named algorithms so
the hard asymptotics (`Γ(2n)/nⁿ`, the Gruntz oscillation), the incomplete-gamma
integral family, and relational `assuming` simplification reach SymPy parity.

## What Changes

- **Gruntz-style limit asymptotics** — the heuristic `limit` gains composable
  leading-term stages: power-as-exp rewrite, dominant-term `∞−∞` resolution,
  constant-exponent power continuity, log-exp reduction for nested
  transcendentals, small-angle / Maclaurin-head substitution for `0·∞` trig and
  exp/cos forms (numerically verified), harmonic-number asymptotics, and a full
  gamma-ratio asymptotic engine (Stirling slope-1, Legendre/Gauss multiplication,
  constant-base exponential rates, unbalanced gamma powers). Two non-terminating
  inputs were also made to terminate.
- **Incomplete gamma functions** — new `lowergamma`/`uppergamma` real function
  classes (FunctionId, class with rebuild/ask/diff, builders, parser registration)
  with positive-integer and half-integer (erf/erfc) closed forms; wired into the
  integrator (`∫xˢ⁻¹e⁻ˣ = γ(s,x)`, the Gamma integral `∫₀^∞ = Γ(s)`).
- **Relational `assuming` context** — a thread-local scoped assumption stack that
  the central `ask` consults, propagating scoped facts through `Mul`/`Add`/`Pow`,
  plus assumption-aware `refine` rules (`√(x²) → |x|`, `|x·y|` distribution).
- Several correctness/robustness fixes uncovered along the way: malformed
  `∞`-arithmetic limit output, the infinite geometric series for negative/symbolic
  ratios, the generalized binomial for numeric upper index, `polygamma` at `∞`.

## Capabilities

### New Capabilities
- `gruntz-limits`: most-rapidly-varying / leading-term limit asymptotics at `±∞`,
  covering the gamma-ratio family, nested transcendentals, `0·∞` analytic-head
  forms, and harmonic numbers — with numeric verification where a substitution
  drops higher-order terms.
- `incomplete-gamma`: `lowergamma`/`uppergamma` as first-class functions
  (eval/diff/ask/str/parser) and their use as closed forms in integration.
- `assuming-context`: scoped relational assumptions (`assuming(x, positive)`) that
  drive `ask`/`refine` globally and retract at scope exit.

### Modified Capabilities
<!-- No prior openspec/specs/ existed before this change; these are all new. -->

## Impact

- `src/calculus/limit.cpp` (the bulk — new asymptotic stages + verification),
  `src/calculus/summation.cpp` (geometric/exp-series fixes),
  `src/functions/combinatorial.cpp` (incomplete gamma, binomial, polygamma),
  `src/integrals/integrate.cpp` (incomplete-gamma rules),
  `src/core/{assumption_context,queries,refine}.{hpp,cpp}` (assuming context),
  `src/core/{mul,add,pow}.cpp` (context-aware child queries),
  `include/sympp/core/function_id.hpp` (new FunctionIds), `src/parsing/parser.cpp`.
- No public-API breakage; all additions are new functions/classes. Regression
  coverage in `tests/calculus/series_limit_test.cpp`, `tests/integrals/`,
  `tests/functions/combinatorial_test.cpp`, `tests/core/assumptions_test.cpp`.
