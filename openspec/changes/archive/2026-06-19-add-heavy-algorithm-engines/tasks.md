# Tasks

All tasks are implemented on branch `feat/assumptions-and-engines`; each was a
single zero-regression increment, verified against the SymPy oracle with a
regression test and a green full suite. Known-issue tags in
`docs/09-known-issues.md` are noted in parentheses.

## 1. Gruntz-style limit asymptotics

- [x] 1.1 Power-as-exp rewrite for general powers `f(x)^g(x)` (LIMIT-POW-AS-EXP-1)
- [x] 1.2 Dominant-term resolution of unequal-order `∞−∞` sums (LIMIT-DOMINANT-SUM-1)
- [x] 1.3 Constant-exponent power continuity `lim bʳ = (lim b)ʳ` (LIMIT-POW-CONTINUITY-1)
- [x] 1.4 Log-exp reduction for nested transcendentals (LIMIT-LOG-EXP-REDUCTION-1)
- [x] 1.5 Buried-`∞` guard + `polygamma(∞)` values (LIMIT-BURIED-INF-1, POLYGAMMA-INF-1)
- [x] 1.6 Legendre duplication for the central binomial (LIMIT-GAMMA-DUP-1)
- [x] 1.7 Constant-base exponential rate in the gamma-ratio asymptotic (LIMIT-GAMMA-EXPRATE-1)
- [x] 1.8 Unbalanced gamma power → fixes the ~80 s `Γ(2x)/Γ(x)` case (LIMIT-GAMMA-UNBALANCED-1)
- [x] 1.9 Gauss multiplication for slope-k gammas (LIMIT-GAMMA-GAUSS-1)
- [x] 1.10 Harmonic-number asymptotics (LIMIT-HARMONIC-1)
- [x] 1.11 Small-angle substitution + the Gruntz oscillation (LIMIT-SMALL-ANGLE-1)
- [x] 1.12 Nested-radical conjugate differences (LIMIT-NESTED-RADICAL-1)
- [x] 1.13 Terminate `(x^(1/x)−1)·x/log x` instead of hanging (LIMIT-NOHANG-1; now resolves to `1`)
- [x] 1.14 Rational functions of constant-base exponentials `cˣ` (LIMIT-CONST-BASE-RATIO-1)
- [x] 1.15 Full Stirling prefactor for cancelling products, `n!/(nⁿe⁻ⁿ)` (LIMIT-STIRLING-PRODUCT-1)
- [x] 1.16 Flatten combined denominators in the gamma-ratio asymptotic (LIMIT-GAMMARATIO-1)
- [x] 1.17 Log-Stirling expansion for `log(n!)`, `log Γ(n)` (LIMIT-LOGGAMMA-1)
- [x] 1.18 Rewrite hyperbolics of a diverging argument to exp form (LIMIT-HYPERBOLIC-1)
- [x] 1.19 Gruntz leading-term-by-series for `1^∞` corrections (LIMIT-SERIES-1)
- [x] 1.20 Finite-point + transcendental-base/divergent-term series (LIMIT-SERIES-1)
- [x] 1.21 Unit-tending power difference via `exp(t)−1 ~ t` (LIMIT-UNIT-POWER-1)
- [x] 1.22 MRV rewrite for differences of exponentials, `e^{x+e⁻ˣ}−eˣ` (LIMIT-EXP-DIFF-1)
- [x] 1.23 Inverse-hyperbolic asymptotics `asinh/acosh ~ log(2g) ± 1/(4g²)` (LIMIT-INVHYP-1)

## 2. Incomplete gamma functions

- [x] 2.1 `lowergamma`/`uppergamma` function classes + positive-integer forms (FUNC-INCGAMMA-1)
- [x] 2.2 Indefinite integral `∫ xˢ⁻¹e⁻ˣ = γ(s,x)` (INT-INCGAMMA-1)
- [x] 2.3 The Gamma-function integral `∫₀^∞ xˢ⁻¹e⁻ᶜˣ = Γ(s)/cˢ` (INT-GAMMA-1)
- [x] 2.4 Half-integer erf/erfc closed forms (FUNC-INCGAMMA-HALF-1)

## 3. Relational assuming context

- [x] 3.1 Thread-local `assuming` RAII context + `refine_abs` (ASSUMING-CONTEXT-1)
- [x] 3.2 `refine(√(x²)) → |x|` gated on reality (REFINE-SQRT-SQUARE-1)
- [x] 3.3 `refine(|·|)` distributes over products (REFINE-ABS-MUL-1)
- [x] 3.4 Context propagates through `Mul`/`Add` queries (ASSUMING-PROPAGATE-1)
- [x] 3.5 Context propagates through `Pow` queries

## 4. Correctness fixes surfaced along the way

- [x] 4.1 Gamma-spelling of the exponential series `Σ xᵏ/Γ(k+1)` (SUM-EXP-GAMMA-1)
- [x] 4.2 Infinite geometric series for negative/symbolic ratios (SUM-GEOM-INF-1)
- [x] 4.3 Generalized binomial for numeric upper index (BINOM-GEN-1)

## 5. Validation

- [x] 5.1 Each increment has a regression test and oracle cross-check
- [x] 5.2 Full suite green (5373 assertions / 1534 cases) after every commit
