# gruntz-limits Specification

## Purpose
TBD - created by archiving change add-heavy-algorithm-engines. Update Purpose after archive.
## Requirements
### Requirement: Gamma-ratio asymptotics at infinity

The limit engine SHALL evaluate limits of products of gamma/factorial factors,
power factors `x^q`, and constant-base exponentials `c^(k·x)` at `+∞` using the
Stirling asymptotic `Γ(x+a)/Γ(x) ~ x^a`, returning a determinate value (a finite
constant, `0`, or signed `±∞`) rather than `nan`.

#### Scenario: Balanced same-rate gamma ratio
- **WHEN** `limit(Γ(x+1/2)/Γ(x)/√x, x, ∞)` is requested
- **THEN** the result is `1`

#### Scenario: Doubled-rate gamma via Legendre/Gauss multiplication
- **WHEN** `limit(Γ(2n)/Γ(n)², n, ∞)` or `limit(Γ(3n)/Γ(n)³, n, ∞)` is requested
- **THEN** the gamma is split by the multiplication formula and the result is `∞`

#### Scenario: Central binomial over its scale
- **WHEN** `limit(Γ(2x+1)/Γ(x+1)²/4ˣ, x, ∞)` is requested
- **THEN** the result is `0`, and `·√x` gives `1/√π`

#### Scenario: Exponential rate dominates a gamma ratio
- **WHEN** `limit(2ˣ·Γ(x+1/2)/Γ(x), x, ∞)` is requested
- **THEN** the result is `∞` (the exponential rate dominates the polynomial part)

#### Scenario: Unbalanced gamma power
- **WHEN** `limit(Γ(2x)/Γ(x), x, ∞)` is requested
- **THEN** the result is `∞`, computed without falling into the slow gamma-growth
  fallback

### Requirement: Leading-term substitution for indeterminate products

The limit engine SHALL resolve an indeterminate zero-times-infinity product at
`±∞` by replacing an analytic function of a vanishing argument with its Maclaurin
head (`sin g, tan g, … → g`; `eᵍ → 1+g+g²/2`; `cos g → 1−g²/2`), re-taking the
limit, and accepting the result only after a numeric check against the original
confirms it.

#### Scenario: Canonical Gruntz oscillation
- **WHEN** `limit(eˣ·(sin(1/x + e⁻ˣ) − sin(1/x)), x, ∞)` is requested
- **THEN** the result is `1`

#### Scenario: Exponential and cosine heads
- **WHEN** `limit(eˣ·(e^{e⁻ˣ} − 1), x, ∞)` and `limit(eˣ·(cos e⁻ˣ − 1), x, ∞)` are requested
- **THEN** the results are `1` and `0`

#### Scenario: A non-vanishing argument is not substituted
- **WHEN** `limit(x·sin x, x, ∞)` is requested
- **THEN** the result is `nan` (the argument does not tend to 0)

### Requirement: Nested-transcendental and harmonic asymptotics

The engine SHALL resolve positive product/power forms whose growth ranking is
opaque via the log-exp reduction `lim e = exp(lim log e)`, and SHALL expand a
harmonic number `H(g)` with `g → +∞` to `log g + γ + 1/(2g) − 1/(12g²)`.

#### Scenario: Nested transcendental ratio
- **WHEN** `limit(log x / exp(√(log x · log log x)), x, ∞)` is requested
- **THEN** the result is `0`

#### Scenario: Harmonic number limits
- **WHEN** `limit(H(n), n, ∞)`, `limit(H(n)/log n, n, ∞)`, `limit(H(n) − log n, n, ∞)` are requested
- **THEN** the results are `∞`, `1`, `γ` respectively

### Requirement: Robust termination

The limit engine SHALL terminate (returning a determinate value or `nan`) on
inputs that previously did not, and SHALL never return an `∞` buried inside an
arithmetic expression.

#### Scenario: General power inside a difference
- **WHEN** `limit((x^(1/x) − 1)·x/log x, x, ∞)` is requested
- **THEN** the call terminates with a determinate value (does not hang)

#### Scenario: No ∞-arithmetic noise
- **WHEN** a divergent gamma ratio drives L'Hôpital to an `∞·(∞+…)⁻¹` form
- **THEN** the engine returns a clean `±∞` or `nan`, never the buried-`∞` expression

### Requirement: Rational functions of constant-base exponentials

The limit engine SHALL evaluate rational functions of constant-base
exponentials `cˣ` (c a positive number) at `±∞` by dividing through by the
denominator's dominant term, resolving forms that L'Hôpital loops on and that
per-factor folding mis-reads.

#### Scenario: Constant-base exponential ratio
- **WHEN** `limit((2ˣ+3ˣ)/3ˣ, x, ∞)`, `limit((2ˣ+1)/(2ˣ−1), x, ∞)` are requested
- **THEN** the results are `1` and `1`

### Requirement: Full Stirling and log-Stirling asymptotics

The engine SHALL use the full leading Stirling form `n! ~ √(2πn)·(n/e)ⁿ` for
factorial/gamma products whose `(n/e)ⁿ` factor cancels, and the log-Stirling
expansion `log Γ(z) ~ (z−½)·log z − z + ½·log 2π` for a log of a divergent
factorial/gamma. Gamma ratios written as ordinary division (a combined
denominator) SHALL be flattened so the asymptotic still applies.

#### Scenario: Stirling-prefactor product
- **WHEN** `limit(n!/(nⁿ·e⁻ⁿ), n, ∞)` and `limit(n!·eⁿ·n^(−n−1/2), n, ∞)` are requested
- **THEN** the results are `∞` and `√(2π)`

#### Scenario: Log-Stirling
- **WHEN** `limit(log(n!)/(n·log n), n, ∞)` and `limit(log(n!)/n, n, ∞)` are requested
- **THEN** the results are `1` and `∞`

#### Scenario: Combined-denominator gamma ratio
- **WHEN** `limit(Γ(n+½)/(√n·Γ(n)), n, ∞)` is requested
- **THEN** the result is `1`

### Requirement: Gruntz leading-term-by-series for unit-tending powers

The engine SHALL resolve the asymptotic correction of a `1^∞` power inside a
difference, at `±∞` and at a finite point, by substituting the most-rapidly-
varying coordinate, pre-expanding the exponent to a polynomial, and reading the
limit off the leading series term (numerically verified). It SHALL use the
exponential expansion `f^g − 1 = exp(t) − 1 = t + t²/2 + …` for a power tending
to 1 whose base tends to `0` or `∞`.

#### Scenario: Asymptotic correction of (1+a/x)^x
- **WHEN** `limit(x·((1+1/x)ˣ − e), x, ∞)` and `limit(((1+x)^(1/x) − e)/x, x, 0)` are requested
- **THEN** the results are `−e/2` and `−e/2`

#### Scenario: Unit-tending power with a log-singular exponent
- **WHEN** `limit((xˣ − 1)/(x·log x), x, 0)` and `limit((x^(1/x) − 1)/(log x/x), x, ∞)` are requested
- **THEN** both results are `1`

### Requirement: Hyperbolic and inverse-hyperbolic asymptotics

The engine SHALL rewrite hyperbolic functions of a diverging argument to their
exponential form (resolving `sinh`/`cosh` combinations) and SHALL use the
two-term asymptotics `asinh(g), acosh(g) = log(2g) ± 1/(4g²) + O(g⁻⁴)` for an
inverse hyperbolic of a diverging argument (numerically verified).

#### Scenario: Hyperbolic combination
- **WHEN** `limit((sinh x + cosh x)/eˣ, x, ∞)` is requested
- **THEN** the result is `1`

#### Scenario: Inverse-hyperbolic over log
- **WHEN** `limit(asinh(x)/log x, x, ∞)` and `limit((acosh x − asinh x)·x², x, ∞)` are requested
- **THEN** the results are `1` and `−1/2`

### Requirement: MRV rewrite for differences of exponentials

The engine SHALL resolve a difference of asymptotically-equal exponentials by
factoring out a common `exp(b)` (an exact identity), collapsing the difference
to a unit term — Gruntz's flagship most-rapidly-varying example.

#### Scenario: Difference of asymptotically-equal exponentials
- **WHEN** `limit(e^{x+e⁻ˣ} − eˣ, x, ∞)` is requested
- **THEN** the result is `1`

