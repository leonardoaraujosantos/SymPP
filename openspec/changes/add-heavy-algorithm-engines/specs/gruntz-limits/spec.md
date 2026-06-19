## ADDED Requirements

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
