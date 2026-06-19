# gruntz-limits Specification

## Purpose

Evaluate hard asymptotic limits at `±∞` and at finite points that the bare
substitution / L'Hôpital heuristics leave as `nan` or spin on, using a composable
stack of Gruntz-style named asymptotic stages (most-rapidly-varying rewrites,
dominant-term selection, special-function asymptotic series, and leading-term-by-
series), always terminating with a determinate value or an honest `nan`.

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
to a unit term — Gruntz's flagship most-rapidly-varying example. It SHALL also
resolve such a difference wrapped in a product, `M·(c·eᵖ − c·eᵍ)` with `p − q → 0`,
via the exact equivalence `eᵖ − eᵍ ~ eᵍ·(p − q)` (a bare constant counts as `e⁰`,
so `eᵘ − 1` is the same rule), read off algebraically without a numeric guard
(which would underflow on the `e^{−e⁻ˣ} − 1` cancellation).

#### Scenario: Difference of asymptotically-equal exponentials
- **WHEN** `limit(e^{x+e⁻ˣ} − eˣ, x, ∞)` is requested
- **THEN** the result is `1`

#### Scenario: Product times a unit exponential difference
- **WHEN** `limit(eˣ·(e^{1/x − e⁻ˣ} − e^{1/x}), x, ∞)`, `limit(eˣ·(e^{1/x} − 1), x, ∞)` are requested
- **THEN** the results are `−1` and `∞`

### Requirement: Special-function asymptotic series

The engine SHALL expand a special function of a diverging argument by its named
asymptotic series so that limits where a leading term has been cancelled resolve to
a determinate value: the complementary error function `erfc`/`erf` (Gaussian tail,
DLMF 7.12.1), the exponential integral `Ei` (DLMF 6.12.2), the Riemann zeta `ζ(s)`
(`→ 1` with `2⁻ˢ` correction), the digamma and higher polygamma (DLMF 5.15.8), the
log-Stirling `loggamma`, the Tricomi–Erdélyi gamma-ratio subleading term, and the
arctangent / arccotangent (`atan(g) ~ ±π/2 − 1/g + 1/(3g³) − …`, DLMF 4.24.3).

#### Scenario: Complementary error function and exponential integral
- **WHEN** `limit(x·eˣ²·erfc(x), x, ∞)` and `limit(x·e⁻ˣ·Ei(x), x, ∞)` are requested
- **THEN** the results are `1/√π` and `1`

#### Scenario: Riemann zeta tends to 1
- **WHEN** `limit(ζ(x), x, ∞)` and `limit(2ˣ·(ζ(x) − 1), x, ∞)` are requested
- **THEN** the results are `1` and `1`

#### Scenario: Arctangent subleading asymptotics
- **WHEN** `limit(x·(atan x − π/2), x, ∞)`, `limit(x³·(atan x − π/2 + 1/x), x, ∞)`, and the expanded `limit(x·atan x − π·x/2, x, ∞)` are requested
- **THEN** the results are `−1`, `1/3`, and `−1`

### Requirement: Logarithm of a dominated sum

The engine SHALL extract the dominant summand of a `log` of a var-dependent sum,
`log(Σ tᵢ) = log(t*) + log(Σ/t*)`, for both an exponential-dominated sum (treating
a polynomial summand as a sub-exponential, rate-0 coefficient) and a sum whose
dominant summand has logarithmic rate, so the residual `log(1 + small)` vanishes
and ratios/differences of such logs resolve.

#### Scenario: Ratio of logs of exponential-dominated polynomial sums
- **WHEN** `limit(log(x + eˣ)/log(x + e²ˣ), x, ∞)` and `limit(log(x³ + 5eˣ)/x, x, ∞)` are requested
- **THEN** the results are `1/2` and `1`

#### Scenario: Log of a logarithmic-rate sum
- **WHEN** `limit(log(log x + log log x) − log log x, x, ∞)` is requested
- **THEN** the result is `0`

### Requirement: Series core times a non-polynomial multiplier

The engine SHALL resolve a vanishing `1^∞` series core multiplied by a factor the
direct coordinate substitution cannot expand (`eˣ → exp(1/u)`, `√x`/`log x` poles)
by peeling the core's leading monomial `c₀·x⁻ᵐ` and re-taking `M·c₀·x⁻ᵐ`.

#### Scenario: Series core times an exponential / fractional-power multiplier
- **WHEN** `limit(eˣ·((1+1/x)ˣ − e), x, ∞)`, `limit(√x·((1+1/x)ˣ − e), x, ∞)` are requested
- **THEN** the results are `−∞` and `0`

### Requirement: Tower products via log-exp reduction

The engine SHALL resolve a product of variable-exponent powers (a "tower product"
substituting to the indeterminate `∞·0`) by `lim e = exp(lim log e)`, hoisted ahead
of the `f(x)^{g(x)}` power stages that otherwise spin on it.

#### Scenario: x^x / (x+1)^(x+1)
- **WHEN** `limit(xˣ/(x+1)^{x+1}, x, ∞)` and `limit((x+1)^{x+1}/xˣ, x, ∞)` are requested
- **THEN** the results are `0` and `∞`

### Requirement: Different-rate radical differences resolve by dominant term

For an `∞−∞` difference whose terms have different growth rates, the engine SHALL
apply the dominant-term rule before the conjugate-rationalization stage, so a
radical difference of unequal rate (`log log x − √(log x·log log x)`) follows its
dominant term rather than recursing on a worse `∞−∞`. Equal-rate radical
differences SHALL still resolve via the conjugate.

#### Scenario: Different-rate radical difference
- **WHEN** `limit(log log x − √(log x·log log x), x, ∞)` is requested
- **THEN** the result is `−∞`, resolved in milliseconds (not by tens of seconds of conjugate recursion)

#### Scenario: Equal-rate radical difference unchanged
- **WHEN** `limit(√(x²+x) − x, x, ∞)` is requested
- **THEN** the result is `1/2` (still the conjugate path)

