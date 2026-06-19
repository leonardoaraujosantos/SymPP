# incomplete-gamma Specification

## Purpose

Provide the incomplete gamma functions `lowergamma`/`uppergamma` as real
two-argument function classes with their closed forms and exact derivatives, wire
them into integration, and — more broadly — guarantee that gamma-family special
functions return correct derivatives (a closed form where one exists, an honest
unevaluated `Derivative` where the direction is non-elementary), never a silent `0`.

## Requirements
### Requirement: lowergamma / uppergamma function classes

SymPP SHALL provide `lowergamma(s, x)` and `uppergamma(s, x)` as real
two-argument function classes — each with a `FunctionId`, `rebuild`, `ask`,
`diff_arg`, a builder factory, and parser registration — that reduce to closed
forms for a positive-integer or half-integer first argument and stay symbolic
otherwise.

#### Scenario: Positive-integer order closed form
- **WHEN** `uppergamma(2, x)` is evaluated
- **THEN** the result is `(x + 1)·e⁻ˣ`, and `lowergamma(2, x)` is `1 − (x + 1)·e⁻ˣ`

#### Scenario: Half-integer order via erf/erfc
- **WHEN** `uppergamma(1/2, x)` and `lowergamma(1/2, x)` are evaluated
- **THEN** the results are `√π·erfc(√x)` and `√π·erf(√x)`

#### Scenario: Special points
- **WHEN** `uppergamma(s, 0)`, `lowergamma(s, 0)`, `uppergamma(s, ∞)` are evaluated for a positive-integer order
- **THEN** the results are `Γ(s)`, `0`, `0`

#### Scenario: Exact x-derivatives, symbolic round-trip
- **WHEN** `∂/∂x γ(s, x)` is taken and `uppergamma(s, x)` is parsed from its `str()`
- **THEN** the derivative is `xˢ⁻¹·e⁻ˣ` and the parse round-trips to the original node

### Requirement: Incomplete gamma in integration

The integrator SHALL return incomplete-gamma closed forms for the
`xᵖ·e⁻ᶜˣ` family.

#### Scenario: Indefinite incomplete-gamma integral
- **WHEN** `∫ xˢ⁻¹·e⁻ˣ dx` is requested (symbolic exponent)
- **THEN** the result is `γ(s, x)`; a non-negative integer power keeps the
  elementary by-parts result

#### Scenario: The Gamma-function integral
- **WHEN** `∫₀^∞ xˢ⁻¹·e⁻ᶜˣ dx` is requested for a positive rate `c`
- **THEN** the result is `Γ(s)/cˢ`

### Requirement: Correct gamma-family parameter derivatives

Differentiating a gamma-family special function SHALL return its elementary closed
form where one exists, and an honest unevaluated `Derivative` for a non-elementary
direction (one needing a Meijer-G or symbolic-order Hurwitz zeta SymPP does not
represent) — never a silently-wrong `0`, which would drop the chain-rule term when
the order/parameter depends on the differentiation variable.

#### Scenario: Beta-function derivative closed form
- **WHEN** `∂/∂a Β(a, b)` is taken
- **THEN** the result is `Β(a, b)·(ψ(a) − ψ(a+b))` (symmetric in `b`)

#### Scenario: Harmonic-number derivative closed form
- **WHEN** `∂/∂x H(x)` is taken
- **THEN** the result is `polygamma(1, x + 1)` (the trigamma; the Hurwitz form `ζ(2, x+1)`)

#### Scenario: Non-elementary order/parameter derivative stays unevaluated
- **WHEN** `∂/∂s γ(s, x)`, `∂/∂s Γ(s, x)`, `∂/∂n ψ⁽ⁿ⁾(x)` are taken
- **THEN** each is an unevaluated `Derivative` of the function (not `0`)

