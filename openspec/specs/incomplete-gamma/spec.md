# incomplete-gamma Specification

## Purpose
TBD - created by archiving change add-heavy-algorithm-engines. Update Purpose after archive.
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

