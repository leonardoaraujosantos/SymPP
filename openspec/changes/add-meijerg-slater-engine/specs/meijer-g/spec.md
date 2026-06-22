## ADDED Requirements

### Requirement: Generic-case Slater reduction of Meijer-G

SymPP SHALL evaluate the Meijer-G function `G^{m,n}_{p,q}(z | a₁…aₚ; b₁…b_q)`
in the generic lower-parameter case — when no two of the lower poles
`b₁…b_m` differ by an integer — by Slater's theorem:

```
G = Σ_{k=1}^{m} A_k · z^{b_k} · pF_{q−1}( 1+b_k−a_j (j=1..p) ;
                                          1+b_k−b_j (j=1..q, j≠k) ;
                                          (−1)^{p−m−n} z )
A_k = [ Π_{j≠k}^{m} Γ(b_j−b_k) · Π_{j=1}^{n} Γ(1+b_k−a_j) ]
      / [ Π_{j=m+1}^{q} Γ(1+b_k−b_j) · Π_{j=n+1}^{p} Γ(a_j−b_k) ]
```

The reduction SHALL be invoked from `hyperexpand`, which then collapses each
resulting hypergeometric to a closed form where one is known. When two lower
parameters differ by an integer (the confluent case, a later phase), the
Meijer-G node SHALL be left unevaluated rather than returning a wrong value.

#### Scenario: Exponential — G^{1,0}_{0,1}
- **WHEN** `hyperexpand(meijerg([], [], [0], [], z))` is evaluated
- **THEN** the result is `e^{−z}`

#### Scenario: Power times exponential
- **WHEN** `hyperexpand(meijerg([], [], [1/2], [], z))` is evaluated
- **THEN** the result equals `√z·e^{−z}` (matching `sympy.hyperexpand`)

#### Scenario: Rational — G^{1,1}_{1,1}
- **WHEN** `hyperexpand(meijerg([0], [], [0], [], z))` is evaluated
- **THEN** the result is `1/(z+1)`

#### Scenario: Algebraic with a Gamma coefficient
- **WHEN** `hyperexpand(meijerg([1/2], [], [0], [], z))` is evaluated
- **THEN** the result equals `√π/√(z+1)`

#### Scenario: Two lower poles → hyperbolic combination
- **WHEN** `hyperexpand(meijerg([], [], [0, 1/2], [], z))` is evaluated
- **THEN** the result equals `√π·cosh(2√z) − √π·sinh(2√z)`

#### Scenario: Confluent case left unevaluated
- **WHEN** a Meijer-G whose lower parameters `b₁…b_m` contain two values
  differing by an integer is passed to `hyperexpand`
- **THEN** the Meijer-G node is returned unchanged (no wrong closed form)

### Requirement: Meijer-G Mellin transform and definite integral

SymPP SHALL provide the Mellin transform of a Meijer-G as the Gamma-ratio master
formula, and the definite integral `∫₀^∞` as that transform at `s = 1`, rejecting
divergent (Gamma-pole) results and non-Meijer / non-`η·var` inputs.

#### Scenario: Mellin transform
- **WHEN** `meijerg_mellin_transform(G^{1,0}_{0,1}([],[],[1/2],[], x), x, s)` is evaluated
- **THEN** the result is `Γ(s + 1/2)` (matching `sympy.mellin_transform`)

#### Scenario: Convergent definite integral
- **WHEN** `meijerg_integrate_0_inf(G^{1,0}_{0,1}([],[],[1/2],[], x), x)` is evaluated
- **THEN** the result is `√π/2`

#### Scenario: Divergent integral rejected
- **WHEN** `meijerg_integrate_0_inf(G^{1,1}_{1,1}([0],[],[0],[], x), x)` is evaluated
- **THEN** no value is returned (`std::nullopt`), since `∫₀^∞ 1/(1+x) dx` diverges

### Requirement: Function-to-Meijer-G recognition

SymPP SHALL recognize a starter table of elementary functions as Meijer-G forms
such that `hyperexpand(to_meijerg(f)) = f`, and route `∫₀^∞ f dx` through the
Mellin master formula for recognized `f`.

#### Scenario: Recognition round-trips
- **WHEN** `hyperexpand(to_meijerg(x²·e^{−x}, x))` is evaluated
- **THEN** the result is `x²·e^{−x}`

#### Scenario: Gamma integral via Meijer-G
- **WHEN** `meijerg_integrate_0_inf_of(xᵃ·e^{−x}, x)` is evaluated
- **THEN** the result is `Γ(a + 1)`

### Requirement: General integrate dispatches ∫₀^∞ through Meijer-G

The general `integrate(expr, var, 0, ∞)` SHALL route a Meijer-G integrand — or a
function recognized as one — through the Mellin master formula, but only when the
exponential-decay scale `η` is provably positive, leaving the integral
unevaluated when convergence depends on an unknown sign.

#### Scenario: Bare Meijer-G through the general integrator
- **WHEN** `integrate(G^{1,0}_{0,1}([],[],[1/2],[], x), x, 0, ∞)` is evaluated
- **THEN** the result is `√π/2`

#### Scenario: Recognized function through the general integrator
- **WHEN** `integrate(xᵃ·e^{−x}, x, 0, ∞)` is evaluated
- **THEN** the result is `Γ(a + 1)`

#### Scenario: Unknown-sign scale stays unevaluated
- **WHEN** `integrate(e^{−a·x}, x, 0, ∞)` is evaluated for a generic `a` (sign unknown)
- **THEN** the result is an unevaluated `Integral`, not `1/a`

### Requirement: Power-argument definite integration

SymPP SHALL evaluate `∫₀^∞ C·xᵈ·core(var) dx` for `core ∈ {exp, sin, cos}` with a
monomial argument `k·varᵉ` via the substitution
`∫₀^∞ C·xᵈ·G(η·xᶜ) dx = C·(1/c)·M_G((1+d)/c)`, gated by a convergence test: the
scale MUST be positive, and for oscillatory cores the strip
`−1 < (d+1)/e < 1` (sin) / `0 < (d+1)/e < 1` (cos) MUST hold. This routes the
Gaussian, Dirichlet and Fresnel integrals through Meijer-G.

#### Scenario: Gaussian integral
- **WHEN** `integrate(e^{−x²}, x, 0, ∞)` is evaluated
- **THEN** the result is `√π/2` (and `∫ x·e^{−x²} = 1/2`, `∫ x²·e^{−x²} = √π/4`)

#### Scenario: Dirichlet integral
- **WHEN** `integrate(sin(x)/x, x, 0, ∞)` is evaluated
- **THEN** the result is `π/2`

#### Scenario: Fresnel integrals
- **WHEN** `integrate(cos(x²), x, 0, ∞)` and `integrate(sin(x²), x, 0, ∞)` are evaluated
- **THEN** both are `√(2π)/4`

#### Scenario: Divergent oscillatory integral rejected
- **WHEN** `integrate(sin(x), x, 0, ∞)` is evaluated (oscillates, no limit)
- **THEN** the Meijer-G route abstains (strip condition fails), not returning a value
