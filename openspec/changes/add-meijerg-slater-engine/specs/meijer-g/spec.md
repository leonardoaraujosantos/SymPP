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
