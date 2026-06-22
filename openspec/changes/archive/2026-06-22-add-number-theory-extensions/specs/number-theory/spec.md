## ADDED Requirements

### Requirement: Chinese Remainder Theorem

SymPP SHALL provide `crt(moduli, residues)` that returns the unique solution
`x` modulo `M = lcm(moduli)` of the system `x ≡ residues[i] (mod moduli[i])`,
normalized to `0 ≤ x < M`, or no solution when the congruences are inconsistent.
Moduli need NOT be pairwise coprime.

#### Scenario: Coprime moduli
- **WHEN** `crt([3, 5, 7], [2, 3, 2])` is evaluated
- **THEN** the result is `(23, 105)` and `23 mod 3 = 2`, `23 mod 5 = 3`, `23 mod 7 = 2`

#### Scenario: Non-coprime but consistent moduli
- **WHEN** `crt([4, 6], [0, 2])` is evaluated
- **THEN** the result is `(8, 12)` with `12 = lcm(4, 6)`

#### Scenario: Inconsistent congruences
- **WHEN** `crt([6, 10], [2, 3])` is evaluated
- **THEN** no solution is returned (`std::nullopt`)

### Requirement: Discrete logarithm

SymPP SHALL provide `discrete_log(n, a, b)` returning the least non-negative
integer `x` with `bˣ ≡ a (mod n)`, or no solution when none exists. The search
SHALL be baby-step/giant-step (`O(√n)`).

#### Scenario: Solvable discrete log
- **WHEN** `discrete_log(41, 15, 7)` is evaluated
- **THEN** the result is `3` (since `7³ ≡ 15 (mod 41)`)

#### Scenario: Verified against SymPy
- **WHEN** `discrete_log(23, 10, 5)` is evaluated
- **THEN** the result equals `sympy.ntheory.discrete_log(23, 10, 5)` (= `3`)

#### Scenario: No solution
- **WHEN** `discrete_log(n, a, b)` has no exponent satisfying the congruence
- **THEN** no solution is returned (`std::nullopt`)

### Requirement: Linear Diophantine equations

SymPP SHALL provide `diop_linear(a, b, c)` solving `a·x + b·y = c` over the
integers, returning a parametric solution family or no solution when
`gcd(a, b)` does not divide `c`.

#### Scenario: Solvable equation
- **WHEN** `diop_linear(2, 3, 5)` is evaluated
- **THEN** a parametric solution `(x₀, y₀, dx, dy)` is returned such that
  `a·(x₀ + dx·t) + b·(y₀ + dy·t) = c` for every integer `t`

#### Scenario: No solution when gcd does not divide c
- **WHEN** `diop_linear(2, 4, 5)` is evaluated
- **THEN** no solution is returned (`std::nullopt`), since `gcd(2, 4) = 2 ∤ 5`
