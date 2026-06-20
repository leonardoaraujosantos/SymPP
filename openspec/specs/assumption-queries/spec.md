# assumption-queries Specification

## Purpose

Bring SymPP's assumption query layer to SymPy parity on the predicate vocabulary:
add the `commutative` predicate (default-true, structural) and a boolean /
SAT-style `ask` that resolves arbitrary `&&` / `||` / `!` combinations of
predicate literals — including exhaustive disjunctions that no single literal can
decide — by three-valued evaluation plus refutation against the assumption
closure.

## Requirements
### Requirement: Extended-signed predicates

SymPP SHALL expose `extended_positive`, `extended_negative`,
`extended_nonnegative` and `extended_nonpositive` predicates that describe the
sign of a value on the extended real line ℝ ∪ {±∞}. A strict sign SHALL imply its
extended counterpart (`positive ⇒ extended_positive`); the extended signs SHALL
NOT imply finiteness, and `extended_positive ∧ finite ⇒ positive`. They SHALL be
declarable via `AssumptionMask::set_extended_*` and queryable via
`is_extended_positive` etc.

#### Scenario: Strict sign implies extended sign
- **WHEN** `p` is positive and `n` is negative
- **THEN** `is_extended_positive(p)`, `is_extended_nonnegative(p)`,
  `is_extended_negative(n)` and `is_extended_nonpositive(n)` are all `true`

#### Scenario: Infinities carry an extended sign
- **WHEN** `is_extended_positive(+∞)` and `is_extended_negative(−∞)` are queried
- **THEN** both are `true` and `is_extended_real(+∞)` is `true`

#### Scenario: Extended-positive refines to positive under finiteness
- **WHEN** a symbol is declared `extended_positive`
- **THEN** it is `extended_real` and `nonzero`, `is_positive` is `unknown`, and
  adding `finite` makes `is_positive` `true`

### Requirement: Hermitian and antihermitian predicates

SymPP SHALL expose `hermitian` and `antihermitian` predicates with
`real ⇒ hermitian`, `imaginary ⇒ antihermitian`, and
`hermitian ∧ antihermitian ⇒ zero`.

#### Scenario: Reality and imaginarity imply conjugation symmetry
- **WHEN** `r` is real and `m` is imaginary
- **THEN** `is_hermitian(r)` is `true` and `is_antihermitian(m)` is `true`

#### Scenario: Both symmetries force zero
- **WHEN** a symbol is declared both `hermitian` and `antihermitian`
- **THEN** `is_zero` of that symbol is `true`

### Requirement: Commutative predicate

SymPP SHALL expose an `AssumptionKey::Commutative` predicate, declarable via
`AssumptionMask::set_commutative` and queryable via `is_commutative`. It SHALL
default to `true` for every number and symbol (matching SymPy) and SHALL be
`false` for any expression that contains a subexpression explicitly declared
`commutative=false`.

#### Scenario: Default commutativity
- **WHEN** `is_commutative(x)` is queried for an ordinary symbol `x` or any number
- **THEN** the result is `true`

#### Scenario: Non-commutative propagation
- **WHEN** `a` is declared `commutative=false` and `c` is an ordinary symbol
- **THEN** `is_commutative(a)`, `is_commutative(a·c)`, `is_commutative(a + c)` and
  `is_commutative(a²)` are all `false`, while `is_commutative(c·3)` is `true`

### Requirement: Boolean / SAT-style ask

SymPP SHALL provide a `Proposition` type built from predicate literals
(`Q(e, k)` / `Qn(e, k)`) combined with `!`, `&&` and `||`, and an
`ask(Proposition)` that returns `true` / `false` / `unknown`. It SHALL evaluate
the proposition three-valued (Kleene) and, when the literals all concern a single
expression, SHALL upgrade the answer by refutation: a proposition is `true` when
asserting its negation yields a logically inconsistent assumption mask, and
`false` when asserting the proposition itself does.

#### Scenario: Exhaustive disjunction over one expression
- **WHEN** `x` is real and nonzero and `ask(Q(x, Positive) || Q(x, Negative))` is asked
- **THEN** the result is `true`, even though neither `is_positive(x)` nor
  `is_negative(x)` is individually decidable

#### Scenario: Trichotomy and parity completeness
- **WHEN** `y` is real and `k` is an integer
- **THEN** `ask(Q(y, Positive) || Q(y, Negative) || Q(y, Zero))` and
  `ask(Q(k, Even) || Q(k, Odd))` are both `true`

#### Scenario: Contradictions are refuted
- **WHEN** `ask(Q(x, Positive) && Q(x, Negative))` is asked for any real `x`
- **THEN** the result is `false`

#### Scenario: Genuine ignorance stays unknown
- **WHEN** `x` is a generic symbol with no assumptions
- **THEN** `ask(Q(x, Positive) || Q(x, Negative))` is `unknown`
