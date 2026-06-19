## ADDED Requirements

### Requirement: Scoped assumption context

SymPP SHALL provide a thread-local `assuming(target, AssumptionMask)` RAII scope
that augments the global assumption context. The central `ask(e, k)` SHALL consult
the active scopes (innermost first) so a scoped fact overrides the symbol's own
mask and feeds the implication chain, and the scope SHALL retract its facts when
it goes out of scope. When no scope is active, `ask` behavior SHALL be unchanged.

#### Scenario: Scoped fact and its consequences
- **WHEN** `assuming(x, positive)` is in scope
- **THEN** `is_positive(x)`, `is_real(x)`, `is_nonzero(x)` are all `true`, and they
  return to `unknown` after the scope exits

#### Scenario: Propagation through compound expressions
- **WHEN** `assuming(x, positive)` and `assuming(y, positive)` are in scope
- **THEN** `is_positive(x·y)`, `is_positive(x + y)`, `is_real(x³)`, `is_positive(√x)` are all `true`

### Requirement: Assumption-aware refine

`refine` SHALL apply assumption-justified rewrites that consult the active
assumption context.

#### Scenario: Sign-aware absolute value and square root
- **WHEN** `refine(|x|)` and `refine(√(x²))` are requested under a sign scope
- **THEN** under `positive` they give `x`, under `negative` they give `−x`, and
  under only `real` `√(x²)` gives `|x|`

#### Scenario: Absolute value distributes over a product
- **WHEN** `refine(|x·y|)` is requested
- **THEN** it is `x·y` under `x>0 ∧ y>0`, `x·|y|` under `x>0` alone, and unchanged
  with no sign facts
