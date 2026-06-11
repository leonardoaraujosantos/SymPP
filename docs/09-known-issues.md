# Known issues & SymPy-parity bug log

This file tracks correctness gaps found by cross-checking SymPP against
the SymPy 1.14 oracle. Each entry records the failing input, the
SymPy-expected result, current status, and the regression test (once
fixed). When `api.github.com` is reachable, each open entry should also
be filed as a GitHub issue; the local log is the source of truth in the
meantime.

## Fixed

### SUM-1 — geometric summation dropped non-trivial exponents
- **Input:** `summation(2**(-k), k, 0, n)`, `summation(2**(2*k), k, 0, 3)`,
  `summation(1/2**k, k, 0, n)`
- **Was:** returned the summand unchanged (e.g. `2**(-k)`), leaking the
  bound variable `k` into the result.
- **Expected (SymPy):** `2 - 2*2**(-n - 1)`, `85`, `2 - 2*2**(-n - 1)`.
- **Cause:** the geometric branch only matched when the exponent was
  *exactly* the summation variable, so `base^(c*var + d)` forms
  (negated / scaled / shifted index) fell through.
- **Fix:** detect any exponent linear in the index via `diff`, rewriting
  `base^(c*var + d)` as `base^d * (base^c)^var`.
- **Regression test:** `tests/calculus/series_limit_test.cpp`
  — `[summation][regression]` cases (Σ 2^(-k), Σ 2^(-k) 0..3, Σ 2^(2k)).
- **Commit:** see git log for `fix(calculus): geometric summation …`.

## Open

### LIM-1 — limit of the classic `e` definition returns garbage
- **Input:** `limit((1 + 1/x)**x, x, oo)`
- **Is:** `(oo**(-1) + 1)**oo` (unevaluated, nonsensical).
- **Expected (SymPy):** `E`.
- **Notes:** the limit engine substitutes the target naively instead of
  detecting the `1^∞` indeterminate form (needs `exp(limit(x*log(1+1/x)))`
  rewrite). Partly blocked by the missing `Infinity` singleton — `oo` is
  currently parsed as a plain symbol.
- **Status:** open. Regression test to be added with the fix.

### SIMP-1 — `simplify` can return a *more complex* expression
- **Input:** `simplify((x**2 - 1)/(x - 1))`
- **Is:** `((x - 1)*x**2 - (x - 1))*(x - 1)**(-2)` (worse than input).
- **Expected (SymPy):** `x + 1`.
- **Notes:** rational-function cancellation (`cancel`) is not invoked /
  not effective inside the `simplify` orchestrator for this shape; the
  result should never be structurally larger than the input.
- **Status:** open. Regression test to be added with the fix.

### Missing auto-simplifications (lower priority)
- `exp(log(x))` does not reduce to `x` (SymPy auto-evaluates this).
- Radical perfect-square extraction: `sqrt(8)` stays `8**(1/2)` instead
  of `2*sqrt(2)`; `(-4)**(1/2)` stays instead of `2*I`. Surfaces as
  unsimplified roots from `solve` (e.g. `1/2*8**(1/2)` for `x**2-2`
  rather than `sqrt(2)`).
- **Status:** open, tracked as parity gaps rather than correctness bugs
  (values are correct, only the canonical form differs).
