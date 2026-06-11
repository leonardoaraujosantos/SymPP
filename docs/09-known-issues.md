# Known issues & SymPy-parity bug log

This file tracks correctness gaps found by cross-checking SymPP against
the SymPy 1.14 oracle. Each entry records the failing input, the
SymPy-expected result, current status, and the regression test (once
fixed). Each entry is filed as a GitHub issue; the local log is the source of
truth and links the issue number.

> Networking note: this environment's DNS resolves `api.github.com` and
> `github.com` to a dead Azure block (`4.228.31.x`), so `gh` and plain
> `curl` to the API time out. Git push over HTTPS still works. To reach
> the REST API, force a real GitHub IP, e.g.
> `curl --resolve api.github.com:443:140.82.121.6 ...`.

## Fixed

### SUM-1 — geometric summation dropped non-trivial exponents ([#1], PR [#4])
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

### SIMP-1 — `simplify` could return a *more complex* expression ([#3])
- **Input:** `simplify((x**2 - 1)/(x - 1))`
- **Was:** `((x - 1)*x**2 - (x - 1))*(x - 1)**(-2)` (worse than input).
- **Expected (SymPy):** `x + 1`.
- **Cause:** the `simplify` pipeline never reduced `n/d` by polynomial GCD.
- **Fix:** added a rational-cancellation step that calls `cancel()` and
  adopts the result only when strictly fewer nodes. Restricted to
  *univariate* rational functions with a symbol-dependent denominator —
  `cancel()` hangs on transcendental and on multivariate-symbolic input
  (see CANCEL-1), and the strictly-simpler guard means `simplify` never
  returns something larger than the pipeline already produced.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[simplify][regression]` cases ((x²-1)/(x-1), (x²+2x+1)/(x+1),
  multivariate-terminates, sin²+cos² no-hang guard).
- **Scope:** multivariate cancellation (e.g. (x²-y²)/(x-y) → x+y) is
  deliberately deferred until CANCEL-1 is fixed.

### INT-1 — `integrate(exp(x)*sin(x))` segfaulted ([#7])
- **Input:** `integrate(exp(x)*sin(x), x)` (and `exp·cos`, `exp(2x)·sin(3x)`, …).
- **Was:** **SIGSEGV** — integration by parts recursed
  `exp·sin → exp·cos → exp·sin → …` without bound and overflowed the stack.
- **Expected (SymPy):** `exp(x)*sin(x)/2 - exp(x)*cos(x)/2`.
- **Fix:** (1) a dedicated closed-form rule for `c·e^(a x+·)·sin/cos(g x+·)`
  — `∫E·S = E(aS−gC)/(a²+g²)`, `∫E·C = E(aC+gS)/(a²+g²)`; (2) a
  recursion-depth backstop in `integrate()` (limit 64) so any *other*
  cyclic integrand degrades to the unevaluated `Integral(...)` marker
  instead of crashing.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[byparts][regression]` (exp·sin, exp·cos, exp(2x)·sin(3x), and x²·exp
  to confirm polynomial by-parts still works), each verified by
  differentiating the result back to the integrand.

### INT-2 — `integrate(sinh/cosh)` returned the unevaluated marker ([#9])
- **Input:** `integrate(sinh(x))`, `integrate(cosh(3*x))`, …
- **Was:** `Integral(sinh(x), x)` (table fell through).
- **Expected (SymPy):** `cosh(x)`, `sinh(3*x)/3`, …
- **Fix:** added `Sinh`/`Cosh` cases to the affine-argument table:
  `∫sinh(ax+b) = cosh(ax+b)/a`, `∫cosh(ax+b) = sinh(ax+b)/a`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[hyperbolic][regression]`.
- **Scope:** polynomial × hyperbolic (e.g. `∫x·cosh(x)`) still deferred —
  the by-parts target set is `{exp,sin,cos}`.

### INT-3 — `integrate(tan / 1/cos² / 1/sin²)` returned the unevaluated marker
- **Input:** `integrate(tan(x))`, `integrate(1/cos(x)**2)`,
  `integrate(1/sin(x)**2)`, and their affine-argument variants.
- **Was:** `Integral(tan(x), x)`, `Integral(cos(x)**(-2), x)`, … (table fell
  through — only `sin`/`cos`/`exp`/`sinh`/`cosh` of an affine argument were
  tabulated, and the `Pow` branch only handled affine bases).
- **Expected (SymPy):** `-log(cos(x))`, `tan(x)`, `-cot(x)`, …
- **Fix:** added a `Tan` case to the affine-argument function table
  (`∫tan(ax+b) = -log(cos(ax+b))/a`) and a reciprocal-square trig case to the
  `Pow` branch (`∫1/cos²(ax+b) = sin/(a·cos)`, `∫1/sin²(ax+b) = -cos/(a·sin)`).
  SymPP emits the `sin/cos` forms, equivalent to SymPy's `tan`/`-cot`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[trig][regression]`.
- **Scope:** `sec`/`csc`/`cot` are not distinct function types in SymPP, so
  results are spelled with `sin`/`cos`. Inverse-trig antiderivatives
  (`∫1/(1+x²) = atan`, `∫1/√(1-x²) = asin`) remain deferred.

### INT-4 — `integrate(xⁿ·log(x))` returned the unevaluated marker
- **Input:** `integrate(x*log(x))`, `integrate(x**2*log(x))`,
  `integrate((x+1)*log(x))`, `integrate(x*log(2*x+1))`.
- **Was:** `Integral(x*log(x), x)` — integration by parts only ever used
  `sin`/`cos`/`exp` of an affine argument as the `dv` factor, never `log`, so
  a polynomial × `log` product fell through. (Standalone `∫log(ax+b)` already
  worked via its own branch.)
- **Expected (SymPy):** `x²·log(x)/2 − x²/4`, `x³·log(x)/3 − x³/9`, …
- **Fix:** added a by-parts branch with `u = log(ax+b)`, `dv = rest dx`:
  `∫rest·log(ax+b) = log(ax+b)·∫rest − ∫(∫rest)·a/(ax+b)`. The trailing
  integral is rational (∫rest is polynomial, `du = a/(ax+b)`), so
  `try_rational` closes it; the marker/depth guards bail on anything that does
  not reduce. The result is `expand`ed for a distributed polynomial form.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[byparts][log][regression]`.
- **Scope:** `log` powers (`∫log(x)²`, `∫x·log(x)²`) still defer — they are
  `Pow(log, n)`, not a single `Log` factor, and need recursive by-parts.
  For an affine log argument the primitive matches SymPy only up to an
  additive constant (`log(x+1/2)` vs `log(2x+1)`); the derivative is exact.

### SOLVE-1 — `solve()` returned empty for transcendental equations ([#11])
- **Input:** `solve(log(x) - 1, x)`, `solve(exp(x) - 2, x)`, …
- **Was:** `[]` — the vector `solve` was polynomial-only (`Poly.roots()`),
  even though `solveset(log(x)-1)` already returned `{E}`.
- **Expected (SymPy):** `[E]`, `[log(2)]`, `[asinh(1)]`.
- **Fix:** when the polynomial path is empty and the expression contains a
  function of `var`, route through `solveset` and surface a `FiniteSet`
  result as the root vector. The polynomial-only logic was split into a
  `solve_poly` helper, and `solveset`'s internal fallback now calls *that*
  (not the public `solve`) to avoid `solve ↔ solveset` infinite recursion.
- **Regression test:** `tests/solvers/solve_test.cpp`
  — `[transcendental][regression]` (log, exp, sinh) plus a polynomial guard.
- **Scope:** periodic/infinite solution sets (e.g. `sin(x)=0`) remain the
  domain of `solveset`; `solve` yields no finite vector for those.

### DIFF-1 — `diff(Abs(x))` returned 0 instead of `sign(x)` ([#13])
- **Input:** `diff(abs(x), x)`
- **Was:** `0` — `Abs` had no `diff_arg` override, so it fell through to
  `Function::diff_arg`'s default of `0`.
- **Expected (SymPy):** `sign(x)`.
- **Fix:** `Abs::diff_arg(i) = sign(arg)`; `diff()`'s chain rule supplies
  the `arg'` factor (so `diff(abs(2x+1)) = 2*sign(2x+1)`,
  `diff(x*abs(x)) = x*sign(x) + Abs(x)`).
- **Regression test:** `tests/calculus/diff_test.cpp`
  — `[diff][abs][regression]`.
- **Minor follow-on:** `diff(abs(x**2))` gives `2*x*sign(x**2)` (correct
  but unsimplified — `sign` doesn't yet auto-reduce `sign(x**2) → 1` for a
  manifestly-nonnegative argument).

## Open

### CANCEL-1 — `cancel()`/`Poly` GCD hangs on symbolic coefficients ([#5])
- **Input:** `cancel((b - a + 1)*(a + b)/2, a)` (and any multivariate
  polynomial whose coefficients in the cancel variable are symbolic).
- **Is:** infinite loop — never returns.
- **Expected (SymPy):** returns the already-reduced expression immediately.
- **Cause:** the polynomial-remainder-sequence GCD does not terminate when
  the leading coefficient is itself symbolic and the other operand is a
  degree-0 constant.
- **Impact:** reached via `simplify`/`summation`; `simplify` now guards
  against it (univariate-only cancel — see SIMP-1). Also makes
  `factor(x**2 - y**2, x)` (multivariate square-free factorization, which
  uses the same GCD) hang. The root `cancel()`/`Poly::gcd` loop is open.
- **Status:** open. Regression test to be added with the fix.

### LIM-1 — limit of the classic `e` definition returns garbage ([#2])
- **Input:** `limit((1 + 1/x)**x, x, oo)`
- **Is:** `(oo**(-1) + 1)**oo` (unevaluated, nonsensical).
- **Expected (SymPy):** `E`.
- **Notes:** the limit engine substitutes the target naively instead of
  detecting the `1^∞` indeterminate form (needs `exp(limit(x*log(1+1/x)))`
  rewrite). Partly blocked by the missing `Infinity` singleton — `oo` is
  currently parsed as a plain symbol.
- **Status:** open. Regression test to be added with the fix.

### Missing auto-simplifications (lower priority)
- `exp(log(x))` does not reduce to `x` (SymPy auto-evaluates this).
- Radical perfect-square extraction: `sqrt(8)` stays `8**(1/2)` instead
  of `2*sqrt(2)`; `(-4)**(1/2)` stays instead of `2*I`. Surfaces as
  unsimplified roots from `solve` (e.g. `1/2*8**(1/2)` for `x**2-2`
  rather than `sqrt(2)`).
- **Status:** open, tracked as parity gaps rather than correctness bugs
  (values are correct, only the canonical form differs).

[#1]: https://github.com/leonardoaraujosantos/SymPP/issues/1
[#2]: https://github.com/leonardoaraujosantos/SymPP/issues/2
[#3]: https://github.com/leonardoaraujosantos/SymPP/issues/3
[#4]: https://github.com/leonardoaraujosantos/SymPP/pull/4
[#5]: https://github.com/leonardoaraujosantos/SymPP/issues/5
[#7]: https://github.com/leonardoaraujosantos/SymPP/issues/7
[#9]: https://github.com/leonardoaraujosantos/SymPP/issues/9
[#11]: https://github.com/leonardoaraujosantos/SymPP/issues/11
[#13]: https://github.com/leonardoaraujosantos/SymPP/issues/13
