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
- **Scope:** `sec`/`csc`/`cot` now exist as distinct function types (TRIG-RECIP)
  but their antiderivatives are a separate item; inverse-trig antiderivatives
  (`∫1/(1+x²) = atan`, `∫1/√(1-x²) = asin`) are handled by INT-5/INT-6.

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

### INT-5 — `integrate(1/(a·x²+b·x+c))` (irreducible) returned the marker
- **Input:** `integrate(1/(1+x**2))`, `integrate(1/(x**2+4))`,
  `integrate(1/(4*x**2+9))`, `integrate(1/(x**2+2*x+5))`.
- **Was:** `Integral((x**2 + 1)**(-1), x)` — `try_rational` only decomposes
  denominators with *real* roots (via `apart`); an irreducible quadratic
  (negative discriminant, complex roots) had no closed-form path.
- **Expected (SymPy):** `atan(x)`, `atan(x/2)/2`, `atan(2*x/3)/6`,
  `atan((x+1)/2)/2`.
- **Fix:** added `try_arctan_quadratic`, dispatched right after
  `try_rational`. For `1/(a·x²+b·x+c)` with `D = 4ac − b² > 0` it returns
  `2·atan((2ax+b)/√D)/√D`. Requires rational coefficients; `D ≤ 0` (real
  roots) falls through to `try_rational`, so `1/(x²−1)` still yields logs.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[arctan][regression]` (incl. a reducible `1/(x²−1)` guard).
- **Scope:** the `√(quadratic)` reciprocals (`∫1/√(1−x²) = asin`,
  `∫1/√(x²+1) = asinh`) are still deferred — a separate branch on the
  `−1/2` exponent. Symbolic coefficients (`1/(k²+x²)`) are out of scope.

### INT-6 — `integrate(1/sqrt(a·x²+c))` returned the unevaluated marker
- **Input:** `integrate(1/sqrt(1-x**2))`, `integrate(1/sqrt(9-4*x**2))`,
  `integrate(1/sqrt(x**2+1))`, `integrate(1/sqrt(4*x**2+9))`.
- **Was:** `Integral((-x**2 + 1)**(-1/2), x)` — no path handled the `−1/2`
  exponent over a quadratic radicand.
- **Expected (SymPy):** `asin(x)`, `asin(2*x/3)/2`, `asinh(x)`,
  `asinh(2*x/3)/2`.
- **Fix:** added `try_sqrt_quadratic`, dispatched after
  `try_arctan_quadratic`. For a pure quadratic radicand `a·x²+c` (no linear
  term) with `c > 0`: `a > 0 → asinh(x·√(a/c))/√a`,
  `a < 0 → asin(x·√(−a/c))/√(−a)`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[invtrig][regression]`.
- **Scope:** a linear term under the root is out of scope; the `c < 0`
  (acosh/log) family is handled by INT-7. (The perfect-square
  coefficients now print reduced thanks to SQRT-1.)

### SQRT-1 — `sqrt` of a perfect-square *rational* was not reduced
- **Input:** `sqrt(rational(1,4))`, `sqrt(rational(9,4))`.
- **Was:** `(1/4)**(1/2)`, `(9/4)**(1/2)` — left symbolic, even though
  `sqrt(integer(4))` reduced to `2`. The integer and rational paths were
  inconsistent.
- **Expected (SymPy):** `1/2`, `3/2`.
- **Fix:** generalised `try_integer_perfect_root` → `try_perfect_root` in
  `src/core/pow.cpp` to accept a non-negative Rational base, rooting numerator
  and denominator independently (`√(9/4) = √9/√4 = 3/2`). This also cleans up
  the previously unreduced coefficients in the arctan/asin/asinh integration
  results (INT-5, INT-6) — e.g. `asin(2*x/3)/2` now prints directly.
- **Regression test:** `tests/functions/miscellaneous_test.cpp`
  — `[sqrt][regression]`.
- **Scope:** perfect-square *factor extraction* (`√8 → 2√2`, `√(1/2) → √2/2`)
  is a further SymPy behaviour still not implemented; non-perfect-square
  rationals stay a symbolic `Pow`.

### INT-7 — `integrate(1/sqrt(a·x²+c))` with `c < 0` returned the marker
- **Input:** `integrate(1/sqrt(x**2-1))`, `integrate(1/sqrt(4*x**2-9))`.
- **Was:** `Integral((x**2 - 1)**(-1/2), x)` — INT-6 only covered `c > 0`
  (asin / asinh); the `a > 0, c < 0` case was explicitly deferred.
- **Expected (SymPy):** `log(x + sqrt(x**2 - 1))`,
  `log(2*x + sqrt(4*x**2 - 9))/2`.
- **Fix:** extended `try_sqrt_quadratic` with the `a > 0, c < 0` branch:
  `∫1/√(a·x²+c) = log(√a·x + √(a·x²+c))/√a` (acosh-equivalent, the form SymPy
  prints — and `√a` now reduces for perfect-square `a`, cf. SQRT-1/SQRT-2).
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[invtrig][regression]`.
- **Scope:** a linear term under the root, and the `a < 0, c < 0` case
  (radicand never positive), still fall through.

### INT-8 — `integrate(tan(x)**2)` returned the unevaluated marker
- **Input:** `integrate(tan(x)**2)`, `integrate(tan(2*x+1)**2)`.
- **Was:** `Integral(tan(x)**2, x)` — only `sin²`/`cos²` had a trig-reduction
  rewrite; `tan²` fell through.
- **Expected (SymPy):** `-x + sin(x)/cos(x)` (= `tan(x) - x`).
- **Fix:** added a `tan²(u) → 1/cos²(u) − 1` (Pythagorean) rewrite to
  `try_trig_reduction`, guarded to an affine `u` so the recursion lands on the
  tabulated `∫1/cos²(u)` (INT-3). Result: `tan(u)/a − u`, spelled with
  `sin/cos` and confirmed equivalent to SymPy via the oracle.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[trig][regression]`.
- **Scope:** non-affine arguments, and higher powers (`tan⁴`, `sec⁴`), are not
  handled.

### TRIG-1 — `sin`/`cos`/`tan` not evaluated at rational multiples of π
- **Input:** `sin(pi/6)`, `cos(pi/3)`, `tan(pi/4)`, `sin(2*pi/3)`,
  `cos(5*pi/6)`, …
- **Was:** `sin(1/6*pi)`, … — only `0`, `π/2`, `π` were special-cased (the
  `π/2` case via a brittle two-factor `Mul` match).
- **Expected (SymPy):** `1/2`, `1/2`, `1`, `sqrt(3)/2`, `-sqrt(3)/2`, …
- **Fix:** added a `pi_coefficient` helper (recognises `r·π` for rational `r`)
  plus exact-value tables with full period/quadrant reduction:
  `cos_pi`/`sin_pi` (denominators 1,2,3,4,6) and a dedicated `tan_pi` (clean
  `√3/3`, `√3`, `1`). Poles (`tan(π/2)`) and out-of-table denominators
  (`sin(π/12)`) are left unevaluated. The old `π/2` `Mul`-match special cases
  were removed — the helper subsumes them.
- **Regression test:** `tests/functions/trigonometric_test.cpp`
  — `[trig][regression]` (rational + radical values, quadrant images, pole and
  out-of-table guards).
- **Scope:** denominators outside {1,2,3,4,6,12} (e.g. `π/5`, `π/8`) stay
  symbolic. Denominator 12 ships in TRIG-PI12; `π/8` is a genuinely nested
  radical, deferred.

### TRIG-PI12 — exact values at multiples of π/12 (15°)
- **Input:** `cos(π/12)`, `sin(π/12)`, `tan(π/12)`, `cos(5π/12)`, `tan(5π/12)`,
  `cos(7π/12)`, …
- **Was:** denominator 12 fell outside the {1,2,3,4,6} table → left symbolic
  (and an old test wrongly called π/12 a "nested radical").
- **Expected (SymPy):** `cos(π/12) = (√6+√2)/4`, `sin(π/12) = (√6−√2)/4`,
  `tan(π/12) = 2−√3`, `tan(5π/12) = 2+√3`, with the usual quadrant signs.
- **Fix (`src/functions/trigonometric.cpp`):** added the `den = 12` reference
  values (`r = 1/12, 5/12`) to `base_cos_pi` and `base_tan_pi`; all 24 multiples
  reduce in through the existing `cos_pi`/`tan_pi` period + reflection folds, and
  `sin_pi` gets them via the `cos((1/2−r)π)` co-function identity.
- **Regression test:** `tests/functions/trigonometric_test.cpp`
  — `[trig][regression]` (π/12 family + a `π/8` nested-radical no-op guard).
- **Scope:** π/12 (15°) only. π/8 (22.5°) and π/5 (36°) need nested-radical /
  golden-ratio forms and stay deferred.

### TRIG-2 — `asin`/`acos`/`atan` not evaluated at special arguments
- **Input:** `asin(1/2)`, `acos(1/2)`, `atan(sqrt(3))`, `asin(sqrt(2)/2)`,
  `acos(-1/2)`, `atan(sqrt(3)/3)`.
- **Was:** `asin(1/2)`, … — only the trivial `0`, `±1` arguments folded.
- **Expected (SymPy):** `pi/6`, `pi/3`, `pi/3`, `pi/4`, `2*pi/3`, `pi/6`.
- **Fix:** `asin_special` / `atan_special` reverse-lookup tables (matching the
  same radical constants the forward TRIG-1 table emits, so structural
  equality fires), with oddness routed through the factory so negatives fold.
  `acos(x) = π/2 − asin(x)`, adopted only when `asin` produced an exact angle
  (otherwise `acos(x)` stays unevaluated, as SymPy does).
- **Regression test:** `tests/functions/inverse_trig_test.cpp`
  — `[asin]/[acos]/[atan][regression]` (incl. negative args and a
  non-special `acos(1/3)` guard).
- **Scope:** mirrors TRIG-1 — only the {1,2,3,4,6}-denominator angles; other
  arguments (e.g. `asin(1/3)`) stay symbolic.

### EXP-1 — `exp` not evaluated at imaginary multiples of π (Euler identity)
- **Input:** `exp(I*pi)`, `exp(2*I*pi)`, `exp(I*pi/2)`, `exp(-I*pi/2)`.
- **Was:** `exp(pi*I)`, … — left unevaluated.
- **Expected (SymPy):** `-1`, `1`, `I`, `-I`.
- **Fix:** added an `imaginary_pi_coeff` helper (detects `r·I·π`) and the Euler
  rule `exp(r·I·π) = i^(2r)` when `2r` is an integer (`pow(I, n)` already
  cycles through {1, I, −1, −I}). This matches SymPy, which folds only the
  `q ∈ {1,2}` coefficients and keeps `exp(I·π/3)`, `exp(I·π/4)` symbolic.
- **Regression test:** `tests/functions/exponential_test.cpp`
  — `[exp][regression]` (±1/±I values, periodicity, and a `π/3` symbolic guard).
- **Scope:** only half-integer coefficients; `exp(I·x)` for symbolic `x` is not
  expanded to `cos + I·sin` (that is `expand_complex`/`rewrite`, not auto-eval).

### INT-9 — `integrate(1/(a·x²+b·x+c))` (repeated root) returned the marker
- **Input:** `integrate(1/(x**2+2*x+1))`, `integrate(1/(4*x**2+4*x+1))`.
- **Was:** `Integral((x**2 + 2*x + 1)**(-1), x)` — the denominator is a perfect
  square `a·(x−r)²` (discriminant 0); `apart` did not decompose the repeated
  root and the arctan branch (INT-5) needs `D > 0`.
- **Expected (SymPy):** `-1/(x + 1)`, `-1/(2*(2*x + 1))`.
- **Fix:** added the `D = 0` case to `try_arctan_quadratic`:
  `∫1/(a·(x−r)²) = −2/(2ax+b)`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[arctan][regression]`.
- **Scope:** distinct-real-root (`D < 0`) denominators still go through
  `try_rational` (logs); irreducible (`D > 0`) through the arctan branch.

### TRIG-3 — `sin`/`cos`/`tan` did not reduce arguments by π multiples
- **Input:** `sin(x+2*pi)`, `sin(x+pi)`, `cos(x+pi)`, `tan(x+pi)`,
  `sin(x+y+pi)`.
- **Was:** `sin(x + 2*pi)`, … — left unevaluated.
- **Expected (SymPy):** `sin(x)`, `-sin(x)`, `-cos(x)`, `tan(x)`, `-sin(x+y)`.
- **Fix:** a `split_pi_term` helper sums the rational π-coefficient `C` across
  an additive argument; when `C` is a nonzero integer `k`,
  `sin(rest+kπ)=(−1)^k sin(rest)`, `cos` likewise, and `tan(rest+kπ)=tan(rest)`
  (period π).
- **Regression test:** `tests/functions/trigonometric_test.cpp`
  — `[trig][regression]`.
- **Scope:** integer multiples of π only. Half-integer shifts (the co-function
  `sin(x+π/2)=cos(x)`) stay symbolic — a separate follow-up.

### TRIG-4 — `sin`/`cos` did not apply the half-integer-π co-function shift
- **Input:** `sin(x+pi/2)`, `cos(x+pi/2)`, `sin(x-pi/2)`, `cos(x+3*pi/2)`.
- **Was:** `sin(x + pi/2)`, … — TRIG-3 reduced only integer π multiples;
  half-integer shifts stayed symbolic.
- **Expected (SymPy):** `cos(x)`, `-sin(x)`, `-cos(x)`, `sin(x)`.
- **Fix:** extended the `split_pi_term` reduction with the `C = m/2` (m odd)
  case: `sin(rest+(m/2)π) = ±cos(rest)`, `cos(rest+(m/2)π) = ∓sin(rest)`,
  the sign from `m mod 4`.
- **Regression test:** `tests/functions/trigonometric_test.cpp`
  — `[trig][regression]`.
- **Scope:** `tan(x+π/2) = −cot(x)` is left symbolic — SymPP has no `cot`
  function type.

### INT-10 — heurisch missed a u-sub when the inner function was a Mul factor
- **Input:** `integrate(1/(x*log(x)))`, `integrate(1/(x*log(x)**2))`.
- **Was:** `Integral((x*log(x))**(-1), x)` — heurisch only collected function
  *arguments* and `Pow` *bases* as substitution candidates, so `log(x)` buried
  as a factor of the `Pow` base `x·log(x)` was never tried.
- **Expected (SymPy):** `log(log(x))`, `-1/log(x)`.
- **Fix:** the candidate `walk` now also adds the function application itself
  (e.g. `log(x)`), so `u = log(x)` is considered.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[heurisch][regression]`.
- **Scope:** still single-substitution heurisch; integrands needing erf/erfi
  (`∫exp(x²)`) remain unevaluated (no `erfi` function type).

### SUM-2 — arithmetic-geometric `Σ k·r^k` returned the summand unchanged
- **Input:** `summation(k*2**k, k, 0, n)`, `summation(k*3**k, k, 0, n)`.
- **Was:** `k*2**k` — a `Mul` of two var-dependent factors isn't split by the
  constant-extraction path, and the geometric handler only matched a pure
  `base^(linear·k)`.
- **Expected (SymPy):** `2*2**n*n - 2*2**n + 2`, ….
- **Fix:** added an arithmetic-geometric case `Σ k·r^k` for a numeric ratio
  `r = base^c ≠ 1`, using the closed form
  `Σ_{k=0}^{N} k·r^k = r(1 − (N+1)r^N + N·r^{N+1})/(1−r)²` with telescoping for
  general bounds; the `base^d` prefactor factors out.
- **Regression test:** `tests/calculus/series_limit_test.cpp`
  — `[summation][regression]`.
- **Scope:** numeric ratio only (a symbolic `r` would need a Piecewise on
  `r = 1`, as SymPy emits); higher-degree `P(k)·r^k` still defers.

### SUM-3 — an unrecognised sum collapsed to its bare summand
- **Input:** `summation(1/k**2, k, 1, oo)`, `summation(1/k, k, 1, oo)`,
  `summation(1/factorial(k), k, 0, oo)`.
- **Was:** `k**(-2)`, … — the fallback `return expr` handed back the summand,
  so `Σ_{k=1}^∞ 1/k²` came out as `1/k²` (the summation silently dropped).
- **Expected (SymPy):** a closed form where one exists, otherwise an
  unevaluated `Sum`. SymPP does not yet close `ζ(2)` / exponential series, so
  it should at least preserve the sum.
- **Fix:** the no-closed-form fallback now returns an unevaluated
  `Sum(expr, var, lo, hi)` marker (an `UndefinedFunction`, mirroring the
  `Integral(_, _)` marker), never the bare summand. Also added the
  single-term range rule `Σ_{k=a}^{a} f(k) = f(a)`. Infinite *geometric* series
  already close (the `ratio^oo → 0` fold from the Infinity work):
  `Σ (1/2)^k = 2`, and divergent `Σ k = oo`.
- **Regression test:** `tests/calculus/series_limit_test.cpp`
  — `[summation][regression]`.
- **Scope:** convergent non-geometric series (`Σ 1/k² = π²/6`, `Σ 1/k! = E`)
  stay as `Sum` markers — closing them needs `zeta` / series recognition,
  deferred. The fix guarantees correctness (no dropped sum), not closure.

### FUNC-1 — `f(f⁻¹(x))` not simplified to `x`
- **Input:** `sin(asin(x))`, `cos(acos(x))`, `tan(atan(x))`, `sinh(asinh(x))`,
  `cosh(acosh(x))`, `tanh(atanh(x))`.
- **Was:** `sin(asin(x))`, … — left unevaluated.
- **Expected (SymPy):** all `x`.
- **Fix:** an `arg_of` helper in the trig and hyperbolic factories returns the
  inner argument when the forward function wraps its own inverse, collapsing
  `f(f⁻¹(x)) → x`.
- **Regression test:** `tests/functions/inverse_trig_test.cpp` and
  `tests/functions/hyperbolic_test.cpp` — `[regression]`.
- **Scope:** only the `f(f⁻¹)` direction. The reverse `f⁻¹(f(x))`
  (e.g. `asin(sin(x))`) stays unevaluated — it is `x` only on a restricted
  range, matching SymPy.

### FUNC-2 — cross-function inverse compositions not simplified
- **Input:** `cos(asin(x))`, `sin(acos(x))`, `tan(asin(x))`, `cos(atan(x))`,
  `sin(atan(x))`, `tan(acos(x))`.
- **Was:** `cos(asin(x))`, … — left unevaluated.
- **Expected (SymPy):** `√(1−x²)`, `√(1−x²)`, `x/√(1−x²)`, `1/√(1+x²)`,
  `x/√(1+x²)`, `√(1−x²)/x`.
- **Fix:** the `sin`/`cos`/`tan` factories now recognise a different inverse-trig
  argument (via `arg_of`) and emit the corresponding algebraic form. Extends
  FUNC-1 (the same-function `f(f⁻¹)` collapse).
- **Regression test:** `tests/functions/inverse_trig_test.cpp` — `[regression]`.
- **Scope:** the trig × inverse-trig table; hyperbolic cross-compositions
  (`cosh(asinh(x)) = √(x²+1)`, …) are a follow-up.

### FUNC-3 — hyperbolic cross-function inverse compositions not simplified
- **Input:** `cosh(asinh(x))`, `sinh(acosh(x))`, `tanh(asinh(x))`,
  `cosh(atanh(x))`, `sinh(atanh(x))`, `tanh(acosh(x))`.
- **Was:** `cosh(asinh(x))`, … — left unevaluated.
- **Expected (SymPy):** `√(x²+1)`, `√(x−1)·√(x+1)`, `x/√(x²+1)`, `1/√(1−x²)`,
  `x/√(1−x²)`, `√(x−1)·√(x+1)/x`.
- **Fix:** the hyperbolic analogue of FUNC-2 — `sinh`/`cosh`/`tanh` recognise a
  different inverse-hyperbolic argument (via `arg_of`) and emit the algebraic
  form. The `acosh` cases use `√(x−1)·√(x+1)`, the form SymPy prints.
- **Regression test:** `tests/functions/hyperbolic_test.cpp` — `[regression]`.

### DIFF-2 — `diff(erf/erfc/Heaviside)` returned 0
- **Input:** `diff(erf(x), x)`, `diff(erfc(x), x)`, `diff(Heaviside(x), x)`.
- **Was:** `0` — these classes had no `diff_arg` override, so they fell through
  to `Function::diff_arg`'s default of `0` (the same root cause as DIFF-1/Abs).
- **Expected (SymPy):** `2*exp(-x**2)/sqrt(pi)`, `-2*exp(-x**2)/sqrt(pi)`,
  `DiracDelta(x)`.
- **Fix:** added `diff_arg` to `Erf`, `Erfc`, `HeavisideFn` —
  `erf' = 2·exp(−x²)/√π`, `erfc' = −that`, `Heaviside' = DiracDelta(x)`. The
  chain rule supplies the `arg'` factor (so `diff(erf(2x)) = 4·exp(−4x²)/√π`).
- **Regression test:** `tests/functions/special_test.cpp` — `[diff][regression]`.
- **Scope:** `gamma`/`loggamma` derivatives need `digamma`/`polygamma`
  (not yet a function type) and stay at 0; `sign`/`floor`/`re`/`im`/`conjugate`
  match SymPy in keeping an unevaluated/zero derivative.

### INT-11 — `integrate(exp(-a·x²))` (Gaussian) returned the marker
- **Input:** `integrate(exp(-x**2))`, `integrate(exp(-x**2/2))`,
  `integrate(2*exp(-x**2)/sqrt(pi))`.
- **Was:** `Integral(exp(-x**2), x)` — no error-function path.
- **Expected (SymPy):** `sqrt(pi)*erf(x)/2`, `sqrt(2*pi)*erf(sqrt(2)*x/2)/2`,
  `erf(x)`.
- **Fix:** added `try_gaussian`: for `exp(c·x²)` with a concrete negative
  rational `c`, `∫ = √π·erf(√a·x)/(2√a)`, `a = −c`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[erf][regression]`.
- **Scope:** pure `c·x²` exponent (no linear/constant term — completing the
  square is out of scope); positive `c` would need `erfi` (no such function
  type). Pairs with DIFF-2 (the `erf` derivative).

### INT-12 — `integrate(tanh / 1/cosh² / 1/sinh²)` returned the marker
- **Input:** `integrate(tanh(x))`, `integrate(1/cosh(x)**2)`,
  `integrate(1/sinh(x)**2)`, and affine-argument variants.
- **Was:** `Integral(tanh(x), x)`, … — the hyperbolic counterparts of INT-3
  were missing (only `sinh`/`cosh` of an affine argument were tabulated).
- **Expected (SymPy):** `log(cosh(x))`, `tanh(x)`, `-cosh(x)/sinh(x)` (= −coth).
- **Fix:** added a `Tanh` case to the affine-argument function table
  (`∫tanh(ax+b) = log(cosh(ax+b))/a`) and `Cosh`/`Sinh` reciprocal-square cases
  to the `Pow` branch (`∫1/cosh²(ax+b) = tanh(ax+b)/a`,
  `∫1/sinh²(ax+b) = -cosh/(a·sinh)`).
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[hyperbolic][regression]`.
- **Scope:** `coth`/`sech`/`csch` now exist as distinct function types
  (HYP-RECIP); their antiderivatives are a separate item.

### INT-13 — `integrate(poly·cosh / poly·sinh)` returned the marker
- **Input:** `integrate(x*cosh(x))`, `integrate(x*sinh(x))`,
  `integrate(x**2*cosh(x))`, `integrate(x*cosh(2*x+1))`.
- **Was:** `Integral(x*cosh(x), x)`, … — integration by parts only recognised
  `{exp, sin, cos}` of an affine argument as the `dv` factor, so a polynomial
  times `sinh`/`cosh` fell through to the unevaluated marker.
- **Expected (SymPy):** `x*sinh(x) - cosh(x)`, `x*cosh(x) - sinh(x)`, etc.
- **Fix:** added `FunctionId::Sinh` and `FunctionId::Cosh` to the by-parts
  target-set condition in `try_integration_by_parts`
  (`src/integrals/integrate.cpp`). The polynomial `u` is differentiated down
  each step, so the recursion terminates — `sinh`/`cosh` don't cycle the way
  `exp·sin/cos` does.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[parts][hyperbolic][regression]`.

### INT-14 — `integrate(log(x)**n)` / `integrate(poly·log(x)**n)` returned the marker
- **Input:** `integrate(log(x)**2)`, `integrate(log(x)**3)`,
  `integrate(x*log(x)**2)`, `integrate(log(2*x)**2)`.
- **Was:** `Integral(log(x)**2, x)`, … — integration by parts only recognised a
  single power-1 `log(affine)` factor (INT-4), so any `log` raised to an
  integer power fell through to the unevaluated marker.
- **Expected (SymPy):** `x*log(x)**2 - 2*x*log(x) + 2*x`, etc.
- **Fix:** added `is_log_or_log_power` (accepts `log(affine)` or a positive
  integer power of one) in `src/integrals/integrate.cpp`, a standalone
  `log(affine)**n` by-parts branch (`u = log**n, dv = dx, v = x`), and relaxed
  the existing polynomial×log branch to use the same predicate. By parts
  reduces the exponent each step (`(log**n)' = n·log**(n-1)·a/(ax+b)`), so it
  recurses down to the `∫log` case; the marker guard bails on anything that
  does not reduce, so it never loops or emits a wrong closed form.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[parts][log][regression]`.
- **Scope:** affine arguments with a non-zero constant term (`log(2x+1)**2`)
  may stay symbolic — the remaining `∫x·log**(n-1)/(ax+b)` does not always
  close, in which case the marker guard leaves it unevaluated.

### INT-15 — `integrate(exp(x)/x)` hung (non-elementary by-parts)
- **Input:** `integrate(exp(x)/x, x)` (the non-elementary `Ei(x)`).
- **Was:** infinite loop. Integration by parts took `u = x^(-1)`,
  `dv = exp(x) dx`, producing `∫exp(x)/x = exp(x)/x + ∫exp(x)/x²`, then
  `∫exp(x)/x²`, `∫exp(x)/x³`, … — each step *raises* the negative power, so the
  recursion never terminated (the depth guard only bounds a single chain; the
  branching across `try_*` made it effectively hang).
- **Expected (SymPy):** `Ei(x)`. SymPP has no `Ei`, so the correct fallback is
  the unevaluated `Integral` marker — and crucially it must *terminate*.
- **Fix:** the poly×{exp,sin,cos,sinh,cosh} by-parts branch now requires the
  `u` factor to be a polynomial in `var` (`is_polynomial_in`), so its
  derivative chain reaches zero in finitely many steps. A non-polynomial `u`
  like `x^(-1)` (derivatives `x^(-2)`, `x^(-3)`, … grow) is rejected and the
  integral returns the marker (`src/integrals/integrate.cpp`).
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][parts][regression]`.
- **Scope:** poly·exp/trig/hyperbolic by parts is unchanged (`u` is a genuine
  polynomial). Recognising `Ei`/`Si`/`Ci` special-function antiderivatives is a
  separate feature.

### INT-16 — `integrate((linear)/(irreducible quadratic))` returned the marker
- **Input:** `∫(x+1)/(x²+2x+5)`, `∫(2x+3)/(x²+1)`, `∫(3x+5)/(x²+4)`,
  `∫x/(x²+2x+5)`, `∫x/(x²+1)`.
- **Was:** the marker, or — for `∫x/(x²+1)` via the logarithmic-derivative
  heuristic — the spurious `1/2*log(2*(x²+1))` (extra factor inside the log).
  `try_arctan_quadratic` only handled a *constant* numerator (`1/quad`).
- **Expected (SymPy):** `log(x²+2x+5)/2`, `log(x²+1)+3*atan(x)`,
  `3*log(x²+4)/2 + 5*atan(x/2)/2`, `log(x²+2x+5)/2 − atan(x/2+1/2)/2`,
  `log(x²+1)/2`.
- **Fix:** added `try_linear_over_quadratic` (`src/integrals/integrate.cpp`):
  for `(p·x+q)/(a·x²+b·x+c)` over an irreducible quadratic (`4ac−b² > 0`), split
  `p·x+q = (p/2a)(2a·x+b) + (q−pb/2a)` ⇒
  `(p/2a)·log(a·x²+b·x+c) + (q−pb/2a)·∫1/(a·x²+b·x+c)` (reusing the arctangent
  rule for the second term). Dispatched right after `try_arctan_quadratic`, so
  it also produces the clean `log(x²+1)/2` for `∫x/(x²+1)`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[rational][arctan][regression]`.
- **Scope:** irreducible quadratic denominators with a degree-≤1 numerator and
  rational coefficients. Reducible / repeated-root denominators stay with
  `try_rational` (partial fractions → logs); higher-degree denominators such as
  `∫1/(x³+1)` were addressed by APART-1 / INT-17 below.

### APART-1 — `apart` did not decompose over irreducible quadratic factors
- **Input:** `apart(1/(x³+1))`, `apart(1/(x⁴-1))`.
- **Was:** the input fraction unchanged — `apart` only did Heaviside cover-up
  over distinct *linear* (rational-root) factors and bailed when an irreducible
  quadratic (or higher) factor was present.
- **Expected (SymPy):** `1/(3(x+1)) − (x−2)/(3(x²−x+1))`, etc.
- **Fix:** added `partial_fractions_squarefree` (`src/polys/poly.cpp`) —
  factor the denominator (`factor_list`), then solve the undetermined-
  coefficients identity `num = Σ Pᵢ·(den/fᵢ)` (an N×N rational linear system,
  N = deg den) via `Matrix::rref`, giving `Σ Pᵢ/fᵢ` with `deg Pᵢ < deg fᵢ`.
- **Regression test:** `tests/polys/poly_test.cpp` — `[apart][regression]`.
- **Scope:** squarefree denominators (distinct factors, multiplicity 1) with
  rational coefficients, `deg ≤ 10`. Repeated factors `(x−1)²` still defer.

### INT-17 — `integrate(1/(x³+1))` / rational with irreducible quadratic denominator
- **Input:** `∫1/(x³+1)`, `∫1/(x⁴-1)`.
- **Was:** the unevaluated marker — `try_rational` relied on `apart`, which
  could not split an irreducible quadratic factor (APART-1).
- **Expected (SymPy):** `log(x+1)/3 − log(x²−x+1)/6 + √3·atan(...)/3`, etc.
- **Fix:** with APART-1, `apart` now produces `(linear)/(irreducible quadratic)`
  terms which `integrate` closes via INT-16 (log + atan) plus the linear terms
  (logs). No integration code changed — the fix is entirely in `apart`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][rational][regression]`.

### INT-18 — `integrate(sin(x)**m * cos(x)**n)` returned the marker
- **Input:** `∫sin³`, `∫cos³`, `∫cos⁵`, `∫sin³cos²`, `∫sin²cos²`, `∫sin⁴`.
- **Was:** the marker — only `sin²`/`cos²` (single, power-2) had a reduction;
  higher powers and products fell through.
- **Expected (SymPy):** `cos³/3 − cos`, `sin − sin³/3`, `x/8 − sin(4x)/32`, etc.
- **Fix:** added `try_trig_power` (`src/integrals/integrate.cpp`) for
  `sin(g)^m·cos(g)^n` with an affine argument `g` and non-negative integer
  powers. An odd power uses the `u = sin`/`u = cos` substitution into a
  polynomial integral; the both-even case uses half-angle reduction
  (`sin² = (1−cos2g)/2`, `cos² = (1+cos2g)/2`) and recurses via `integrate`
  (degree strictly drops, so it terminates).
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][trig_power][regression]`.
- **Scope:** integer powers of sin/cos with a common affine argument. `tan`/
  hyperbolic powers and trig substitution (`∫√(1−x²)`) are separate items.

### INT-19 — `integrate(tan**n)` / `integrate(sinh**m·cosh**n)` returned the marker
- **Input:** `∫tan³`, `∫tan⁴`, `∫tan⁵`, `∫sinh²`, `∫cosh²`, `∫sinh³`,
  `∫cosh³`, `∫sinh²cosh²`.
- **Was:** the marker — only `tan²` (INT-8) and `sinh/cosh` of an affine
  argument at power 1 (INT-2) were handled; higher powers fell through.
- **Expected (SymPy):** `tan²/2 + log(cos)`, `x + tan³/3 − tan`,
  `sinh(2x)/4 − x/2`, `sinh(2x)/4 + x/2`, `cosh³/3 − cosh`, etc.
- **Fix (`src/integrals/integrate.cpp`):**
  - `try_tan_power`: `∫tanⁿ = tan^(n-1)/((n-1)·g') − ∫tan^(n-2)`, recursing
    through `integrate` down to the `∫tan` table case.
  - `try_hyperbolic_power`: the hyperbolic analogue of `try_trig_power` using
    `cosh²−sinh²=1` — an odd power uses `u = sinh`/`u = cosh` substitution,
    both-even uses the `cosh(2g)` half-angle forms and recurses.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[tan_power][regression]` and `[hyperbolic][regression]`.
- **Scope:** integer powers with a common affine argument. `tanh`/`coth`
  powers and trig substitution remain separate items.

### INT-20 — `integrate(sqrt(a·x²+c))` returned the marker
- **Input:** `∫√(1−x²)`, `∫√(4−x²)`, `∫√(x²+1)`, `∫√(2x²+3)`, `∫√(x²−1)`.
- **Was:** the marker — only the *reciprocal* `1/√(a·x²+c)` (INT-6/7) was
  handled; the radical in the numerator fell through.
- **Expected (SymPy):** `x·√(1−x²)/2 + asin(x)/2`, `x·√(4−x²)/2 + 2·asin(x/2)`,
  `x·√(x²+1)/2 + asinh(x)/2`, `x·√(2x²+3)/2 + 3√2·asinh(√6·x/3)/4`,
  `x·√(x²−1)/2 − log(x + √(x²−1))/2`.
- **Fix (`src/integrals/integrate.cpp`):** `try_sqrt_quadratic` now also matches
  exponent `+1/2`. Integration by parts gives
  `∫√(a·x²+c) = (x/2)·√(a·x²+c) + (c/2)·∫1/√(a·x²+c)`, so it reuses its own
  reciprocal branch (asin / asinh / log) for the second term. A `nullopt` inner
  integral (`c = 0`, or `a < 0` with `c ≤ 0` — no real region) propagates, so
  those still fall through unevaluated.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][invtrig][regression]` (five cases, verified by
  differentiation against the oracle).
- **Scope:** pure quadratic radicand (no linear term), rational coefficients.
  A linear term needs completing-the-square; `∫√(x+1)`-style algebraic u-subs
  are handled by INT-21.

### INT-21 — `integrate(P(x)·(a·x+b)^r)` returned the marker
- **Input:** `∫x·√(x+1)`, `∫x·√(2x+3)`, `∫x²·√(x+1)`, `∫x/√(x+1)`,
  `∫x·(x+1)^(1/3)`.
- **Was:** the marker — a polynomial times a fractional power of a linear
  expression had no rule (the table only handles a bare `(a·x+b)^r`, and
  `try_heurisch` couldn't pick the substitution).
- **Expected (SymPy):** e.g. `∫x·√(x+1) = 2x²√(x+1)/5 + 2x√(x+1)/15 −
  4√(x+1)/15` (SymPP prints the equivalent `−2(x+1)^(3/2)/3 + 2(x+1)^(5/2)/5`).
- **Fix (`src/integrals/integrate.cpp`):** `try_algebraic_linear_sub` matches a
  single `(affine)^(non-integer rational)` factor with the rest of the product
  polynomial in `var`. The substitution `u = a·x+b` (so `x = (u−b)/a`) turns the
  integrand into `Σ cₖ·u^(k+r)`, integrated term-by-term to
  `Σ cₖ·u^(k+r+1)/(k+r+1)` — `r ∉ ℤ` guarantees the denominator is never 0 —
  then back-substituted. Dispatched after `try_sqrt_quadratic`, before
  `try_heurisch`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][algebraic][regression]` (five cases incl. a negative exponent
  and a cube root, verified by differentiation against the oracle).
- **Scope:** one fractional power of an *affine* base. A fractional power of a
  *quadratic* (genuine trig/hyperbolic substitution) and products of two
  distinct algebraic radicals remain out of scope.

### INT-22 — `integrate(rational(exp(x)))` returned the marker
- **Input:** `∫1/(1+exp(x))`, `∫exp(x)/(1+exp(x))`, `∫1/(exp(x)−1)`,
  `∫1/(1+exp(2x))`.
- **Was:** the marker — `try_heurisch` *did* pick `g = exp(x)` and substitute
  it out, but its inner integration was table-only, so the resulting rational
  integrand (e.g. `1/(u·(1+u))`) was never decomposed into partial fractions.
- **Expected (SymPy):** `x − log(exp(x)+1)`, `log(exp(x)+1)`,
  `x − log(...)`-style log combinations (SymPP keeps `log(exp(x))` rather than
  folding it to `x`, but the antiderivatives are equal).
- **Fix (`src/integrals/integrate.cpp`):** in `try_heurisch`, when the table
  can't close the substituted integrand `q_sub`, fall back to
  `try_rational(q_sub, u)`. `try_rational` decomposes via `apart` into strictly
  simpler pieces (so it terminates), and the existing depth guard backstops its
  internal `integrate()` calls. This generalises beyond `exp`: any substitution
  that yields a rational function in `u` now closes.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][heurisch][regression]` (four cases, verified by differentiation
  against the oracle).
- **Scope:** integrands that become a *rational function* of the substituted
  variable. `∫1/(a+b·exp(x)+c·exp(2x))`-style cases work when `apart` can split
  the denominator.

### INT-23 — `integrate(P(x)·exp(a·x)·sin/cos(g·x))` returned the marker
- **Input:** `∫x·eˣ·sin(x)`, `∫x·eˣ·cos(x)`, `∫x²·eˣ·sin(x)`,
  `∫x·e^(2x)·sin(3x)`.
- **Was:** the marker — the pure cyclic case (`∫e·sin/cos`) bails once a
  polynomial factor makes its leftover non-constant, and the single-function
  by-parts bails because `u = x·sin(x)` is not a polynomial.
- **Expected (SymPy):** e.g. `∫x·eˣ·sin(x) = x·eˣ·sin(x)/2 − x·eˣ·cos(x)/2 +
  eˣ·cos(x)/2`.
- **Fix (`src/integrals/integrate.cpp`):** a new by-parts branch in
  `try_integration_by_parts` — when a `Mul` has both an `exp(affine)` and a
  `sin/cos(affine)` factor and the remaining factors form a polynomial `u`, take
  `dv = exp·trig` (antiderivative = the cyclic closed form) and `u = P(x)`.
  Differentiating `u` lowers its degree each step, so `∫v·u'` recurses down to
  the bare cyclic base case (the marker/depth guards backstop it). The product
  `v·u'` is `expand`ed so it splits over its `Add` and `integrate()` recurses
  per term.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][byparts][regression]` (four cases incl. degree-2 and non-unit
  exp/trig rates, verified by differentiation against the oracle).
- **Scope:** a single `exp(affine)` and a single `sin/cos(affine)` with a
  polynomial multiplier. `exp·sinh/cosh` (non-cyclic) and products of two trig
  factors remain separate.

### TRIG-RECIP — `sec`, `csc`, `cot` were not function types
- **Input:** `cot(pi/4)`, `sec(pi/3)`, `csc(pi/6)`, `cot(0)`, `sec(x)`,
  `diff(cot(x))`, `parse("csc(x)")`.
- **Was:** the parser turned `sec`/`csc`/`cot` into generic undefined-function
  nodes — no auto-evaluation, no derivatives, no exact values; SymPy results
  could only be matched after a manual `1/cos`-style rewrite.
- **Now:** three distinct function types (`Cot`/`Sec`/`Csc` in
  `functions/trigonometric.{hpp,cpp}`, enum values already reserved). Each
  factory folds exact values at rational multiples of π via the existing
  `cos_pi`/`sin_pi`/`tan_pi` tables, handles poles → `zoo`
  (`cot(0)=sec(π/2)=csc(0)`), parity (`cot`/`csc` odd, `sec` even), period, the
  inverse compositions (`cot(atan x)=1/x`, …), and numeric `Float` evalf.
  Derivatives: `cot'=-csc²`, `sec'=sec·tan`, `csc'=-csc·cot`. Parser + LaTeX
  printer (`\cot`/`\sec`/`\csc`) updated; `str()`/C/Octave fall back to the
  `name()` spelling (Octave/MATLAB have these natively).
- **Implementation note:** exact values use a `recip_value` helper that inverts
  a clean `coeff·√k` value by parts (`c⁻¹·k⁻¹ᐟ²`) so the radical stays
  rationalised; `cot` routes through `1/tan(rπ)` to avoid multiplying two equal
  radicals (`√2·√2`), which the Mul canonicaliser leaves unfolded.
- **Regression test:** `tests/functions/trigonometric_test.cpp`
  — `[trig][reciprocal]` (canonical angles, poles, parity, inverse comps,
  parse round-trip, derivatives, evalf — verified against the oracle).
- **Scope:** the antiderivatives `∫cot/sec/csc` ship in INT-24;
  `acot`/`asec`/`acsc` inverses are not yet added.

### HYP-RECIP — `coth`, `sech`, `csch` were not function types
- **Input:** `coth(0)`, `sech(0)`, `coth(oo)`, `coth(-x)`, `sech(acosh(x))`,
  `diff(coth(x))`, `parse("csch(x)")`.
- **Was:** the parser made generic undefined-function nodes — no
  auto-evaluation, no derivatives.
- **Now:** three distinct function types (`Coth`/`Sech`/`Csch` in
  `functions/hyperbolic.{hpp,cpp}`, enum values already reserved), the
  hyperbolic analogue of TRIG-RECIP. Factories handle the values at 0 and ±oo
  (`coth(0)=csch(0)=zoo`, `sech(0)=1`, `coth(±oo)=±1`, `sech(±oo)=csch(±oo)=0`),
  parity (`coth`/`csch` odd, `sech` even), inverse compositions
  (`coth(atanh x)=1/x`, …), and numeric `Float` evalf via `mpfr_coth`/`sech`/
  `csch`. Derivatives: `coth'=-csch²`, `sech'=-sech·tanh`, `csch'=-csch·coth`.
  Parser + LaTeX (`\coth`, `\operatorname{sech}`, `\operatorname{csch}`)
  updated; `str()` falls back to the `name()` spelling.
- **Regression test:** `tests/functions/hyperbolic_test.cpp`
  — `[3f][reciprocal]` (values/poles, parity, inverse comps, parse round-trip,
  derivatives, evalf — verified against the oracle).
- **Scope:** the antiderivatives `∫coth/sech/csch` ship in INT-26;
  `acoth`/`asech`/`acsch` inverses are not added.

### INT-26 — `integrate(coth/sech/csch)` and their squares returned the marker
- **Input:** `∫coth`, `∫sech`, `∫csch`, `∫sech²`, `∫csch²`, `∫coth²`,
  `∫coth(2x+1)`.
- **Was:** the marker — `coth/sech/csch` only became function types in
  HYP-RECIP, so the integration table had no entries.
- **Expected (SymPy):** `∫coth=log(sinh)`, `∫sech=atan(sinh)` (Gudermannian),
  `∫csch=log(tanh(x/2))`, `∫sech²=tanh`, `∫csch²=−coth`, `∫coth²=x−coth`.
- **Fix (`src/integrals/integrate.cpp`):**
  - `integrate_term` affine switch: `Coth→log(sinh)/a`, `Sech→atan(sinh)/a`,
    `Csch→log(tanh((ax+b)/2))/a`.
  - `try_trig_reduction` squares: `sech²(u)→1/cosh²(u)`, `csch²(u)→1/sinh²(u)`,
    `coth²(u)→1/sinh²(u)+1` (`coth²−csch²=1`), reusing the tabulated
    `1/cosh²`, `1/sinh²` cases (affine `u`).
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][reciprocal][hyperbolic][regression]` (seven cases incl. affine
  scaling, verified by differentiation against the oracle).
- **Scope:** singles and squares. `∫cothⁿ`/`∫sechⁿ`/`∫cschⁿ` for n ≥ 3 stay
  deferred (no `tanhⁿ` power handler exists either).

### INT-24 — `integrate(cot/sec/csc)` returned the marker
- **Input:** `∫cot(x)`, `∫sec(x)`, `∫csc(x)`, `∫cot(2x+1)`, `∫sec(3x)`.
- **Was:** the marker — `cot/sec/csc` only became real function types in
  TRIG-RECIP, so the integration table had no entries for them.
- **Expected (SymPy):** `∫cot=log(sin(x))`,
  `∫sec=(log(sin+1)−log(sin−1))/2`, `∫csc=(log(cos−1)−log(cos+1))/2`, each
  divided by the affine slope `a`.
- **Fix (`src/integrals/integrate.cpp`):** three new `case` labels in the
  `integrate_term` affine-function switch (alongside Sin/Cos/Tan), reusing the
  closed forms above with the `1/a` argument scaling.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][reciprocal][regression]` (five cases incl. affine arguments,
  verified by differentiation against the oracle).
- **Scope:** an affine argument `a·x+b`. Squares and `cotⁿ` powers ship in
  INT-25; `∫secⁿ`/`∫cscⁿ` for n ≥ 3 (by-parts reduction) and products remain
  separate items.

### INT-25 — reciprocal-trio powers `∫sec²/csc²/cot²` and `∫cotⁿ`
- **Input:** `∫sec²`, `∫csc²`, `∫cot²`, `∫cot³`, `∫cot⁴`, `∫sec(2x)²`.
- **Was:** the marker — only `sin²/cos²/tan²` had power rewrites; the reciprocal
  trio (added in TRIG-RECIP) had none.
- **Expected (SymPy):** `∫sec²=tan`, `∫csc²=−cot`, `∫cot²=−cot−x`,
  `∫cot³=−cot²/2−log(sin)`, etc.
- **Fix (`src/integrals/integrate.cpp`):**
  - `try_trig_reduction` squares: `sec²(u)→1/cos²(u)`, `csc²(u)→1/sin²(u)`,
    `cot²(u)→1/sin²(u)−1` (Pythagorean), reusing the existing `1/cos²`,
    `1/sin²` table cases (affine `u`).
  - `try_tan_power` generalised to `Cot`: `∫cotⁿ = −cot^(n-1)/((n-1)a) −
    ∫cot^(n-2)`, the sign-flipped analogue of `∫tanⁿ`, recursing to the `∫cot`
    table case.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][reciprocal][regression]` (six cases incl. cubes/quartics and
  affine scaling, verified by differentiation against the oracle).
- **Scope:** `sec²/csc²/cot²` and integer `cotⁿ`. `∫secⁿ`/`∫cscⁿ` for n ≥ 3
  ship in INT-27.

### INT-27 — `∫secⁿ` / `∫cscⁿ` (n ≥ 3) returned the marker
- **Input:** `∫sec³`, `∫sec⁴`, `∫csc³`, `∫csc⁴`, `∫sec(2x)³`.
- **Was:** the marker — only the `n = 1` table case (INT-24) and `n = 2` square
  (INT-25) were handled; higher powers of `sec`/`csc` have no simple Pythagorean
  reduction (unlike `tan`/`cot`) and need integration by parts.
- **Expected (SymPy):** e.g. `∫sec³ = sec·tan/2 + log(sec+tan)/2` (SymPP emits
  the equivalent `sin/cos` log form).
- **Fix (`src/integrals/integrate.cpp`):** new `try_sec_csc_power` with the
  by-parts reduction, recursing to the `∫sec`/`∫sec²` base cases:
  - `∫secⁿ =  sec^(n-2)·tan/((n-1)·g') + (n-2)/(n-1)·∫sec^(n-2)`
  - `∫cscⁿ = −csc^(n-2)·cot/((n-1)·g') + (n-2)/(n-1)·∫csc^(n-2)`
  Dispatched right after `try_tan_power`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][reciprocal][regression]` (five cases incl. affine scaling,
  verified by differentiation against the oracle).
- **Scope:** integer `secⁿ`/`cscⁿ`. The hyperbolic `sechⁿ`/`cschⁿ` analogues
  ship in INT-28.

### INT-28 — `∫sechⁿ` / `∫cschⁿ` (n ≥ 3) returned the marker
- **Input:** `∫sech³`, `∫sech⁴`, `∫csch³`, `∫csch⁴`, `∫sech(2x)³`.
- **Was:** the marker — only `n = 1` (table) and `n = 2` (square) of `sech`/
  `csch` were handled (INT-26); higher powers need integration by parts.
- **Expected (SymPy):** e.g. `∫sech³ = sech·tanh/2 + atan(sinh)/2`.
- **Fix (`src/integrals/integrate.cpp`):** new `try_sech_csch_power`, the
  hyperbolic analogue of `try_sec_csc_power`. The Pythagorean sign differs
  (`coth² − csch² = 1`), so the `csch` rest term is **subtracted**:
  - `∫sechⁿ =  sech^(n-2)·tanh/((n-1)·g') + (n-2)/(n-1)·∫sech^(n-2)`
  - `∫cschⁿ = −csch^(n-2)·coth/((n-1)·g') − (n-2)/(n-1)·∫csch^(n-2)`
  Recurses to the `∫sech`/`∫sech²` base cases (INT-26); dispatched after
  `try_sec_csc_power`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][reciprocal][hyperbolic][regression]` (five cases incl. affine
  scaling, verified by differentiation against the oracle).
- **Scope:** integer `sechⁿ`/`cschⁿ`. With INT-27 this closes the
  reciprocal-power integration family (trig + hyperbolic, all six functions).

### GAMMA-1 — `gamma` at a half-integer stayed symbolic
- **Input:** `gamma(1/2)`, `gamma(3/2)`, `gamma(5/2)`, `gamma(7/2)`,
  `gamma(-1/2)`, `gamma(-3/2)`.
- **Was:** `gamma(1/2)`, … — only positive *integer* arguments reduced (to
  `(n-1)!`); rational arguments fell straight through to the symbolic node.
- **Expected (SymPy):** `sqrt(pi)`, `sqrt(pi)/2`, `3*sqrt(pi)/4`,
  `15*sqrt(pi)/8`, `-2*sqrt(pi)`, `4*sqrt(pi)/3`.
- **Fix:** in `gamma` (`src/functions/combinatorial.cpp`), a `Rational` with
  denominator 2 reduces to the base `gamma(1/2) = sqrt(pi)` via the recurrence
  `gamma(z) = (z-1)·gamma(z-1)` (and its inverse `gamma(z) = gamma(z+1)/z` for
  `z < 1/2`), accumulating an exact rational coefficient
  (`half_integer_gamma_coeff`). The numerator is bounded (±100001) so the
  recurrence can never spin.
- **Regression test:** `tests/functions/combinatorial_test.cpp`
  — `[gamma][regression]`.

### SQRT-2 — `sqrt` did not extract square factors or rationalise
- **Input:** `sqrt(8)`, `sqrt(12)`, `sqrt(rational(1,2))`,
  `sqrt(rational(2,3))`, `sqrt(rational(8,9))`.
- **Was:** `8**(1/2)`, `(1/2)**(1/2)`, … — left fully under the root.
- **Expected (SymPy):** `2*sqrt(2)`, `2*sqrt(3)`, `sqrt(2)/2`, `sqrt(6)/3`,
  `2*sqrt(2)/3`.
- **Fix:** added `try_sqrt_factor_extraction` in `src/core/pow.cpp`, dispatched
  after `try_perfect_root` (SQRT-1). It pulls the largest square factor out of
  the radicand and, for a rational `p/q`, rationalises via
  `√(p/q) = √(p·q)/q`. Trial division is bounded (radicand ≤ 1e12) so a huge
  radicand can never hang.
- **Regression test:** `tests/functions/miscellaneous_test.cpp`
  — `[sqrt][regression]` (incl. a prime-radicand `√7` no-op guard).
- **Scope:** square roots only — n-th-root factor extraction (`cbrt(16) →
  2·cbrt(2)`) is generalised in NROOT-1.

### NROOT-1 — `N^(1/n)` did not extract n-th-power factors (n ≥ 3)
- **Input:** `16^(1/3)`, `54^(1/3)`, `24^(1/3)`, `96^(1/5)`, `48^(1/4)`,
  `(2/3)^(1/3)`, `(16/27)^(1/3)`.
- **Was:** only perfect n-th roots (`8^(1/3)=2`) and the square case (SQRT-2)
  folded; a non-perfect cube/n-th root stayed fully under the radical.
- **Expected (SymPy):** `2·2^(1/3)`, `3·2^(1/3)`, `2·3^(1/3)`, `2·3^(1/5)`,
  `2·3^(1/4)`, `18^(1/3)/3`, `2·2^(1/3)/3`.
- **Fix (`src/core/pow.cpp`):** generalised `try_sqrt_factor_extraction` to
  `try_nth_root_factor_extraction` — for a unit `1/n` power it pulls the largest
  `sⁿ` factor (`N = sⁿ·m → s·m^(1/n)`) and rationalises a rational radicand via
  `(p/q)^(1/n) = (p·q^(n-1))^(1/n)/q`. Trial division stays bounded (radicand ≤
  1e12) so a huge radicand can never hang.
- **Regression test:** `tests/functions/miscellaneous_test.cpp`
  — `[sqrt][nthroot][regression]` (incl. a prime-radicand `7^(1/3)` no-op guard).
- **Scope:** unit `1/n` exponents on a non-negative Integer/Rational base.
  Non-unit non-perfect powers (`16^(2/3)`) ship in POW-RAT-2.

### POW-RAT-2 — `b^(p/q)` (p ≥ 2, non-perfect base) not simplified
- **Input:** `16^(2/3)`, `2^(5/2)`, `12^(2/3)`, `48^(5/2)`, `72^(2/3)`,
  `7^(3/2)`.
- **Was:** only perfect q-th powers (`8^(2/3)=4`, POW-RAT) and unit numerators
  (`16^(1/3)`, NROOT-1) reduced; a non-unit power of a non-perfect base stayed
  fully symbolic (`16^(2/3)` → `16**(2/3)`).
- **Expected (SymPy):** `4·2^(2/3)`, `4√2`, `2·2^(1/3)·3^(2/3)`, `9216√3`,
  `12·3^(1/3)`, `7√7`.
- **Fix (`src/core/pow.cpp`):** `try_rational_power_extraction` prime-factorises
  `b = ∏ pᵢ^aᵢ`; each prime's power exponent `aᵢ·p/q` splits into an integer
  part `⌊aᵢp/q⌋` (pulled into the integer coefficient) and a remainder
  `rᵢ = aᵢp mod q` (kept under a per-prime radical `pᵢ^(rᵢ/q)`). Keeping primes
  *separate* matches SymPy's canonical form (`16^(2/3) = 4·2^(2/3)`, not the
  equivalent `4·4^(1/3)`). The residual `pow()` factors are built only after the
  "something was pulled" check, so the recursive call on a bare prime
  (`2^(2/3)`) bails first and can't recurse without bound. Factorisation is
  trial-division bounded (base ≤ 1e12, numerator ≤ 1000).
- **Regression test:** `tests/core/arithmetic_test.cpp`
  — `[pow][regression]` (incl. a `2^(2/3)` no-op guard).
- **Scope:** positive integer base, `p ≥ 2`. Rational bases and negative
  exponents (`16^(-2/3)`) stay deferred.

### SQRT-3 — `sqrt` of a negative number not folded to imaginary
- **Input:** `sqrt(-1)`, `sqrt(-4)`, `sqrt(-8)`, `sqrt(-1/4)`, `sqrt(-2/3)`.
- **Was:** `(-1)**(1/2)`, … — SQRT-1/SQRT-2 deferred negative bases for
  branch handling.
- **Expected (SymPy):** `I`, `2*I`, `2*sqrt(2)*I`, `I/2`, `sqrt(6)*I/3`.
- **Fix:** added `try_sqrt_of_negative` in `src/core/pow.cpp` — for the ½ power
  of a negative Integer/Rational, returns `I·√|base|`, reusing the
  perfect-root / factor-extraction paths so the magnitude comes back fully
  reduced.
- **Regression test:** `tests/functions/miscellaneous_test.cpp`
  — `[sqrt][regression]`.
- **Scope:** only the principal square root (½ power); other fractional powers
  of a negative base (`(-8)^(1/3)`) need full branch-cut handling and stay
  symbolic.

### SQRT-4 — inverse square root left the denominator irrational
- **Input:** `3**(-1/2)`, `2**(-1/2)`, `12**(-1/2)`, `(2/3)**(-1/2)`,
  `atan(1/sqrt(3))`.
- **Was:** `3**(-1/2)`, … — `try_sqrt_factor_extraction` only handled the `+½`
  power, so a `−½` power stayed unrationalised. SymPy rationalises
  `n^(-1/2) = sqrt(n)/n`.
- **Expected (SymPy):** `sqrt(3)/3`, `sqrt(2)/2`, `sqrt(3)/6`, `sqrt(6)/2`,
  and the knock-on `atan(1/sqrt(3)) = pi/6`.
- **Fix:** added `try_inverse_sqrt` in `src/core/pow.cpp` — for a `−½` power of
  a positive Integer/Rational, returns `pow(reciprocal, 1/2)`, which the
  existing `+½` factor-extraction path then rationalises. The inverse-trig
  table already recognises the resulting `sqrt(3)/3` form, so `atan(1/sqrt(3))`
  now folds to `pi/6`.
- **Regression test:** `tests/functions/miscellaneous_test.cpp`
  — `[sqrt][regression]`.
- **Scope:** the direct power form `n^(-1/2)`. `1/sqrt(12)` — where `sqrt(12)`
  first evaluates to `2·sqrt(3)`, leaving a product `1/(2·sqrt(3))` — needs
  `radsimp`-style product-denominator rationalisation and stays as written.

### POW-RAT — `a^(p/q)` of a perfect power was not folded (p ≠ 1)
- **Input:** `8**(2/3)`, `16**(3/4)`, `4**(3/2)`, `32**(2/5)`,
  `(8/27)**(2/3)`, `8**(-2/3)`.
- **Was:** `8**(2/3)`, … — `try_perfect_root` bailed unless the exponent
  numerator was 1, so only `1/q` roots (`27**(1/3)=3`) folded.
- **Expected (SymPy):** `4`, `8`, `8`, `4`, `4/9`, `1/4`.
- **Fix:** `try_perfect_root` (`src/core/pow.cpp`) now handles any numerator
  `p`: it takes the exact `q`-th root of a non-negative Integer/Rational base
  and raises it to `p` (`a^(p/q) = (a^(1/q))^p`); `pow` then folds the
  integer/rational power (negative `p` → Rational). A non-exact root stays an
  irreducible `Pow`.
- **Regression test:** `tests/core/arithmetic_test.cpp` — `[pow][regression]`.
- **Scope:** the base must be a perfect `q`-th power. `12**(2/3)` (no exact
  cube root) stays symbolic — SymPy extracts `2*18**(1/3)` via n-th-root factor
  extraction, a separate feature not yet implemented.

### LOG-1 — `log` of a negative / imaginary argument not evaluated
- **Input:** `log(-1)`, `log(-2)`, `log(-E)`, `log(I)`, `log(-I)`, `log(2*I)`.
- **Was:** `log(-1)`, … — left unevaluated.
- **Expected (SymPy):** `I*pi`, `log(2) + I*pi`, `1 + I*pi`, `I*pi/2`,
  `-I*pi/2`, `log(2) + I*pi/2`.
- **Fix:** in the `log` factory (the inverse of EXP-1): `log(x) = log(|x|) + Iπ`
  for a negative real `x` (guarded by `is_real`/`is_negative`, so it also folds
  `−E`), and `log(b·I) = log(|b|) + sign(b)·Iπ/2` for a nonzero rational `b`
  via an `imaginary_coeff` helper.
- **Regression test:** `tests/functions/exponential_test.cpp`
  — `[log][regression]` (negative reals, imaginary axis, and a positive/symbolic
  no-op guard).
- **Scope:** principal branch; general complex `log(a+b·I)` (off the axes) is
  not auto-evaluated, matching SymPy.

### ABS-1 — `Abs(c·x)` did not pull out a numeric coefficient
- **Input:** `abs(-2*x)`, `abs(2*x)`, `abs(x/2)`, `abs(-x/3)`, `abs(-2*x*y)`.
- **Was:** `Abs(-2*x)`, … — only a leading `−1` was stripped (`Abs(-x)=Abs(x)`);
  any other numeric coefficient stayed inside.
- **Expected (SymPy):** `2*Abs(x)`, `2*Abs(x)`, `Abs(x)/2`, `Abs(x)/3`,
  `2*Abs(x*y)`.
- **Fix:** in the `abs` factory, `Abs(c·rest) = |c|·Abs(rest)` for a numeric
  leading factor `c` (canonical Mul sorts a number first). Subsumes the old
  `−1` rule and matches SymPy. `|·|` is multiplicative, so it is valid for any
  coefficient.
- **Regression test:** `tests/functions/miscellaneous_test.cpp`
  — `[abs][regression]`.
- **Scope:** the imaginary unit `I` sorts last in a Mul, so `Abs(I·x)` stays
  `Abs(x·I)` (correct, equal to `Abs(x)`, just not folded). `Sign`/`Re`/`Im`
  keep their existing `−1`-only handling.

### ABS-2 — `Abs` of a numeric complex number stayed symbolic
- **Input:** `Abs(3+4*I)`, `Abs(1+I)`, `Abs(2+3*I)`, `Abs(2*I)`, `Abs(I)`,
  `Abs(-3-4*I)`.
- **Was:** `Abs(4*I + 3)`, … — only real numbers reduced; a complex literal
  fell through to the symbolic node.
- **Expected (SymPy):** `5`, `sqrt(2)`, `sqrt(13)`, `2`, `1`, `5`.
- **Fix:** in the `abs` factory, a value that parses as `a + b·I` with rational
  real and imaginary parts (`rational_complex`) returns the modulus
  `sqrt(a² + b²)`. `rational_imag_coeff` extracts the coefficient of a
  pure-imaginary term; the existing `sqrt` then reduces perfect squares
  (`sqrt(25)=5`). Purely real / symbolic inputs are untouched.
- **Regression test:** `tests/functions/miscellaneous_test.cpp`
  — `[abs][regression]`.
- **Scope:** rational real/imaginary parts only — a symbolic or irrational
  component (`Abs(x+I)`, `Abs(sqrt(2)+I)`) stays unevaluated.

### ATAN2-1 — `atan2` only reduced on the axes
- **Input:** `atan2(1,1)`, `atan2(-1,1)`, `atan2(1,-1)`, `atan2(-1,-1)`,
  `atan2(1,sqrt(3))`, `atan2(2,1)`.
- **Was:** `atan2(1, 1)`, … — only the axis cases (`y=0` or `x=0`) reduced; a
  general quadrant stayed unevaluated.
- **Expected (SymPy):** `pi/4`, `-pi/4`, `3*pi/4`, `-3*pi/4`, `pi/6`, `atan(2)`.
- **Fix:** in `atan2` (`src/functions/trigonometric.cpp`), when `x` has a known
  sign and `y` is real, rewrite `atan2(y, x) = atan(y/x)` with a quadrant
  correction (`+pi` for `x<0, y≥0`; `-pi` for `x<0, y<0`). `atan` then folds the
  special values (`atan(1)=pi/4`, `atan(sqrt(3))=pi/3`). The rewrite is faithful
  even when `atan` cannot fold the argument (`atan2(2,1)=atan(2)`).
- **Regression test:** `tests/functions/inverse_trig_test.cpp`
  — `[atan2][regression]`.
- **Scope:** applies when `x`'s sign is decidable and `y` is real; fully
  symbolic arguments stay unevaluated.

### REIM-1 — `re`/`im`/`conjugate` of a numeric complex stayed unevaluated
- **Input:** `re(3+4*I)`, `im(3+4*I)`, `conjugate(3+4*I)`, `conjugate(2*I)`,
  `conjugate(1/2+I/3)`.
- **Was:** `re(4*I + 3)`, … — only real arguments reduced (`re(x)=x` for real
  `x`); a numeric complex fell through (the code noted "no Complex type yet").
- **Expected (SymPy):** `3`, `4`, `3 - 4*I`, `-2*I`, `1/2 - I/3`.
- **Fix:** `re`/`im`/`conjugate` (`src/functions/miscellaneous.cpp`) now use the
  `rational_complex` helper (shared with ABS-2) to split `a + b·I` into rational
  real/imaginary parts: `re → a`, `im → b`, `conjugate → a − b·I`. Real and
  symbolic arguments are untouched.
- **Regression test:** `tests/functions/miscellaneous_test.cpp`
  — `[complex][regression]`.
- **Scope:** rational real/imaginary parts only — a symbolic component
  (`re(x+I)`) stays unevaluated.

### BINOM-1 — `binomial(n, 1)` not simplified to `n`
- **Input:** `binomial(n, 1)`.
- **Was:** `binomial(n, 1)` — kept symbolic (only `binomial(n,0)=1` and the
  numeric / `n==n` cases were handled).
- **Expected (SymPy):** `n` (valid for any `n`).
- **Fix:** added `binomial(n, 1) = n` to the factory.
- **Regression test:** `tests/functions/combinatorial_test.cpp`
  — `[binomial][regression]` (incl. a `binomial(n,2)` stays-symbolic guard).
- **Note:** `binomial(n, n)` for a plain symbol is *not* auto-simplified — SymPy
  keeps it too, so SymPP matches by leaving it (it only folds when `n` is a
  known nonnegative integer).

### PARSE-1 — parser rejected the capitalised names `str()` emits
- **Input:** `parse("Abs(-3)")`, `parse(abs(x)->str())` (= `parse("Abs(x)")`),
  same for `Heaviside`, `DiracDelta`.
- **Was:** an *undefined function* `Abs(...)` — the parser table held only the
  lowercase aliases (`abs`, `heaviside`, `dirac_delta`), but `str()` prints the
  SymPy-canonical capitalised names, so `parse(e->str())` did not round-trip.
- **Expected (SymPy):** `Abs(-3) → 3`, and `parse(e->str()) == e`.
- **Fix:** added `Abs` / `Heaviside` / `DiracDelta` aliases to the parser's
  one-argument function table (the lowercase spellings still work).
- **Regression test:** `tests/parsing/parser_test.cpp`
  — `[parser][regression]` (capital-name eval + str round-trip).
- **Scope:** the other functions (`sign`, `floor`, `re`, `im`, `conjugate`,
  `gamma`, `erf`, …) already print lowercase, matching both the parser and
  SymPy, so they round-trip unchanged.

### PARSE-2 — parser did not recognise `Min`/`Max`
- **Input:** `parse("Min(3, 5)")`, `parse(min(x,y)->str())` (= `"Min(x, y)"`).
- **Was:** an undefined function `Min(...)` — the parser's two-argument table
  had no `Min`/`Max` entry, so `parse(e->str())` did not round-trip.
- **Expected (SymPy):** `Min(3,5) → 3`, `Max(3,5) → 5`, and `parse(e->str()) == e`.
- **Fix:** added `Min`/`Max` (the names `str()` emits) to the two-argument
  function table, bound to the binary `min`/`max` overloads.
- **Regression test:** `tests/parsing/parser_test.cpp`
  — `[parser][regression]`.
- **Scope:** the binary form only — 3+-argument `Min`/`Max` were addressed in
  PARSE-3 below.

### PARSE-3 — parser did not fold 3+-argument `Min`/`Max`
- **Input:** `parse("Max(3, 7, 1)")`, `parse("Min(3, 7, 1)")`,
  `parse("Max(1, 2, 3, 4)")`, `parse("Max(x, 3, 1)")`.
- **Was:** `Max(3, 7, 1)`, … — the parser only dispatched 1- and 2-argument
  `Min`/`Max`; a 3+-argument call fell through to an undefined function and
  stayed unevaluated, even though the variadic `min`/`max` already fold.
- **Expected (SymPy):** `7`, `1`, `4`, `Max(3, x)`.
- **Fix:** `apply_function` now routes any-arity `Min`/`Max` to the variadic
  `min(args)`/`max(args)` (which combine the numeric args into one extreme and
  keep the symbolic ones), instead of only the 2-argument table entries
  (`src/parsing/parser.cpp`).
- **Regression test:** `tests/parsing/parser_test.cpp`
  — `[parser][regression]`.

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

### LIM-1 — limits at infinity returned garbage; no `Infinity` type ([#2])
- **Input:** `limit((1 + 1/x)**x, x, oo)`, `limit(x**2/(x+1), x, oo)`,
  `limit(exp(x)/x, x, oo)`, `limit(x*sin(1/x), x, oo)`.
- **Was:** `(oo**(-1) + 1)**oo`, … — `oo` parsed as a plain symbol and the
  limit engine only did finite-point direct substitution + 0/0 L'Hôpital, so
  every limit at infinity was wrong.
- **Expected (SymPy):** `E`, `oo`, `oo`, `1`.
- **Fix:** added real `Infinity` / `NegativeInfinity` / `ComplexInfinity` /
  `NaN` atoms (`src/core/infinity.hpp`/`.cpp`) with `S::Infinity()` etc., wired
  them through `add`/`mul`/`pow` (oo+finite=oo, oo-oo=nan, oo*0=nan, 1/oo=0,
  2^oo=oo, 1^oo=nan, …) and through `exp`/`log` (exp(oo)=oo, exp(-oo)=0,
  log(oo)=oo, log(0)=zoo). The parser maps `oo`/`zoo`/`nan` (and `-oo`). The
  limit engine (`src/calculus/limit.cpp`) now resolves the indeterminate forms:
  `1^∞`/`∞^0`/`0^0` via `a^b = exp(b·log a)`, `0·∞` by rewriting to a `0/0`
  quotient, and `∞/∞` (and `0/0`) by L'Hôpital with `together()`-based
  re-rationalisation each step.
- **Regression test:** `tests/core/infinity_test.cpp` (`[infinity]`) and
  `tests/calculus/series_limit_test.cpp` (`[limit][infinity][regression]`).
- **Scope:** still unresolved — a finite-point pole (`limit(1/x**2, x, 0)`,
  needs one-sided handling) and bare `∞ − ∞` polynomial forms
  (`limit(x - x**2, x, oo)`, needs dominant-term extraction); both stay
  unevaluated rather than returning a wrong value.

### FUNC-INF — elementary functions did not evaluate at ±oo
- **Input:** `atan(oo)`, `tanh(oo)`, `sinh(oo)`, `cosh(-oo)`, `asinh(-oo)`,
  `acosh(oo)`, `erf(oo)`, `erfc(-oo)`.
- **Was:** `atan(oo)`, … — left unevaluated (the builders only handled finite
  arguments), so e.g. `limit(atan(x), x, oo)` returned `atan(oo)`.
- **Expected (SymPy):** `pi/2`, `1`, `oo`, `oo`, `-oo`, `oo`, `1`, `2`.
- **Fix:** added the infinite-argument limits to the function factories —
  `atan(±oo)=±pi/2` (`trigonometric.cpp`); `sinh(±oo)=±oo`, `cosh(±oo)=oo`,
  `tanh(±oo)=±1`, `asinh(±oo)=±oo`, `acosh(±oo)=oo` (`hyperbolic.cpp`);
  `erf(±oo)=±1`, `erfc(oo)=0`, `erfc(-oo)=2` (`special.cpp`). This also makes
  the corresponding limits resolve directly (`limit(atan(x),x,oo)=pi/2`).
- **Regression test:** `tests/functions/{inverse_trig,hyperbolic,special}_test.cpp`
  and `tests/calculus/series_limit_test.cpp` (`[infinity][regression]`).
- **Scope:** oscillatory `sin(oo)`/`cos(oo)` stay unevaluated (no real limit —
  SymPy returns `AccumBounds`, not modeled here).

### EVALF-1 — `evalf` did not recurse into Add/Mul/Pow/Function
- **Input:** `evalf(2*pi)`, `evalf(sqrt(2))`, `evalf(sin(pi))`, `evalf(exp(1))`.
- **Was:** `2*pi`, `2**(1/2)`, … — the `evalf` switch only handled numeric
  atoms (Integer/Rational/Float/NumberSymbol/RootOf); compound expressions hit
  a `default` identity branch (a `// Phase 1f` TODO), so a numeric constant
  expression never reduced to a `Float`.
- **Expected (SymPy):** `6.2831853…`, `1.4142135…`, `≈0`, `2.7182818…`.
- **Fix:** `evalf` (`src/core/float.cpp`) now recurses — Add/Mul/Pow rebuild
  through `add`/`mul`/`pow` over evalf'd args (numeric folding does the rest),
  and a Function evalfs its arguments then `rebuild`s (a now-Float argument
  folds via the function's numeric path, e.g. `sin(pi) → sin(Float) ≈ 0`).
- **Regression test:** `tests/core/float_test.cpp` — `[evalf][regression]`.

### FLOOR-CONST — `floor`/`ceiling` of a real constant stayed symbolic
- **Input:** `floor(pi)`, `ceiling(pi)`, `floor(2*pi)`, `floor(-pi)`,
  `floor(sqrt(2))`, `floor(pi**2)`.
- **Was:** `floor(pi)`, … — only Integer/Rational/Float (and integer-tagged
  symbols) folded; a symbolic real constant stayed wrapped.
- **Expected (SymPy):** `3`, `4`, `6`, `-4`, `1`, `9`.
- **Fix:** `floor`/`ceiling` (`src/functions/integers.cpp`) now evalf a
  free-symbol-free argument (enabled by EVALF-1) and round the resulting Float
  to an exact Integer, with a boundary guard that refuses to fold when the
  value sits within ~1e-40 of an integer (so a disguised integer cannot be
  mis-rounded). A complex (`floor(I)`) or infinite (`floor(oo)`) argument does
  not evalf to a Float and is left unevaluated.
- **Regression test:** `tests/functions/integers_test.cpp`
  — `[floor][ceiling][regression]`.

### CANCEL-1 — `cancel()`/`Poly` GCD hung on symbolic coefficients ([#5])
- **Input:** `cancel((b - a + 1)*(a + b)/2, a)`, `factor(x**2 - y**2, x)` (and
  any polynomial whose coefficients in the working variable are symbolic).
- **Was:** infinite loop — never returned.
- **Expected (SymPy):** the reduced expression (`cancel`), or — for true
  multivariate factorization — a factored/unfactored result.
- **Cause:** in `Poly::divmod` the leading-term cancellation
  `lc − (lc/lc_other)·lc_other` did not fold to a structural zero when the
  coefficients were symbolic Adds: `(b+b²) − (b+b²)` stayed an unmerged Add
  (the bare Add flattens but the `−1·Add` subtrahend does not), so the
  remainder degree never dropped and the Euclidean GCD spun forever.
- **Fix:** `divmod` now `expand`s each coefficient subtraction (so the
  cancellation folds to `0`) and has a degree-decrease backstop that stops if
  a coefficient cannot cancel. `cancel` is therefore safe on multivariate
  input — `cancel((x²−y²)/(x−y), x) = x+y`. `factor` is `ℚ`-coefficient only,
  so it now bails to the unfactored input when a coefficient is symbolic
  rather than entering the integer-coefficient machinery.
- **Regression test:** `tests/polys/poly_test.cpp`
  — `[cancel][factor][regression]`.
- **Scope:** `cancel` reduces multivariate fractions; `simplify` still applies
  cancel only in the univariate case (auto-applying it multivariate regressed
  a downstream ODE form — a separate quality task). True multivariate
  *factorization* (`x²−y² → (x−y)(x+y)`) is not yet implemented.

> All previously-filed parity bugs are now resolved. The earlier
> "Missing auto-simplifications" cluster (`exp(log(x))`, `sqrt(8)`,
> `(-4)**(1/2)`) was closed by EXP-1 / SQRT-2 / SQRT-3 and is no longer open.

[#1]: https://github.com/leonardoaraujosantos/SymPP/issues/1
[#2]: https://github.com/leonardoaraujosantos/SymPP/issues/2
[#3]: https://github.com/leonardoaraujosantos/SymPP/issues/3
[#4]: https://github.com/leonardoaraujosantos/SymPP/pull/4
[#5]: https://github.com/leonardoaraujosantos/SymPP/issues/5
[#7]: https://github.com/leonardoaraujosantos/SymPP/issues/7
[#9]: https://github.com/leonardoaraujosantos/SymPP/issues/9
[#11]: https://github.com/leonardoaraujosantos/SymPP/issues/11
[#13]: https://github.com/leonardoaraujosantos/SymPP/issues/13
