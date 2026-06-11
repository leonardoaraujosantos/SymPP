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
- **Scope:** denominators outside {1,2,3,4,6} (e.g. `π/12`, `π/5`) stay
  symbolic — SymPy expands those to nested radicals, not yet implemented.

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
- **Scope:** `coth`/`sech`/`csch` are not distinct function types, so results
  are spelled with `cosh`/`sinh`.

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
- **Scope:** square roots only — n-th-root factor extraction
  (`cbrt(16) → 2·cbrt(2)`) is not yet implemented. Radicands above the trial-
  division bound stay symbolic.

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
