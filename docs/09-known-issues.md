# Known issues & SymPy-parity bug log

This file tracks correctness gaps found by cross-checking SymPP against
the SymPy oracle (validated against 1.13.x and 1.14; CI pins 1.14, and a
version guard in the smoke tests fails loudly on any unvetted release). Each
entry records the failing input, the
SymPy-expected result, current status, and the regression test (once
fixed). Each entry is filed as a GitHub issue; the local log is the source of
truth and links the issue number.

> Networking note: this environment's DNS resolves `api.github.com` and
> `github.com` to a dead Azure block (`4.228.31.x`), so `gh` and plain
> `curl` to the API time out. Git push over HTTPS still works. To reach
> the REST API, force a real GitHub IP, e.g.
> `curl --resolve api.github.com:443:140.82.121.6 ...`.

## Fixed

### HYP-DBLRATIO-1 — `simplify(sinh(2x)/sinh(x))` did not reduce to `2·cosh(x)`
- **Problem:** the circular double-angle ratio reduction (TRIG-DBLRATIO-1) had no
  hyperbolic counterpart, so `sinh(2x)/sinh(x)`, `sinh(2x)/cosh(x)` and the `csch`/`sech`
  forms stayed unevaluated where SymPy returns `2·cosh(x)` / `2·sinh(x)`.
- **Fix:** added `hyp_double_angle_ratio_mul` in `src/simplify/simplify.cpp`, run in the
  `trigsimp_node` Mul pipeline after `hyp_ratio_mul`. It mirrors the circular handler with
  `sinh(2u) = 2·sinh(u)·cosh(u)` (no sign flip): for a `sinh(a)` factor it cancels a
  `1/sinh(u)` or `1/cosh(u)` denominator (negative power or `csch`/`sech`) when `a = 2u`,
  giving `sinh(2u)/sinh(u) → 2·cosh(u)`, `sinh(2u)/cosh(u) → 2·sinh(u)`. Coefficients carry
  through and deeper doubling cancels one level (`sinh(4x)/sinh(2x) → 2·cosh(2x)`).
  `cosh(2u)/cosh(u)` and `sinh(3u)/sinh(u)` are left alone, matching SymPy.

### EXPAND-LOG-FRACPOW-1 — `log(√x)` was not pulled to `log(x)/2` for a generic base
- **Problem:** `expand`/`simplify` extracted `log(bᵉ) → e·log(b)` only when the base `b`
  was known positive. For a generic symbol, `log(√x)`, `log(x^(2/3))` stayed unexpanded,
  whereas SymPy returns `log(x)/2`, `2·log(x)/3` even without a positivity assumption.
- **Fix:** in `expand_log_arg` (`src/core/expand.cpp`), added the branch-safe case: when the
  exponent is a rational `e` with `−1 < e < 1` (`e ≠ 0`), extract `e·log(b)` regardless of
  the base's sign. That range is exactly where `e·arg(b) ∈ (−π, π)` keeps the identity on
  the principal branch (the same bound SymPy uses). `log(√x) → log(x)/2`,
  `log(x^(1/3)) → log(x)/3`, `log(x^(−1/2)) → −log(x)/2`. Exponents with `|e| ≥ 1`
  (`log(x²)`, `log(1/x)` at the `e = −1` boundary, `log(x^(3/2))`) and symbolic exponents
  are left intact, matching SymPy. `simplify` picks it up through its `expand` step.

### SOLVE-CPLXFORM-1 — complex polynomial roots came back as `½·(…)` not `a + b·I`
- **Problem:** Cardano (and the quadratic formula) build a complex root as a rational
  prefactor times a sum — `½·(2·I − 2)`, `1/16·(4·I·√3 − 4)` — and `solve` returned it
  undistributed. The values were correct but read nothing like SymPy's `−1 + I`,
  `−1 + √3·I`. `solve(x³−8)` gave `[2, ½·(2I√3−2), …]` where SymPy gives `[2, −1+√3·I, …]`.
- **Fix:** in `src/solvers/solve.cpp`, after the var-filter and before dedup, map each
  root through `expand` (RootOf exempt — it embeds its defining polynomial). `expand`
  distributes the rational prefactor and collapses the factor of two, so a complex root
  reads as `a + b·I` and a real root as its distributed surd (`½·√5 + ½`). Done before
  dedup so roots differing only by distribution collapse to one. All solve tests assert by
  numeric `oracle.equivalent`, so no value changes — only the surface form, now matching
  SymPy. `x²−2x+5 → 1 ± 2I`, `x⁴+4 → ±1 ± I`, `x³−8 → 2, −1 ± √3·I`.

### SUM-SHIFT-1 — infinite sums starting at an index ≠ 1 missed the closed-form handlers
- **Problem:** most closed-form `summation` handlers (arithmetic p-series, ζ, cotangent)
  key on `lo == 1`. So the *standard* odd p-series written from zero,
  `Σ_{n=0}^∞ 1/(2n+1)² = π²/8`, and any shifted-start variant
  (`Σ_{n=0}^∞ 1/(2n+1)⁴`, `Σ_{n=2}^∞ 1/(2n+1)²`, `Σ_{n=2}^∞ 1/n²`) returned an
  unevaluated `Sum(…)` even though the `lo=1` form evaluates fine.
- **Fix:** added an index-shift fallback in `src/calculus/summation.cpp`, tried only after
  every direct handler fails. For an infinite sum with an integer start `lo ≠ 1` it
  re-expresses `Σ_{n=lo}^∞ f(n) = Σ_{m=1}^∞ f(m + lo − 1)` (via `subs`) and recurses; the
  shifted call has `lo = 1` so it cannot loop. Its result is adopted only when it is a
  genuine closed form (var-free) — an unevaluated `Sum` still carries the bound variable
  and is rejected, so nothing that previously stayed symbolic changes. General over the
  summand (not just p-series): `Σ_{n=0}^∞ 1/(2n+1)² → π²/8`,
  `Σ_{n=2}^∞ 1/(2n+1)² → π²/8 − 10/9`. Matches SymPy.

### TOGETHER-NESTED-1 — `together`/`simplify` left compound (nested) fractions uncombined
- **Problem:** `together` decomposed only the top level via `as_numer_denom`, which (by
  design, for `integrate`'s sake) does not recurse. So a reciprocal of a sum of fractions
  stayed compound: `together(1/(1+1/x))` returned `(1/x+1)⁻¹` and `simplify(1/(1+1/x))`
  left it unchanged, where SymPy gives `x/(x+1)`. Nested and mixed-symbol forms
  (`1/(1+1/(1+1/x))`, `a/(b+c/d)`) were likewise stuck.
- **Fix:** made `together` recursive (`together_recursive` in `src/polys/poly.cpp`) without
  touching `as_numer_denom`. It combines each `Add`/`Mul`/`Pow` child into a single
  fraction bottom-up before recombining the current level; function arguments are left
  intact (shallow inside `sin`/`exp`/…, matching SymPy's default). The `Pow` case inverts a
  fractional base explicitly — `(bn/bd)^e = bn^e·bd^(−e)` — so `1/((x+1)/x)` flips to
  `x/(x+1)`. `1/(1+1/x) → x/(x+1)`, `1/(1+1/(1+1/x)) → (x+1)/(2x+1)`,
  `a/(b+c/d) → a·d/(b·d+c)`. Matches SymPy; plain sums of reciprocals are unchanged.

### LIMIT-DBG-1 — a debug `fprintf` leaked to stderr on every radical limit at +∞
- **Problem:** `try_algebraic_inf` in `src/calculus/limit.cpp` (the leading-asymptotic-term
  evaluator for `√`-difference limits at `+∞`) carried a leftover
  `std::fprintf(stderr, "DBG alginf …")` from its development. Every algebraic limit at
  infinity — `x − √(x²−x)`, `√(x²+x) − x`, … — printed a diagnostic line to stderr. The
  computed value was correct; only the stray output was wrong. SymPy emits nothing.
- **Fix:** removed the `fprintf` (and its now-unused transitive `<cstdio>` reliance). The
  handler is unchanged otherwise; the existing `LIMIT-RADICAL-INF-1` regression block plus
  a new `√(x²+2x) − x → 1` case keep the value path covered.

### EXPAND-TRIG-HYP-1 — `expand_trig(sinh(x+y))` left hyperbolic functions unexpanded
- **Problem:** `expand_trig` expanded the circular trio (sin/cos/tan angle-addition and
  multiple angles) but returned `sinh`/`cosh`/`tanh` of a sum or multiple angle
  untouched: `expand_trig(sinh(x+y))` stayed `sinh(x+y)` where SymPy gives
  `sinh(x)·cosh(y) + cosh(x)·sinh(y)`. The hyperbolic angle-addition identities were
  simply missing.
- **Fix:** extended `expand_trig_node` in `src/simplify/simplify.cpp` to dispatch on
  `Sinh`/`Cosh`/`Tanh` as well, reusing the existing Add / multiple-angle argument split
  and adding the Osborn-rule formulas: `sinh(a+b)=sinh a·cosh b + cosh a·sinh b`,
  `cosh(a+b)=cosh a·cosh b + sinh a·sinh b`, `tanh(a+b)=(tanh a+tanh b)/(1+tanh a·tanh b)`
  (note the `+1` denominator, vs `−1` for `tan`). Multiple angles reduce through the same
  `n·g = g + (n−1)·g` split and the `expand_trig` fixpoint: `sinh(2x)→2·sinh x·cosh x`,
  `cosh(3x)→3·sinh²x·cosh x + cosh³x`. Matches SymPy up to identity equivalence.

### TRIG-DBLRATIO-1 — `simplify(sin(2x)/sin(x))` did not reduce to `2·cos(x)`
- **Problem:** `simplify` collapsed the *product* `2·sin(x)·cos(x) → sin(2x)` but not the
  inverse *ratio*: `sin(2x)/sin(x)`, `sin(2x)/cos(x)`, and the `csc`/`sec` forms stayed
  unevaluated. SymPy returns `2·cos(x)` / `2·sin(x)`.
- **Fix:** added `trig_double_angle_ratio_mul` in `src/simplify/simplify.cpp`, run in the
  `trigsimp_node` Mul pipeline (after `trig_ratio_mul`, before `trigsimp_mul`). For a
  `sin(a)` factor it looks for a denominator `1/sin(u)`/`1/cos(u)` — written as a negative
  power or as `csc(u)`/`sec(u)` — with `a = 2u` (checked via `expand(a − 2u) == 0`), and
  folds the pair using `sin(2u) = 2·sin(u)·cos(u)`: `sin(2u)/sin(u) → 2·cos(u)`,
  `sin(2u)/cos(u) → 2·sin(u)`. Only the doubled sine factors into a single sin/cos pair,
  so `cos(2u)/cos(u)` and `sin(3u)/sin(u)` are left alone — exactly as SymPy does.
  Coefficients carry through (`3·sin(2x)/sin(x) → 6·cos(x)`) and deeper doubling cancels
  one level (`sin(4x)/sin(2x) → 2·cos(2x)`).

### SUM-INVQUAD-1 — `Σ_{n=1}^∞ 1/(n²+b)` was unevaluated (irreducible-quadratic denominator)
- **Problem:** convergent rational sums went through `apart`, which only splits
  denominators with rational roots. For an irreducible quadratic denominator with a
  positive constant — `Σ 1/(n²+1)`, `Σ 1/(n²+4)`, `Σ 1/(2n²+1)` — the poles are a
  complex-conjugate pair, so `apart` is a no-op and the sum stayed an unevaluated
  `Sum(…)`. SymPy returns the cotangent closed form (`-1/2 + π·coth(π)/2`, …).
- **Fix:** added `sum_inverse_quadratic` in `src/calculus/summation.cpp`, dispatched
  before the `apart` path. It peels a var-free coefficient `c` off a `(a·n²+b)^(-1)`
  factor, builds `Poly(denom, n)`, requires the linear term to vanish and
  `B = b/a > 0`, then returns the Mittag-Leffler / cotangent result
  `Σ_{n=1}^∞ 1/(n²+B) = (π·√B·coth(π·√B) − 1)/(2B)`, scaled by `c/a`. The `B > 0`
  guard keeps it off the `cot`/digamma cases (`n²−a²`) and off true p-series (`b=0`,
  still ζ). `Σ1/(n²+1)=(π·coth π−1)/2`, `Σ1/(n²+4)=−1/8+π·coth(2π)/4`. Matches SymPy.

### ILAPLACE-REPQUAD-1 — `iL{N(s)/(s²+a²)²}` was unevaluated (repeated irreducible quadratic)
- **Problem:** the inverse Laplace handled simple poles and the irreducible quadratic
  `(αs+β)/((s−a)²+b²)`, but a *repeated* irreducible quadratic denominator
  `(s²+a²)²` (which SymPP expands to a quartic, e.g. `s²/(8s²+s⁴+16)`) had no handler:
  `s/(s²+4)²`, `1/(s²+1)²`, `s/(s²+1)²`, `(s²−1)/(s²+1)²` all returned an unevaluated
  `InverseLaplaceTransform(…)`. This is the inverse of the LAPLACE-TMULT-1 rule
  (`L{t·sin/t·cos}` lands exactly on these), so the pair was asymmetric.
- **Fix:** added `inverse_laplace_repeated_quad` in `src/integrals/transforms.cpp`,
  dispatched before the generic `inverse_laplace_term`. It splits `F = N·D^(-1)` with `D`
  a quartic, builds `Poly(expand(D), s)`, requires the odd coefficients to vanish and the
  even ones to form a perfect square `(s²+a²)²` (`a²=p/2`, `a²²==q`, `a²>0`), then
  decomposes the degree-≤2 numerator over the three basis inverses
  `iL{1/(s²+a²)²}=(sin at − a·t·cos at)/(2a³)`, `iL{s/(s²+a²)²}=t·sin at/(2a)`,
  `iL{s²/(s²+a²)²}=sin at/(2a)+t·cos at/2`. `s/(s²+4)²→t·sin 2t/4`,
  `1/(s²+1)²→(sin t−t·cos t)/2`, `(s²−1)/(s²+1)²→t·cos t`. Matches SymPy.

### LAPLACE-TMULT-1 — `L{t·cos t}` was unevaluated (multiplication-by-tⁿ rule)
- **Problem:** the Laplace transform handled `tⁿ` and the s-shift `L{e^(at)·g}=G(s−a)`,
  so `t·e^t` worked, but `t·cos t`, `t·sin t`, `t²·cos t`, `t·sinh t` (a `t` factor times
  a trig/hyperbolic with no exponential) returned an unevaluated `LaplaceTransform(…)`.
- **Fix:** added the multiplication-by-`t` rule in `src/integrals/transforms.cpp`:
  `L{tⁿ·g(t)} = (−1)ⁿ·dⁿ/dsⁿ L{g(t)}`. In the Mul handler (no-exp path) it splits the
  positive integer powers of `t` from the rest `g`, transforms `g`, and differentiates
  its transform `n` times w.r.t. `s` (sign `(−1)ⁿ`). `L{t·cos t}=(s²−1)/(s²+1)²`,
  `L{t·sin t}=2s/(s²+1)²`, `L{t·sinh t}=2s/(s²−1)²`, `L{t²·cos t}=(2s³−6s)/(s²+1)³`. The
  exp cases still go through the s-shift; the two compose for `t·e^(at)·cos t`. Matches
  SymPy.

### LIMIT-NROOT-INF — `(x³+x²)^(1/3) − x → 1/3` (n-th-root difference) returned nan
- **Problem:** the leading-term conjugate in `leading_pos_inf` only handled *square*
  roots (`t₁+t₂ = (t₁²−t₂²)/(t₁−t₂)`), so cube/fourth-root differences whose leading
  terms cancel — `(x³+x²)^(1/3)−x`, `(x⁴+x³)^(1/4)−x`, `(x³+x²)^(1/3)−(x³−x²)^(1/3)` —
  stayed `nan` (the square conjugate leaves the cube root in the numerator). The
  reciprocal-substitution fallback also can't reach them (their substituted form
  doesn't resolve at `t→0`).
- **Fix:** generalized the conjugate in `src/calculus/limit.cpp` to the n-th root —
  `u − v = (uⁿ − vⁿ)/Σ_{i=0}^{n-1} u^(n−1−i)·vⁱ`, with `n` the LCM of the radical
  exponent denominators (new `radical_order` helper). `uⁿ`/`vⁿ` raise the radicals to
  integer powers, clearing them from the numerator; the denominator has no leading
  cancellation, so `leading_pos_inf` recurses cleanly. `n=2` is the original √ case.
  `(x³+x²)^(1/3)−x → 1/3`, `(8x³+x²)^(1/3)−2x → 1/12`,
  `(x³+x²)^(1/3)−(x³−x²)^(1/3) → 2/3`. Matches SymPy.

### LIMIT-RECIP-INF-1 — asymptotic limits at ∞ with a transcendental subleading term
- **Problem:** limits at +∞ whose value comes from a subleading asymptotic term —
  `x − x²·log(1+1/x) → 1/2`, `x²(1−cos(1/x)) → 1/2`, `x·(e^(1/x)−1) → 1`,
  `x²(e^(1/x)−1−1/x) → 1/2` — returned `nan`. The direct ∞ methods (polynomial leading
  term, conjugate, `try_algebraic_inf`, L'Hôpital) all abandon them.
- **Fix:** added a reciprocal-substitution fallback in `src/calculus/limit.cpp` (after
  L'Hôpital, for ±∞ targets only): substitute `x = ±1/t` to map the limit to `t → 0`,
  where the series/L'Hôpital machinery applies (e.g. `x²(1−cos 1/x) → (1−cos t)/t²`).
  The candidate is accepted only after a **numeric convergence check** — sampling the
  original at `x = 10³, 10⁶, 10⁹`, requiring the diff to not diverge and the largest
  sample to land within `1e-4` — so a one-sided/two-sided mismatch or a wrong `t`-limit
  cannot slip through (it leaves such cases `nan` rather than guessing). Also: L'Hôpital
  returning a `nan` value no longer short-circuits this fallback. (Algebraic n-th-root
  differences like `(x³+x²)^(1/3)−x` whose *substituted* form `(1/t³+1/t²)^(1/3)−1/t`
  still doesn't resolve at `t→0` remain a follow-up.) Matches SymPy.

### ASSUME-MULSIGN-1 / ASSUME-REALFINITE-1 — last assumption-propagation gaps vs SymPy
- **Problem:** a 39-query SymPP-vs-SymPy assumption battery left two genuine gaps
  (the others were SymPP being *more* correct — `Abs(x)` is always real/nonnegative,
  where SymPy returns unknown): (1) `Mul` did not propagate non-strict signs, so
  `positive·nonnegative` → `is_nonnegative` was unknown; (2) `real` did not imply
  `finite`, so `exp(r)` for real `r` was not known finite (SymPy: `real ⇒ finite`).
- **Fix:** `Mul::ask` (`src/core/mul.cpp`) now decides Nonnegative/Nonpositive by the
  parity of the ≤0 factors when every factor has a known sign direction (≥0 or ≤0); a
  provably-zero factor makes the product 0. Conservative — only the definite direction
  is reported. And `close_assumptions` (`src/core/assumption_mask.cpp`) adds `real ⇒
  finite` (the unbounded ±∞ are the separate Infinity atoms; consistent with the
  existing positive/negative/zero ⇒ finite). The battery now matches or exceeds SymPy
  on all 39 queries.

### SIMPLIFY-LOGRATIO-1 — `simplify(log(4)/log(2))` stayed unevaluated (should be 2)
- **Problem:** `log(b)/log(a)` for integer `a, b` that are powers of a common base —
  `log(4)/log(2)`, `log(8)/log(2)`, `log(2)/log(8)` — was left as
  `log(2)⁻¹·log(4)` rather than the rational `log_a(b)` (2, 3, 1/3). SymPy's `simplify`
  reduces these. It also left exponential-equation roots (SOLVE-EXPBASE-SUM-2) as
  `log(4)/log(2)` instead of `2`.
- **Fix:** added a `log_ratio` pass in `src/simplify/simplify.cpp`. On a `Mul`
  carrying a `log(b)` factor and a `log(a)⁻¹` factor (`a, b` positive integers ≥ 2),
  it takes the primitive base `c` of each (smallest `c` with `n = cᵏ`); when both share
  `c` (`b = cʲ`, `a = cⁱ`) it replaces the pair with `j/i`. Incommensurate args
  (`log(2)/log(3)`) and non-power args (`log(6)/log(2)`) are left alone; other factors
  pass through (`x·log(8)/log(2) → 3x`). As a bonus the exponential-quadratic roots now
  render exactly: `4ˣ−5·2ˣ+4 → {0, 2}`, `16ˣ−6·4ˣ+8 → {1/2, 1}`.

### SOLVE-EXPBASE-SUM-2 — `solve(4ˣ−2ˣ−2)` returned `[]` (composite exponential base)
- **Problem:** an exponential quadratic written with a *composite* base — `4ˣ−2ˣ−2=0`,
  `9ˣ−4·3ˣ+3=0` — returned `[]`, even though `2^(2x)−5·2ˣ+4` (same base, scaled exponent)
  solved. The `u = 2ˣ` substitution requires every base to be a power of a common base;
  with bases `{4, 2}` the rate selection picked `log 4` as the unit, giving the
  non-integer ratio `log 2 / log 4 = 1/2`, so the commensurate-rate path bailed.
  (`simplify` doesn't reduce `log 4 / log 2 → 2`, so the rate comparison can't recover.)
- **Fix:** `solve_const_base_exp_sum` (`src/solvers/solve.cpp`) now normalizes the bases
  first: if every constant integer base is a power of the smallest one `b`, it rewrites
  `aˣ = b^(k·x)` (`4ˣ → 2^(2x)`) and re-solves. The existing `u = bˣ` machinery then
  closes it. `4ˣ−2ˣ−2 → 1`, `9ˣ−4·3ˣ+3 → {0,1}`, `16ˣ−6·4ˣ+8 → {1/2,1}`, matching
  SymPy's real roots. (Some roots render as `log 4 / log 2`, value-correct but
  unsimplified — the log-ratio reduction is a separate cosmetic gap.)

### ASSUME-NONNEG-1 — `nonnegative`/`nonpositive` could not be declared (silently lost)
- **Problem:** `AssumptionMask` stored only the *primary* sign facts (positive, negative,
  zero); `nonnegative`/`nonpositive` were derived (positive∨zero / negative∨zero) but had
  no field, so `set(Nonnegative, true)` only recorded `negative=false` — which doesn't
  reconstruct nonnegativity. A symbol declared nonnegative reported `is_nonnegative →
  None`, `is_real → None`, and `√(x²)` would not simplify to `x`. SymPy expresses this as
  `Symbol('x', nonnegative=True)`.
- **Fix:** added `nonnegative`/`nonpositive` as stored fields on `AssumptionMask`
  (`include/sympp/core/assumption_mask.hpp`) with `set_nonnegative`/`set_nonpositive`/
  `set_nonzero` builders, included in `hash()`/`empty()`. `close_assumptions`
  (`src/core/assumption_mask.cpp`) gained the rules: `nonnegative ⇒ real ∧ finite ∧
  ¬negative`, refining to `positive` (when `≠0`) or `zero` (when `¬positive`); the
  primaries imply them back (`positive ⇒ nonnegative`, `negative ⇒ nonpositive`); the
  strict signs exclude the opposite; and `¬real ⇒ ¬nonnegative ∧ ¬nonpositive`. The
  change is additive — existing masks leave the new fields `nullopt`, so behavior is
  unchanged. Declared nonnegativity now flows into simplification (`√(x²)→x`, `|x|→x`)
  and the MATLAB `assume(x,"nonnegative")` facade. Matches SymPy.

### TRIG-ANGLE-ADD-2 — N-term angle addition; `dsolve(y''+4y=sin x)` had a messy particular solution
- **Problem:** the angle-addition simplifier `sin(a)cos(b)±cos(a)sin(b)→sin(a±b)` only
  ran on a *two-term* Add. So a longer trig-product combination — e.g. the
  variation-of-parameters particular solution of `y''+4y=sin(x)` — was left as a large
  unsimplified mess rather than `sin(x)/3` (SymPy gives `C1 sin2x + C2 cos2x + sin(x)/3`).
- **Fix:** in `src/simplify/simplify.cpp`, generalized `trigsimp_angle_addition` to an
  Add of any size: it greedily collapses every collapsible pair of two-trig products
  and lets `add()` collect the folded terms (so `1/12 sin x + 1/4 sin x → sin(x)/3`).
  Also guarded `as_cos_double_term` so it no longer grabs `cos(2x)` out of a genuine
  product `sin(x)·cos(2x)` (treating `sin(x)` as the coefficient and folding it into
  `cos²−sin²`), which was mangling the Add before the angle-addition pass could run.
  `dsolve` already simplifies its result, so the higher-order nonhomogeneous trig
  solutions now come out clean (matching SymPy).

### SIMPLIFY-RADCOEFF-1 — `simplify(√(4a²))` didn't pull out the perfect square
- **Problem:** SymPP factored perfect powers out of *pure-number* radicals (`√8 → 2√2`,
  auto-evaluated) but not when a symbolic factor was present: `√(4a²)`, `√(8x²)`,
  `(8x³)^(1/3)` stayed as `(4a²)^(1/2)` etc., where SymPy gives `2√(a²)`, `2√2·√(x²)`,
  `2·x^(1/3)`.
- **Fix:** added a `radical_coeff` pass in `src/simplify/simplify.cpp`. For a `Pow`
  with a non-integer rational exponent over a `Mul` base with a *positive* numeric
  coefficient `c`, it pulls out the perfect-power part of `c` (via `c^exp`, which
  auto-factors), keeping the non-perfect remainder under the radical with the symbolic
  factors: `√(8x²) → 2√(2x²)`. Valid because `c > 0` makes `(c·X)^e = c^e·X^e` hold on
  the principal branch regardless of `X`'s sign. It runs *after* the anti-bloat guard
  (the extraction can raise the node count yet is a genuine simplification, so the
  guard must not revert it). Equivalent to SymPy (up to the `√c·√X ↔ √(c·X)` regroup).

### SUM-EXP-NOLEAK — `Σcos(k·x)/k!` returned a bogus `e·cos(k·x)` (bound-variable leak)
- **Problem:** `sum_exponential_series` built `Poly(numerator, k)` without checking the
  resulting coefficients are var-free. `Poly()` treats a non-polynomial factor
  (`cos(k)`, `cos(k·x)`, `√k`, …) as an opaque degree-0 coefficient, so the handler
  pulled it out and returned a *wrong* closed form containing the summation variable:
  `Σ cos(k·x)/k! → e·cos(k·x)`, `Σ cos(k)/k! → e·cos(k)`. A summation result must never
  contain the bound variable.
- **Fix:** after building the coefficient list in `src/calculus/summation.cpp`, reject
  any coefficient that still depends on `k` (`for (cf : pcoeffs) if (has(cf, var))
  return nullopt`). The sum now stays unevaluated — matching SymPy, which also returns
  an unevaluated `Sum` here — instead of a wrong answer. Legitimate polynomial
  numerators (`Σ k/k! = e`, `Σ k²/k! = 2e`, `Σ k·xᵏ/k! = x·eˣ`) are unaffected.

### SUM-RATIONAL-1 — `Σ1/(k²(k+1)) = π²/6 − 1` (general rational sum) was unevaluated
- **Problem:** a convergent rational sum mixing a ζ part and a telescoping part —
  `Σ1/(k²(k+1)) = π²/6 − 1`, `Σ1/(k(k+1)²) = 2 − π²/6`, `Σ1/(k²(k+2)) = π²/12 − 3/8` —
  stayed unevaluated. The 2-term-apart telescoping (SUM-TELESCOPE-3) only fires when
  every partial-fraction term cancels into a single `g(k)−g(k+1)`.
- **Fix:** added `sum_rational_via_apart` in `src/calculus/summation.cpp`. It `apart()`s
  the summand and groups the terms: each pole of order `j ≥ 2` sums on its own (the
  arithmetic p-series path gives the ζ-value, e.g. `Σ1/k² = π²/6`), while the simple
  poles (`j = 1`) are recombined into one fraction and handed to `telescope_rational`
  (their residues sum to zero for a convergent series, so they telescope). The two
  parts are added. The recombined `j = 1` fraction is `simplify`'d first so its
  numerator collapses to the var-free constant `telescope_rational` requires. Infinite
  range only (a finite `j ≥ 2` part would need harmonic numbers). Matches SymPy.

### SUM-TELESCOPE-3 — `Σ(2k+1)/(k²(k+1)²)=1` (repeated-root telescoping) was unevaluated
- **Problem:** rational summands that telescope only after partial fractions —
  `(2k+1)/(k²(k+1)²) = 1/k² − 1/(k+1)²`, `(3k²+3k+1)/(k³(k+1)³) = 1/k³ − 1/(k+1)³` —
  were unevaluated. The explicit-difference check sees `g(k)−g(k+1)` only when the
  summand is already written that way, and `telescope_rational` skips repeated roots.
- **Fix:** in `src/calculus/summation.cpp`, before the expand-and-recurse, `apart()`
  the rational summand; if it becomes a 2-term `g(k) − g(k+1)`, return the telescoping
  closed form `g(lo) − g(hi+1)`. A pole guard (no integer root of the denominator
  ≥ `lo`, plus a finite-endpoint check) prevents a bogus singular value when a pole
  lies in the range (`Σ_{k=1} 1/(k(k−1))` stays unevaluated; `Σ_{k=2}` now closes to 1).
  A residual ζ part (`1/(k²(k+1)) → −1/k+1/(k+1)+1/k²`, a 3-term apart) falls through.

### SOLVE-RADISOLATE-2 — `solve(√x + √(x+1) = 3)` returned `[]`
- **Problem:** the isolate-and-square radical solver handled exactly **one** square
  root, so equations with two — `√x + √(x+1) = 3` (→ 16/9), `√(2x+1) − √x = 1`
  (→ 0, 4), `√(x−1) + √(x+4) = 5` (→ 5) — returned `[]`.
- **Fix:** `solve_radical_isolate` (`src/solvers/solve.cpp`) now also accepts two
  radicals. Writing `expr = A1·√g1 + A2·√g2 + P` (A1, A2, P radical-free), it isolates
  and squares once — `A1²·g1 = A2²·g2 + 2·A2·P·√g2 + P²` — leaving a single radical
  that the same path (size 1) then clears. Candidates are filtered against the
  *original* equation (numeric back-substitution) to drop the roots squaring
  introduces. A √g1·√g2 cross term or radical-dependent coefficient falls through.
  Matches SymPy.

### SUM-COSH-SINH-1 — `Σ x^(2k)/(2k)! = cosh x` (even/odd factorial series) was unevaluated
- **Problem:** the exponential-series handler matched only a `k!` denominator, so the
  even/odd bisection of the exponential series — `Σ z^(2k)/(2k)! = cosh z`,
  `Σ z^(2k+1)/(2k+1)! = sinh z`, and the `(−1)^k`-signed `cos z`/`sin z` — stayed
  unevaluated (`Σ1/(2k)!`, `Σx^(2k)/(2k)!`, …).
- **Fix:** added `sum_cosh_sinh_series` in `src/calculus/summation.cpp`. It recognizes
  a `(2k+b)!` denominator (`b ∈ {0,1}`), an optional `(−1)^k` sign, and a numerator
  `z^(2k+b)` (constant numerator → `z = 1`); the result is `cosh z`/`sinh z` (no sign,
  by `b`) or `cos z`/`sin z` (with sign). A lower bound `lo > 0` subtracts the finite
  head `Σ_{k=0}^{lo−1}`. Dispatched before the expand-and-recurse. Matches SymPy.

### SUM-EXP-SHIFT-1 — `Σ1/(k+1)! = e−2` and shifted-factorial e-sums were unevaluated
- **Problem:** the exponential-series handler matched only a bare `k!` denominator, so
  the e-valued sums over a *shifted* factorial — `Σ1/(k+1)!=e−2`, `Σ(2k+1)/(k+1)!=e`,
  `Σk/(k+2)!=3−e`, `Σ1/(k+2)!=e−2` — stayed unevaluated. (These are the non-telescoping
  companions of SUM-FACT-TELESCOPE-1.)
- **Fix:** `sum_exponential_series` (`src/calculus/summation.cpp`) now re-indexes a
  shifted factorial: `(k+m)!` with `j=k+m` turns `Σ_{k=lo}^∞ P(k)/(k+m)!` into
  `Σ_{j=lo+m}^∞ P(j−m)/j!`, the `m=0` case it already closes (`Q(1)·e` via the
  falling-factorial transform, minus the omitted head). Implemented as a `subs(var →
  var−m)` with the lower bound shifted to `lo+m`; the recursion bottoms out at the bare
  `factorial(var)`. A non-unit var coefficient (`(2k)!` → `cosh 1`) is left alone.
  Matches SymPy.

### SUM-FACT-TELESCOPE-1 — `Σ k/(k+1)! = 1` (factorial telescoping) was unevaluated
- **Problem:** sums like `Σ_{k=1}^∞ k/(k+1)! = 1` and `Σ (k²−1)/(k+1)! = 1` were left as
  a partially-split unevaluated `Sum`. The exponential-series handler only matches a
  `k!` denominator (`Σ P(k)/k!`), and the generic expand-and-recurse splits the
  numerator — destroying the telescoping structure (these telescope as a whole, not
  term by term).
- **Fix:** added `sum_factorial_telescope` in `src/calculus/summation.cpp` — Gosper's
  algorithm specialized to a factorial denominator `(k+m)!`. The antidifference, if it
  exists, is `g(k) = Q(k)/(k+m−1)!` with `P(k)/(k+m)! = g(k) − g(k+1)`; multiplying by
  `(k+m)!` gives the polynomial identity `Q(k)·(k+m) − Q(k+1) = P(k)`, solved top-down
  for `Q` of degree `deg(P)−1`. The constant-term equation is a consistency check that
  fails for non-telescoping terms (`1/(k+1)! → e−2` is left unevaluated, correctly).
  The sum is then `g(lo) − g(hi+1)` (`g(∞)=0`). Dispatched before the expand-and-recurse
  so the numerator stays intact. Handles infinite and finite ranges
  (`Σ_{k=1}^n k/(k+1)! = 1 − 1/(n+1)!`). Matches SymPy.

### SOLVE-LAMBERT-2 — `solve(eˣ = x + 2)` returned `[]`
- **Problem:** the additive Lambert-W solver required the bare-`var` term to have a
  unit coefficient (`t == var`), so `eˣ − x − 2 = 0` (coefficient −1 on `x`) fell
  through and `solve` returned `[]`, even though SymPy gives `−2 − W(−e⁻²)`. The same
  blocked `eˣ − 2x − 1`, `2x + eˣ`, `2x + log(x)`, …
- **Fix:** generalized the additive branch of `solve_lambert`
  (`src/solvers/solve.cpp`) to a free coefficient `a` on the var term (recovered as
  `t/var`): `a·var + eᵛᵃʳ + c = 0 → var = −W(e^(−c/a)/a) − c/a`, and the log analogue
  `a·var + log(var) + c = 0 → var = W(a·e^(−c))/a`. The unit-coefficient cases are the
  `a = 1` special case, so existing results (`x + eˣ − 1 → 0`, `x + log(x) − 1 → 1`)
  are unchanged. Matches SymPy across `eˣ − x − 2`, `2x + eˣ`, `2x + log(x)`, etc.

### SUM-ARITH-PSERIES-1 — `Σ1/(2k−1)² = π²/8` and arithmetic p-series were unevaluated
- **Problem:** the p-series handler only recognized `1/kˢ` (base exactly the index),
  so the classic Basel relatives `Σ1/(2k−1)²=π²/8`, `Σ1/(2k)²=π²/24`, `Σ1/(2k−1)⁴=π⁴/96`
  stayed unevaluated even though `ζ(even)` was already known.
- **Fix:** added an arithmetic-argument p-series handler in
  `src/calculus/summation.cpp` for `Σ_{k=1}^∞ c/(a·k+b)ˢ`, `s ≥ 2` integer, `a ∈ {1,2}`.
  The denominator runs over one residue class, so the value is the matching slice of
  `ζ(s)` minus the finitely many leading terms: `a=1,b≥0` → `ζ(s) − Σ_{n=1}^{b} n⁻ˢ`;
  `a=2` odd `b` → `(1−2⁻ˢ)ζ(s) − Σ(2j−1)⁻ˢ` (odd n); `a=2` even `b` → `2⁻ˢζ(s) − Σ(2j)⁻ˢ`
  (even n). `ζ(even)` closes to a `πˢ` rational; odd `s` stays a symbolic `ζ(s)`
  (`Σ1/(2k−1)³ = 7ζ(3)/8`), as SymPy does. `a ≥ 3` needs Hurwitz `ζ` and falls through.

### SUM-TELESCOPE-2 — `Σ1/(k(k+1)(k+2))` (degree ≥ 3 telescoping) was unevaluated
- **Problem:** the rational telescoping handler only covered a *quadratic*
  denominator, so `Σ1/(k(k+1)(k+2)) = 1/4`, `Σ1/(k(k+1)(k+2)(k+3)) = 1/18` and
  `Σ1/((2k−1)(2k+1)(2k+3)) = 1/12` stayed unevaluated even though the 2-factor cases
  (`Σ1/(k(k+1))`, `Σ1/(4k²−1)`) worked.
- **Fix:** generalized `telescope_rational` (`src/calculus/summation.cpp`) to any
  denominator of degree ≥ 2 whose roots are rational and pairwise differ by integers.
  Partial fractions give `c/D = Σ Aᵢ/(k−rᵢ)` with `Aᵢ = c/(lead·∏_{j≠i}(rᵢ−rⱼ))`;
  taking the largest root as a reference, each `1/(k−rᵢ) = u(k+mᵢ)` (`mᵢ = ref−rᵢ ≥ 0`),
  so the summand is `Σ Aᵢ(u(k+mᵢ)−u(k))` (the `−u(k)` parts cancel since `ΣAᵢ = 0` for a
  constant numerator over degree ≥ 2). Each piece telescopes to
  `Σ Aᵢ[Σ_{j=1}^{mᵢ}u(hi+j) − Σ_{j=0}^{mᵢ−1}u(lo+j)]`, exact for finite or infinite
  `hi`. The pole guard (no integer root ≥ `lo`) and var-free-numerator restriction are
  retained; non-integer root gaps (which need digamma) safely fall through.

### LIMIT-RADICAL-INF-1 — `lim √(x²+x)−x` (nonzero) returned `nan`
- **Problem:** √-difference limits at +∞ with a *nonzero* finite value returned `nan`
  (a wrong answer): `√(x²+x)−x → 1/2`, `x−√(x²−x) → 1/2`, `√(x²+x)−√(x²−x) → 1`,
  `x·(√(x²+1)−x) → 1/2`. The conjugate handler clears the ∞−∞ but leaves a residual
  ∞/∞ ratio (e.g. `x/(√(x²+x)+x)`) that L'Hôpital abandons on radicals — repeated
  differentiation balloons the nested radical and never stabilises. (The zero-valued
  cases like `√(x²+1)−x → 0` already worked, because there the conjugate numerator is
  constant, giving const/∞ = 0 with no ∞/∞.)
- **Fix:** added a leading-asymptotic-term evaluator `leading_pos_inf` (the leading
  slice of Gruntz/MRV restricted to polynomials and their roots) plus a
  `try_algebraic_inf` handler in `src/calculus/limit.cpp`, dispatched in the nan/+∞
  branch after the conjugate. It returns `e ~ c·x^d` (degree may be rational, since √
  halves it); the limit is `c` when `d=0`, `±∞` when `d>0`, `0` when `d<0`. On a
  leading cancellation it applies the conjugate identity `t₁+t₂=(t₁²−t₂²)/(t₁−t₂)` and
  recurses, so it also handles the 0·∞ product `x·(√(x²+1)−x)`. Restricted to +∞ (the
  evaluator assumes `x>0` to pull `x` out of a radical); −∞ remains a follow-up.

### SUM-POLYGEOM-SYM-1 — `Σ_{k=1}^n k·xᵏ` (symbolic ratio) was unevaluated
- **Problem:** the polynomial × geometric closed form `Σ P(k)·rᵏ` was gated to a
  *numeric* ratio (`Σk·2ᵏ` worked), so the generating-function identity
  `Σ_{k=1}^n k·xᵏ = x(1−(n+1)xⁿ+n·xⁿ⁺¹)/(x−1)²` — and `Σk²·xᵏ`, `Σk·aᵏ`, … — stayed
  unevaluated for a symbolic base.
- **Fix:** `sum_poly_geometric` (`src/calculus/summation.cpp`) no longer requires the
  geometric base/ratio to be a number — only that the base is var-free and the ratio
  ≠ 1. The antidifference recurrence and finite boundary evaluation work unchanged
  symbolically. A finite sum now yields the clean closed form (matching SymPy's
  general branch; like finite geometric, no `x=1` Piecewise is emitted). An infinite
  sum with a symbolic ratio fails the `|r| < 1` convergence test and is left
  unevaluated rather than emitting `x**∞` terms — consistent with the existing
  numeric-ratio convergence handling.

### SOLVE-ROOTOF-1 — `solve(x⁵−x−1)` returned `[]` (claiming "no solutions")
- **Problem:** an irreducible polynomial of degree ≥5 is not solvable by radicals,
  so the closed-form solver (Cardano/Ferrari for ≤4, rational roots above) produced
  nothing and `solve()` returned an empty list — implying *no solutions* for, e.g.,
  `x⁵−x−1` (which has a real root ≈1.1673) or `2x⁵−10x+5` (three real roots). An
  empty list is a silently wrong answer, worse than an unevaluated result.
- **Fix:** `solve_poly` (`src/solvers/solve.cpp`) now supplements the radical roots:
  when the polynomial is degree 5..12 and roots are missing, it `factor_list`s and,
  for each irreducible factor of degree ≥5, emits `RootOf(factor, var, k)` (rendered
  `CRootOf(poly, k)`, matching SymPy) for each real root — detected by `try_evalf`
  returning a value (it yields `nullopt` past the last real root). The degree window
  avoids paying for (exponential Kronecker) factorization on the common low-degree
  path and bounds worst-case cost. The `solve()` post-filter that drops var-dependent
  candidates now exempts `RootOf`, which legitimately embeds its defining polynomial.
- **Known limitation (partial parity):** SymPP's `RootOf` is **real-root-only**, so
  the *complex* roots of these factors are not yet returned (SymPy returns all via
  `CRootOf`). A polynomial with only complex roots (e.g. `x⁶+x+1`) still yields `[]`.
  Complex-root isolation is the planned follow-up (SOLVE-ROOTOF-2).

### INT-TRIGPROD-1 — `∫sin²(x)cos(2x)` and trig products of mixed arguments were unevaluated
- **Problem:** products of sin/cos powers whose arguments are not all equal —
  `∫sin²(x)cos(2x)`, `∫cos²(x)cos(2x)`, `∫sin²(x)sin(2x)`, `∫sin³(x)cos(2x)`,
  `∫sin²(2x)cos(x)` — were unevaluated. The single-factor product-to-sum
  (`∫sin(x)cos(2x)`) worked, but `try_trig_reduction`'s Mul/half-angle branch
  deliberately defers any trig×trig product, and `try_trig_power` only handles a
  *same-argument* `sinᵐ·cosⁿ`, so mixed-argument products fell through.
- **Fix:** added `try_trig_product_expand` in `src/integrals/integrate.cpp`
  (dispatched after `try_sin_cos_quotient`). Any product of `sin/cos(affine)^k`
  linearizes — by repeated product-to-sum and power reduction — into a sum of single
  `sin(affine)`/`cos(affine)` terms (closed under ±/×2 of affine arguments), each of
  which the table integrates. Gated to ≥2 factors with at least two *distinct*
  arguments so same-argument products still go to `try_trig_power`. A high-precision
  numeric diff-back guards the result. Note SymPy's own `simplify` can't reliably
  reduce a trig product (`sin³(x)cos(2x)` etc.), so the regression test verifies by
  numeric sampling via the oracle's `evalf_is_zero` rather than `equivalent`.

### INT-SINCOS-QUOT-1 — `∫cos²/sin`, `∫sin³/cos²` and sin/cos quotients were unevaluated
- **Problem:** sin/cos quotients such as `∫cos²/sin`, `∫sin²/cos`, `∫sin³/cos`,
  `∫cos³/sin`, `∫sin³/cos²`, `∫cos²/sin³` all returned unevaluated. `try_trig_power`
  had the right `u = sin`/`u = cos` substitution machinery but `parse_sin_cos_powers`
  only accepted non-negative exponents, so every quotient fell through.
- **Fix:** `parse_sin_cos_powers` now accepts negative integer exponents (the degree
  guard uses `|m|+|n|`), and a new `try_sin_cos_quotient` in
  `src/integrals/integrate.cpp` (dispatched after `try_tan_sec_product`) handles the
  quotient case. When at least one exponent is odd, substituting `u = sin(g)` (cos
  odd) or `u = cos(g)` (sin odd) turns the integrand into a *rational* function of
  `u`, which `try_rational` closes. A numeric diff-back self-check gates the result,
  so a mis-step fails to a marker rather than a wrong answer. `try_trig_power` keeps
  its positive-only path (an early `m<0||n<0` guard) so existing sec/csc/tan handlers
  are not shadowed. (Both-even quotients such as `cos⁴/sin²` are handled by
  INT-SINCOS-QUOT-2 below.)

### INT-SINCOS-QUOT-2 — even/even sin/cos quotients (`∫cos⁴/sin²`, `∫cos²/sin²`) were unevaluated
- **Problem:** sin/cos quotients with *both* exponents even — `∫cos⁴/sin²`,
  `∫cos²/sin²` (=`∫cot²`), `∫sin⁴/cos²`, `∫cos²/sin⁴`, … — were unevaluated. The odd
  case (INT-SINCOS-QUOT-1) substitutes `u=sin`/`u=cos`, but with both powers even that
  leaves a `√` and never rationalizes.
- **Fix:** extended `try_sin_cos_quotient` with a both-even branch using `t=tan(g)`:
  `sinᵐcosⁿ dx = (1/a)·tᵐ/(1+t²)^((m+n)/2+1) dt`, which is rational in `t` for even
  `m,n`, so `try_rational` closes it. The same numeric diff-back gate applies; it now
  also accepts an exactly-zero residual (`simplify(diff−integrand)==0`), which had been
  mis-counted as "unverifiable" and wrongly rejected some correct antiderivatives.
  Results carry an `atan(tan(x))` term (a valid antiderivative; SymPy renders it `x`).

### INT-TANSEC-1 — `∫tan³(x)·sec(x)` and tan^m·sec^n products were unevaluated
- **Problem:** `∫tan³·sec`, `∫tan²·sec`, `∫tan³·sec³` and the cot/csc analogues
  were unevaluated. `∫tan·sec³` worked (heurisch with `u = sec`), but higher tan
  powers need `tan² = sec²−1`, which heurisch doesn't apply.
- **Fix:** added `try_tan_sec_product` in `src/integrals/integrate.cpp` (dispatched
  after the pure sec/csc-power handler). For `tan(g)^m·sec(g)^n` (g affine): when
  `m` is odd, `u = sec(g)` turns it into a polynomial `(u²−1)^((m−1)/2)·u^(n−1)`;
  when `m` is even, `tan^m = (sec²−1)^(m/2)` is expanded to pure sec powers and
  reduced via `try_sec_csc_power`. The cot/csc analogue carries the `d(csc) =
  −csc·cot` sign.
- **Verified:** `∫tan³·sec = sec³/3 − sec`, `∫tan²·sec`, `∫tan³·sec³`, `∫tan²·sec²`,
  `∫cot³·csc`, `∫cot²·csc` — all diff-back to the integrand, matching SymPy; the
  existing `∫tan·sec³` is unchanged.
- **Regression test:** `INT-TANSEC-1` in `tests/integrals/integrate_test.cpp`.

### INT-PROD2SUM-1 — `∫sin(2x)·sin(3x)` and sin·sin / cos·cos products were unevaluated
- **Problem:** `∫sin(2x)·sin(3x)`, `∫cos(2x)·cos(5x)`, `∫cos x·cos 2x·cos 3x`,
  `∫x·sin 2x·cos 3x` were unevaluated. The product-to-sum block in
  `try_trig_reduction` only handled the `sin p·cos q` pairing, not `sin·sin` or
  `cos·cos`, and only a single pair.
- **Fix:** generalized the block to collapse the first two sin/cos factors of
  distinct var-dependent arguments via the matching identity (`sin p·sin q =
  (cos(p−q) − cos(p+q))/2`, `cos p·cos q = (cos(p−q) + cos(p+q))/2`,
  `sin p·cos q = (sin(p+q) + sin(p−q))/2`), then `expand` and recurse — so a
  three-way product reduces pair by pair and a polynomial factor distributes for
  per-term integration.
- **Verified:** `∫sin 2x·sin 3x = sin x/2 − sin 5x/10`, `∫cos 2x·cos 5x`,
  `∫cos x·cos 2x·cos 3x`, `∫sin x·sin 2x·sin 3x`, `∫x·sin 2x·cos 3x` — all
  diff-back to the integrand, matching SymPy; the existing `sin·cos` case unchanged.
- **Regression test:** `INT-PROD2SUM-1` in `tests/integrals/integrate_test.cpp`.

### SUM-TELESCOPE-DIFF-1 — `Σ(1/k − 1/(k+1))` (explicit telescoping difference) was unevaluated
- **Problem:** an explicit telescoping difference `Σ(g(k) − g(k+1))` was not
  recognized: `Σ(1/k − 1/(k+1))`, `Σ(1/k! − 1/(k+1)!)`, `Σ(1/k² − 1/(k+1)²)` all
  returned an unevaluated `Sum`. Linearity split the Add into two sums, neither of
  which has a closed form (harmonic / factorial), so the telescoping was lost. (The
  existing telescoping handler only sees the *combined* rational form like
  `1/(k(k+1))`.)
- **Fix:** in `src/calculus/summation.cpp`, a 2-term Add is checked for the pattern
  `g(k) − g(k+1)` (via `t1 + g(k+1) == 0`) *before* the linearity split, returning
  `g(lo) − g(hi+1)`. This also closes factorial differences the rational
  partial-fraction path can't.
- **Verified:** `Σ(1/k − 1/(k+1)) = 1 − 1/(n+1)`, `Σ(1/k! − 1/(k+1)!) = 1 − 1/(n+1)!`,
  `Σ(1/k² − 1/(k+1)²) = 1 − 1/(n+1)²`, matching SymPy; a non-telescoping difference
  (`1/2^k − 1/3^k`, both geometric) is unaffected (falls through to the geometric
  handlers, → 1/2).
- **Regression test:** extended `SUM-TELESCOPE-1` in
  `tests/calculus/series_limit_test.cpp`.

### SUM-BINOMIAL-1 — `Σ_{k=0}^n C(n,k)` (binomial theorem) stayed unevaluated
- **Problem:** binomial-theorem sums `Σ_{k=0}^n C(n,k)·rᵏ = (1+r)ⁿ` were unevaluated:
  `Σ C(n,k) = 2ⁿ`, `Σ(−1)ᵏC(n,k) = 0`, `Σ2ᵏC(n,k) = 3ⁿ`, `ΣxᵏC(n,k) = (1+x)ⁿ`, and
  even the concrete `Σ_{k=0}^5 C(5,k) = 32`.
- **Fix:** added `sum_binomial_theorem` in `src/calculus/summation.cpp`. For a
  summand `const·binomial(n,k)·base^(a·k+b)` over `k = 0…n` — where `n` is exactly
  the binomial's first argument and the geometric factor is optional — it returns
  `const·base^b·(1 + base^a)ⁿ`, with `(1−1)ⁿ = 0` for the alternating case.
- **Verified:** `Σ C(n,k) = 2ⁿ`, `Σ(−1)ᵏC(n,k) = 0`, `Σ2ᵏC(n,k) = 3ⁿ` (which SymPy
  itself leaves unevaluated), `ΣxᵏC(n,k) = (1+x)ⁿ`, `Σ_{k=0}^5 C(5,k) = 32`; a
  mismatched argument `Σ C(m,k)` over `k=0…n` is correctly left unevaluated.
- **Regression test:** `SUM-BINOMIAL-1` in `tests/calculus/series_limit_test.cpp`.

### LIMIT-ESSENTIAL-PT-1 — `lim_{x→0} x/(exp(1/x)−1)` returned `nan`
- **Problem:** limits at a finite point with an *essential* singularity —
  `exp(−1/x²) → 0`, `x/(exp(1/x)−1) → 0`, `x²/(exp(1/x²)−1) → 0` — returned `nan`.
  Direct substitution evaluates `exp(1/x)` at `x = 0` as `exp(zoo) = nan`, and no
  method recovered.
- **Fix:** added a reciprocal-substitution fallback in `src/calculus/limit.cpp`:
  at a finite target `a` whose direct value is non-finite and which carries a
  reciprocal singularity (a negative power of a factor vanishing at `a`),
  substitute `u = 1/(x − a)` and take the `u → +∞` and `u → −∞` one-sided limits.
  They agree iff the two-sided limit exists, so the result is returned only then —
  genuinely two-sided-divergent cases (`exp(1/x)`, `1/x`) keep their `nan`/`zoo`.
- **Verified:** `exp(−1/x²) → 0`, `x/(exp(1/x)−1) → 0`, `x²/(exp(1/x²)−1) → 0`,
  matching SymPy; `exp(1/x)` stays `nan` and `1/x` stays `zoo` (two-sided DNE),
  and ordinary pole limits (`1/x² → ∞`, `1/(x−1)² → ∞`) are unchanged.
- **Regression test:** `LIMIT-ESSENTIAL-PT-1` in
  `tests/calculus/series_limit_test.cpp`.

### INT-INVTRIG-SQRT-SQ-1 — `∫asin(x)²` and √-derivative inverse-function squares were unevaluated
- **Problem:** `∫asin(x)²` (= `x·asin² − 2x + 2√(1−x²)·asin`), `∫x·asin(x)²`,
  `∫acos²`, `∫asinh²`, `∫asin³` were unevaluated, though elementary. (An earlier
  attempt returned *wrong* answers — blocked by the `try_x_over_sqrt_quadratic`
  coefficient bug, fixed in `INT-XSQRTQUAD-NUM-1`.)
- **Fix:** in `src/integrals/integrate.cpp`, extended the inverse-trig by-parts to
  the √-derivative functions (asin/acos/asinh/acosh): the standalone block now
  handles a bare power `f^n` (`u = f^n`, `dv = dx`), and the Mul block admits a
  `dv = P(x)/√(quadratic)` (so the residual `∫P/√Q` and its recursion close). Each
  return is gated by a **numeric diff-back self-check** — the broadened recursion
  threads several integrators, so verifying `d/dx == integrand` ensures any
  remaining mis-step fails to a clean marker rather than a wrong answer.
- **Verified:** `∫asin² = x·asin² − 2x + 2√(1−x²)·asin`, `∫x·asin²`, `∫acos²`,
  `∫asinh²`, `∫asin³` — all matching SymPy exactly; the non-elementary `∫atan²`
  bare stays an unevaluated marker.
- **Regression test:** extended `INT-32` in `tests/integrals/integrate_test.cpp`.

### INT-XSQRTQUAD-NUM-1 — `∫asin(x)/√(1−x²)` returned the wrong `asin(x)²` (should be `asin²/2`)
- **Problem:** a *wrong* (not merely unevaluated) answer: `∫asin(x)/√(1−x²) → asin(x)²`
  (correct is `asin²/2` — a factor-of-2 error), `∫asin²/√ → asin³`, `∫acos/√ → acos·asin`.
  `try_x_over_sqrt_quadratic` builds `Poly(numerator, var)` and reads its constant
  coefficient; for a *non-polynomial* numerator like `asin(x)`, Poly treats the whole
  thing as an opaque degree-0 "coefficient", so the handler pulled the var-dependent
  `asin(x)` out of the integral as if constant: `asin·∫1/√Q = asin·asin = asin²`.
- **Fix:** in `src/integrals/integrate.cpp`, `try_x_over_sqrt_quadratic` now rejects a
  numerator whose Poly coefficients depend on var (the check `try_poly_over_sqrt_quadratic`
  already had). The integrals then route to heurisch (`u = asin`), which gives the
  correct `asin²/2` — and a new **numeric diff-back self-check** was added to
  `try_heurisch` so any future mis-integration there fails to a clean marker rather
  than a wrong answer.
- **Verified:** `∫asin/√(1−x²) = asin²/2`, `∫asin²/√ = asin³/3`, `∫acos/√ = −acos²/2`
  all diff-back to the integrand, matching SymPy; the legitimate `∫x/√(1−x²) = −√(1−x²)`
  and `∫(2x+1)/√(1−x²)` are unchanged.
- **Regression test:** `INT-XSQRTQUAD-NUM-1` in `tests/integrals/integrate_test.cpp`.

### INT-INVTRIG-SQ-1 — `∫x·atan(x)²` (polynomial × inverse-trig squared) was unevaluated
- **Problem:** `∫x·atan(x)²` (= `x²·atan²/2 − x·atan + atan²/2 + log(x²+1)/2`) and
  `∫x·acot(x)²` were left unevaluated, though elementary. The inverse-trig by-parts
  block only matched a bare `f(affine)` (power 1) and required a *polynomial* `dv`.
- **Fix:** in `src/integrals/integrate.cpp`, the block now (a) matches a positive
  integer power `f^k` as the by-parts factor — `u = f^k` lowers the power by one
  each step, recursing to `f^1`; (b) for the rational-derivative functions
  (atan/acot/atanh/acoth) admits a *rational* `dv`, so the parts residual
  (`x²·atan/(1+x²)` for `∫x·atan²`) stays rational and closes; and (c) `expand`s the
  residual `v·f'` so a form like `(x−atan x)/(1+x²)` distributes for term-by-term
  integration. A recursive marker check bails (no partial garbage) when a branch
  doesn't reduce.
- **Verified:** `∫x·atan(x)² `, `∫x·acot(x)²` diff-back to the integrand, matching
  SymPy; bare `∫atan(x)²` (non-elementary) stays an unevaluated marker, and the
  earlier `∫atan/x²` / `∫x²·atan` cases are unchanged.
- **Regression test:** extended `INT-32` in `tests/integrals/integrate_test.cpp`.

### INT-RECIPTRIG-PARTS-1 — `∫x·sec²(x)` (= `∫x/cos²x`) and reciprocal-square trig were unevaluated
- **Problem:** `∫x/cos²(x)` (= `∫x·sec²x = x·tan x + log cos x`) and the family
  `∫x/sin²x`, `∫x/cosh²x`, `∫x/sinh²x` were left unevaluated. The polynomial × trig
  by-parts whitelist (`is_byparts_target`) only accepted *positive* integer powers
  of sin/cos/sinh/cosh, so a reciprocal (negative) power never matched — even though
  the antiderivative of the target (`∫1/cos² = tan`, …) is tabulated.
- **Fix:** in `src/integrals/integrate.cpp`, widened the whitelist to any non-zero
  integer power. Because an *odd* reciprocal power (e.g. `sec = 1/cos`) gives a
  non-elementary `∫v·u'` whose result is an `Add` with buried `Integral(...)` terms,
  the marker check was made **recursive** so the block bails to a clean marker
  instead of returning partial garbage.
- **Verified:** `∫x/cos²x = x·tan x + log cos x`, `∫x/sin²x`, `∫x/cosh²x`,
  `∫x/sinh²x` all diff-back to the integrand, matching SymPy; the non-elementary
  `∫x/cos x` stays an unevaluated marker (no garbage); positive-power cases
  (`∫x·cos²x`) unchanged.
- **Regression test:** `INT-RECIPTRIG-PARTS-1` in
  `tests/integrals/integrate_test.cpp`.

### INT-POLYSQRTQUAD-1 — `∫x²·√(1−x²)` (even power × √quadratic) was unevaluated
- **Problem:** `∫xⁿ·√(1−x²)` for *even* n — `∫x²·√(1−x²)`, `∫x⁴·√(1−x²)`,
  `∫x²·√(4−x²)` — was left unevaluated. The `u = Q` substitution closes the odd
  powers (`∫x·√(1−x²) = −(1−x²)^(3/2)/3`) but not the even ones.
- **Fix:** added `try_poly_times_sqrt_quadratic` in `src/integrals/integrate.cpp`
  (dispatched after the u-substitution handlers, so odd powers keep their cleaner
  form). It detects `P(x)·(quadratic)^(m/2)` with odd `m`, rewrites
  `P·Q^(m/2) = (P·Q^((m+1)/2))/√Q` — a polynomial over `√Q` — and hands it to the
  existing polynomial-over-√(quadratic) reduction.
- **Verified:** `∫x²·√(1−x²)`, `∫x⁴·√(1−x²)`, `∫x²·√(4−x²)`, `∫x²·(1−x²)^(3/2)` all
  diff-back to the integrand, matching SymPy; the odd-power `∫x·√(1−x²)` and bare
  `∫√(1−x²)` keep their existing forms.
- **Regression test:** `INT-POLYSQRTQUAD-1` in `tests/integrals/integrate_test.cpp`.

### LIMIT-LOGSUMEXP-1 — `(2^x+3^x)^(1/x) = 3` and log-of-exponential-sum limits failed
- **Problem:** limits of `log` of an exponential-dominated sum returned `nan` or
  an unevaluated ∞-arithmetic mess: `x − log(cosh x) = log 2`,
  `log(2^x+3^x)/x = log 3`, `(2^x+3^x)^(1/x) = 3` (the max of the bases). The engine
  folded the inner `log(∞)` directly and had no asymptotic for the sum.
- **Fix:** added `try_log_exp_asymptotic` in `src/calculus/limit.cpp`, run before
  direct substitution at `+∞`. For `log(g)` with `g` a sum of exponential terms
  (`cosh`/`sinh` and `a^x` first rewritten into `exp`), it picks the fastest-growing
  exponent `e_dom` (max coefficient of `x`) and rewrites
  `log(g) = e_dom + log(g·e^(−e_dom))`, where the residual `g·e^(−e_dom)` tends to a
  finite positive constant — so the residual `log` has a finite limit. The whole
  expression is expanded and re-limited (so a wrapping `log(g)/x` distributes
  instead of staying an `∞·0` product).
- **Verified:** `x − log(cosh x) → log 2`, `x − log(sinh x) → log 2`,
  `log(1+e^x) − x → 0`, `log(2^x+3^x)/x → log 3`, `(2^x+3^x)^(1/x) → 3`,
  `(3^x+5^x+2^x)^(1/x) → 5`, `log(e^(2x)+e^x)/x → 2`, all matching SymPy.
- **Regression test:** `LIMIT-LOGSUMEXP-1` in
  `tests/calculus/series_limit_test.cpp`.

### LIMIT-LHOPITAL-NEST-1 — `lim x·(π/2 − atan x) = 1` returned `nan`
- **Problem:** `0·∞` limits whose L'Hôpital derivative ratio has a denominator
  derivative that is itself a fraction — `x·(π/2 − atan x)`, `x·atan(1/x)`,
  `x·tan(1/x)` (all → 1) — returned `nan`. After differentiating, `d/dx(1/x) = −x⁻²`
  goes into the denominator, and the re-rationalisation step used `together()`,
  which does not flatten a nested reciprocal like `(−x⁻²)⁻¹`; the leftover negative
  power survived into the next substitution and produced `nan`. (`x·sin(1/x)`
  worked because the stray `x⁻²` happened to cancel against a matching factor.)
- **Fix:** in `src/calculus/limit.cpp`, `lhopital_nd` now rationalises each
  derivative ratio with `flatten_fraction(together(num'/den'))` — `together` cancels
  common factors (keeping `x·sin(1/x)` working) and a new recursive
  `flatten_fraction` helper (`(p/q)^(−k) = q^k/p^k`, descending into `Pow` bases)
  clears any residual nested reciprocal `together` leaves behind.
- **Verified:** `x·(π/2 − atan x) → 1`, `x·atan(1/x) → 1`, `x·tan(1/x) → 1`,
  matching SymPy; `x·sin(1/x) → 1`, `x·log(1+1/x) → 1`, and the existing rational /
  radical L'Hôpital limits are unchanged.
- **Regression test:** extended the `0·∞` test in
  `tests/calculus/series_limit_test.cpp`.

### INT-TRIGSQ-POWER-1 — `∫sin²(x)/xⁿ` and squared-trig over a power were unevaluated
- **Problem:** `∫sin²(x)/x²`, `∫cos²(x)/x`, `∫sin²(x)/x³`, … were left unevaluated.
  `try_trig_reduction` applied the half-angle identity only to a *standalone*
  `sin²(u)`, not to a `sin²(u)/cos²(u)` factor inside a product.
- **Fix:** in `src/integrals/integrate.cpp`, `try_trig_reduction` now also rewrites
  a `sin²(u)`/`cos²(u)` factor in a `Mul` via the half-angle identity
  (`sin²(u) = (1−cos 2u)/2`) and recurses: the integrand becomes `(1∓cos 2u)/(2·rest)`,
  which the linearity + `Si/Ci` power-reduction paths (`INT-EXPINT-POWER-1`) close.
  Gated to fire only when the remaining factors are non-trig (a power of `x`, an
  exponential, …) so a pure trig × trig product like `sin³·cos²` keeps its dedicated
  `sin^m·cos^n` closed form (which the rewrite would otherwise hijack into a messier
  result — a regression caught and fixed).
- **Verified:** `∫sin²(x)/x² = Si(2x) + cos(2x)/(2x) − 1/(2x)`, `∫cos²(x)/x²`,
  `∫sin²(x)/x = log(x)/2 − Ci(2x)/2`, `∫cos²(x)/x`, `∫sin²(x)/x³`, `∫sin²(2x)/x²` —
  all diff-back to the integrand, matching SymPy; `∫sin³·cos²` keeps
  `cos⁵/5 − cos³/3`, and standalone `∫sin²(x)` is unchanged.
- **Regression test:** extended `INT-EXPINT-POWER-1` in
  `tests/integrals/integrate_test.cpp`.

### INT-EXPINT-POWER-1 — `∫sin(x)/x²` and `∫f(c·x)/xⁿ` (n ≥ 2) were unevaluated
- **Problem:** `∫sin(x)/x²`, `∫cos(x)/x²`, `∫exp(x)/x²`, `∫sin(x)/x³`, … were left
  unevaluated. SymPP closed `∫f(c·x)/x = Si/Ci/Ei` (the n = 1 special-integral
  functions) but had no reduction for a higher reciprocal power.
- **Fix:** in `src/integrals/integrate.cpp`, `try_expint_integral` now matches a
  general reciprocal power `x^(−n)` (not just `x⁻¹`). For `n ≥ 2` it integrates by
  parts — `∫f(c·x)/xⁿ = f(c·x)/((1−n)·x^(n−1)) − c/(1−n)·∫f'(c·x)/x^(n−1)` —
  recursing on the residual (which is the same family with `n−1` and `f` replaced
  by its derivative) down to the `n = 1` Si/Ci/Ei base case. The marker guard bails
  if the residual doesn't close.
- **Verified:** `∫sin(x)/x² = Ci(x) − sin(x)/x`, `∫cos(x)/x² = −Si(x) − cos(x)/x`,
  `∫exp(x)/x² = Ei(x) − exp(x)/x`, `∫sin(x)/x³`, `∫sinh(x)/x² = Chi(x) − sinh(x)/x`,
  `∫sin(2x)/x²` — all diff-back to the integrand, matching SymPy; the `n = 1`
  Si/Ci/Ei cases are unchanged.
- **Regression test:** `INT-EXPINT-POWER-1` in `tests/integrals/integrate_test.cpp`.

### SUM-DIRICHLET-BETA-1 — `Σ (−1)^k/(2k+1)` (Leibniz) stayed unevaluated
- **Problem:** the Dirichlet beta series `Σ_{k=0}^∞ (−1)^k/(2k+1)^s` returned an
  unevaluated `Sum`. The Leibniz series `Σ(−1)^k/(2k+1) = π/4` and
  `Σ(−1)^k/(2k+1)² = Catalan` are clean closed forms SymPy produces.
- **Fix:** added a Dirichlet-beta branch in `src/calculus/summation.cpp` (next to
  the alternating p-series). For a summand `C·(−1)^(a·k+b)·(2k+1)^(−s)` (`a` odd,
  `b` integer, base verified to be exactly `2·var+1`) over `k = 0…∞`, it returns
  `π/4` for `s = 1` and Catalan's constant for `s = 2`, with the sign `(−1)^b` and
  leading constant multiplied through. Higher `s` (no elementary form — SymPy
  gives a polylog) are left unevaluated.
- **Verified:** `Σ(−1)^k/(2k+1) = π/4`, `Σ(−1)^(k+1)/(2k+1) = −π/4`,
  `Σ 2(−1)^k/(2k+1) = π/2`, `Σ(−1)^k/(2k+1)² = Catalan`, matching SymPy; `s = 3`
  and non-`(2k+1)` denominators (`3k+1`) stay unevaluated; the alternating
  k-denominator (eta) series and all other sums unchanged.
- **Regression test:** `SUM-DIRICHLET-BETA-1` in
  `tests/calculus/series_limit_test.cpp`.

### SUM-ALT-PSERIES-1 — `Σ (−1)^k/k` and alternating p-series stayed unevaluated
- **Problem:** the alternating p-series `Σ_{k=1}^∞ (−1)^k/k^s` — `Σ(−1)^k/k = −log 2`,
  `Σ(−1)^k/k² = −π²/12`, `Σ(−1)^k/k³ = −¾ζ(3)` — returned an unevaluated `Sum`.
  Only the non-alternating `Σ1/k^s = ζ(s)` was handled.
- **Fix:** added an alternating-p-series branch in `src/calculus/summation.cpp`
  (next to the ζ p-series). It recognizes a summand `C·(−1)^(a·k+b)·k^(−s)` with `a`
  an odd integer (so `(−1)^(a·k) = (−1)^k`) and `b` an integer (constant sign
  `(−1)^b`), and returns the Dirichlet eta value: `−log 2` for `s = 1`, and
  `(2^(1−s) − 1)·ζ(s)` for `s ≥ 2` (closing to a π-power for even `s`). A leading
  constant multiplies through.
- **Verified:** `Σ(−1)^k/k = −log 2`, `Σ(−1)^(k+1)/k = log 2`, `Σ(−1)^k/k² = −π²/12`,
  `Σ(−1)^k/k⁴ = −7π⁴/720`, `Σ(−1)^k/k³ = −¾ζ(3)` (= SymPy's `−η(3)`),
  `Σ 3(−1)^k/k = −3 log 2`, all matching SymPy; non-alternating p-series
  (`Σ1/k² = π²/6`) and divergent/other sums unchanged.
- **Regression test:** `SUM-ALT-PSERIES-1` in `tests/calculus/series_limit_test.cpp`.

### SOLVE-INVFN-SYM-1 — `solve(atan(x) − a)` (inverse fn = symbolic RHS) returned `[]`
- **Problem:** inverting an inverse trig/hyperbolic function against a *symbolic*
  right-hand side returned `[]`: `solve(atan(x) − a) → []`, `asin(x) − a`,
  `acos(x) − a`, … although the numeric case worked (`atan(x) − 1 → tan(1)`).
- **Fix:** in `src/solvers/solve.cpp`, `solve_inverse_func_poly`'s `in_range` check
  no longer rejects a non-numeric angle `c` for the bounded-range functions
  (asin/acos/atan). A concrete out-of-range value is still rejected; a symbolic `c`
  now yields the formal principal-branch inverse, matching SymPy.
- **Verified:** `atan(x) − a → tan(a)`, `asin(x) − a → sin(a)`, `acos(x) − a → cos(a)`,
  `atanh(x) − a → tanh(a)`, `asinh(x) − a → sinh(a)`, `atan(2x) − a → tan(a)/2`,
  all matching SymPy; numeric in-range (`atan(x) − 1 → tan(1)`) and out-of-range
  rejection (`asin(x) − 5 → []`) unchanged. As a knock-on, the ODE `y′ = 1 + y²`
  now solves explicitly to `tan(x + C)` (was an implicit `atan` form).
- **Regression test:** extended `SOLVE-INVFN-1` in `tests/solvers/solve_test.cpp`.

### DSOLVE-SEPARABLE-NONLIN-1 — `dsolve(y′ = y²)` and nonlinear separable ODEs were unsolved
- **Problem:** separable equations `y′ = f(x)·g(y)` with a non-trivial `g(y)` —
  `y′ = y²`, `y′ = √y`, `y′ = x·y²`, `y′ = 1 + y²` — returned an unevaluated
  `Dsolve(…)` marker. `try_separate` decided x-dependence with `has(rhs, x)`, but
  the dependent variable is the *function application* `y(x)`, which literally
  contains the symbol `x` — so `has(y², x)` was wrongly `true` and separation
  failed for every autonomous/`g(y)`-only right-hand side.
- **Fix:** in `src/ode/dsolve.cpp`, `try_separate` now tests x-dependence with `y`
  replaced by a fresh atom (`has(subs(·, y, u), x)` — "depends on x with y held
  fixed"). The explicit-form `solve` fallback also swaps `y(x)` for a plain symbol
  so `solve`'s inverters can isolate `y`. `dsolve_separable` is moved *after* the
  linear/Bernoulli/homogeneous methods, which give cleaner closed forms for the
  equations they recognize (the logistic stays explicit), so separation only
  fills the gaps it uniquely covers.
- **Verified:** `y′ = y² → −1/(x+C)`, `y′ = x·y² → −2/(x²+C)`, `y′ = √y → ((x+C)/2)²`
  (residuals 0, matching SymPy); `y′ = 1+y²`, `y′ = y²−1` are now solved (implicitly
  via atan/log, since `solve` can't invert against a symbolic RHS) rather than
  unevaluated; the logistic `y′ = y(1−y)` stays explicit and all existing dsolve
  results are unchanged.
- **Regression test:** `DSOLVE-SEPARABLE-NONLIN-1` in `tests/ode/dsolve_test.cpp`.

### DSOLVE-RESONANCE-1 — `dsolve(y″ + y = sin x)` returned garbage with `zoo`
- **Problem:** a second-order constant-coefficient ODE whose forcing term is
  itself a homogeneous solution (resonance) — `y″ + y = sin x`, `y″ + 4y = sin 2x`
  — produced garbage like `… + zoo·I·cos(x) + zoo·sin(x)`. Variation of parameters
  used the *complex* basis `e^(±iβx)`, and the cyclic exp·trig integrator
  `∫e^(ax)sin(gx) = e^(ax)(a·sin − g·cos)/(a²+g²)` divides by `a²+g²`, which is `0`
  for `a = −i, g = 1` at resonance → `zoo`.
- **Fix:** in `src/ode/dsolve.cpp`, both `order2_basis` and `dsolve_constant_coeff`
  now emit the **real** basis `e^(αx)cos(βx), e^(αx)sin(βx)` for a complex-conjugate
  root pair `α ± βi` (paired via `simplify(rootⱼ − conj) == 0`, robust to nested
  radical roots), instead of complex exponentials. The real basis keeps the
  variation-of-parameters integrals real, so the `x·(…)` resonance factor falls out
  correctly, and the homogeneous solution now matches SymPy's `C₁cos + C₂sin` form.
- **Verified:** `y″ + y = sin x` → residual 0 (particular part carries `−x·cos x/2`),
  `y″ + y = cos x`, `y″ + 4y = sin 2x` all residual 0 with no `zoo`; `y″ + 4y = 0`
  → `C₀cos 2x + C₁sin 2x` (real, no `I`); `y‴ − y = 0` →
  `e^(−x/2)(C₁cos(√3x/2) + C₂sin(√3x/2)) + C₀eˣ`, matching SymPy.
- **Regression test:** `DSOLVE-RESONANCE-1` in `tests/ode/dsolve_test.cpp`.

### INT-ABS-DEF-1 — definite integral of `|x|` returned garbage `−Integral(1,−1) + Integral(1,1)`
- **Problem:** definite integrals of integrands containing `abs`/`sign` —
  `∫_{-1}^1 |x|`, `∫_0^π |cos x|`, `∫_{-1}^1 (|x|+x²)` — produced garbage. These
  have no elementary antiderivative (SymPy leaves the *indefinite* form too), so
  the Newton–Leibniz path substituted the bounds into the unevaluated
  `Integral(|x|, x)` marker, yielding nonsense like `−Integral(1,−1)+Integral(1,1)`.
- **Fix:** added `try_abs_definite` in `src/integrals/integrate.cpp`, invoked from
  the 4-arg `integrate` when the antiderivative still contains an integral marker
  (detected recursively, since it can be buried in a sum). `|g|` and `sign(g)` are
  piecewise-constant in the sign of `g`, so it splits `[lower, upper]` at the real
  roots of each argument (via `solve`), replaces `abs(g)→±g` / `sign(g)→±1` by the
  numerically-sampled sign on each subinterval, integrates the now-smooth pieces,
  and sums. Finite bounds only; bails unless every piece closes.
- **Verified:** `∫_{-1}^1 |x| = 1`, `∫_{-2}^3 |x| = 13/2`, `∫_0^2 |x−1| = 1`,
  `∫_{-1}^1 x|x| = 0`, `∫_{-1}^2 sign x = 1`, `∫_{-1}^1 (|x|+x²) = 5/3`,
  `∫_0^π |cos x| = 2`, all matching SymPy; a smooth integrand (no interior root)
  reduces to the ordinary integral.
- **Regression test:** `INT-ABS-DEF-1` in `tests/integrals/integrate_test.cpp`.

### SUM-POLYGEOM-INF-1 — `Σ_{k=0}^∞ k/2^k` returned `nan`
- **Problem:** infinite polynomial × geometric sums `Σ_{k=lo}^∞ P(k)·r^k` with
  `|r| < 1` — `Σ k/2^k = 2`, `Σ k²/2^k = 6`, `Σ k/3^k = 3/4` — returned `nan`.
  The closed form sums the antidifference `Q(k)·r^k` and evaluates the upper
  boundary at `k = ∞` as `(polynomial in ∞)·r^∞ = ∞·0 = nan`.
- **Fix:** in `src/calculus/summation.cpp`, `sum_poly_geometric` now treats an
  infinite upper bound specially: when `|r| < 1` (`r² < 1` provable), the boundary
  term `Q(k)·r^k → 0` (geometric decay dominates the polynomial), so the sum is
  `−S(lo)`. A divergent or undecidable ratio is left as an unevaluated `Sum`
  (not a bogus value). The degree-1 arithmetic-geometric block now defers infinite
  bounds to this path so both go through the convergence handling.
- **Verified:** `Σ_{k=0}^∞ k/2^k = 2`, `Σ_{k=1}^∞ k/2^k = 2`, `Σ k²/2^k = 6`,
  `Σ k/3^k = 3/4`, `Σ (k+1)/2^k = 4`, matching SymPy; the divergent `Σ k·2^k`
  stays unevaluated; finite sums (`Σ_{k=1}^n k·2^k`, `Σ_{k=1}^3 k²·2^k`) unchanged.
- **Regression test:** `SUM-POLYGEOM-INF-1` in
  `tests/calculus/series_limit_test.cpp`.

### SOLVE-RADISOLATE-1 — `solve(√(x+1) − x + 1)` and single-radical equations returned `[]`
- **Problem:** equations with a single square root of a non-trivial radicand
  appearing linearly — `√(x+1) − x + 1 = 0`, `√(2x+3) − x = 0`, `√(x+1) + x = 0`,
  `√(x²+1) − x − 1 = 0` — returned `[]`. `solve_radical_poly` only handles a
  polynomial in `x^(1/d)` of the *bare* variable, so a radical of an affine /
  quadratic argument fell through.
- **Fix:** added `solve_radical_isolate` in `src/solvers/solve.cpp`. It detects a
  lone `√(g(x))`, linearizes the equation in it (`A·√g + B = 0`, `A, B`
  radical-free), squares to the polynomial `A²·g − B² = 0`, solves that, and
  filters the roots back through the *original* equation to drop the extraneous
  ones introduced by squaring. The filter is **numeric** (`evalf`, accepting an
  exact `0` or `|·| < 1e-20`): a symbolic check can't denest forms like
  `√((3−√5)/2) = (√5−1)/2`, which would wrongly reject the valid root.
- **Verified:** `√(x+1) − x + 1 → {3}`, `√(2x+3) − x → {3}`,
  `√(x²+1) − x − 1 → {0}`, `√(x+1) + x → {(1−√5)/2}` (the `(1+√5)/2` branch
  correctly dropped), `√(x+1) − x − 1 → {−1, 0}`, all matching SymPy; no-solution
  cases (`√(x+1) + 2`) stay `[]`; `√(x+1) − 2 → {3}` and `x − √x − 2 → {4}`
  unchanged.
- **Regression test:** `SOLVE-RADISOLATE-1` in `tests/solvers/solve_test.cpp`.

### INT-LOGSUB-1 — `∫cos(log x)`, `∫log(log x)/x` and log-composite integrands were unevaluated
- **Problem:** integrands built from `log(x)` — `∫cos(log x)`, `∫sin(log x)`,
  `∫cos(2·log x)`, `∫log(log x)/x` — were left unevaluated, though each is
  elementary under the substitution `u = log(x)`.
- **Fix:** added `try_log_substitution` in `src/integrals/integrate.cpp`
  (dispatched after integration-by-parts, before the Weierstrass path). When
  `log(var)` appears, it substitutes `u = log(x)` (`x = eᵘ`, `dx = eᵘ du`) by
  replacing `log(var) → u` and every remaining bare `var → eᵘ`, leaving
  `∫f(u)·eᵘ du`, which it integrates and back-substitutes. A surviving `var` (e.g.
  `log(2x)`, not the `log(x)` node) means the substitution is incomplete and it
  bails. The gate on `log(var)` keeps ordinary integrands untouched.
- **Verified:** `∫cos(log x) = x(cos(log x)+sin(log x))/2` (a cyclic exp·trig
  integral in `u`), `∫sin(log x)`, `∫cos(2·log x) = x(cos(2log x)+2sin(2log x))/5`,
  `∫log(log x)/x = log x·log(log x) − log x` (becomes `∫log u`) — all diff-back to
  the integrand, matching SymPy; `∫1/x`, `∫x·log x` unchanged.
- **Regression test:** `INT-LOGSUB-1` in `tests/integrals/integrate_test.cpp`.

### SERIES-COMPOSE-1 — `series(log(sin x / x))` stayed unexpanded
- **Problem:** the Taylor series of a composite `f(g(x))` whose inner `g` is finite
  but non-simple at the expansion point — e.g. `g = sin(x)/x`, with its `1/x`
  factor — was not produced. `taylor_series` differentiates `f(g)` directly and
  evaluates each derivative via a `limit`; for such `g` those derivative-limits
  get hard (`lim (log(sin x/x))'' = nan`), so the Taylor path bailed and
  `series(log(sin x / x))` returned the input unexpanded. (This was the underlying
  cause worked around in `LIMIT-POWFORM-COMPOSITE-1`.) Single-function bases like
  `log(cos x)` worked because their derivatives stay simple.
- **Fix:** added `try_composition_series` in `src/calculus/series.cpp` (between the
  Taylor and Laurent paths). It expands the inner `g` with the full `series()`
  engine (so an inner that itself needs Laurent division, e.g. `tan x / x`, still
  expands), requires `g` analytic at `x₀` (`c = g(x₀)` finite — which rejects a
  genuine pole like `csc x` whose inner series diverges), expands the outer about
  the constant `c`, and substitutes `(t − c) → (g_series − c)` — positive
  valuation, so only finitely many terms reach a given order — then truncates. The
  outer operation is a single-argument function `f` *or* a power `g^p` with a
  var-free exponent (covers `√(tan x / x)`). A genuine singularity (`log x`, `√x`,
  where the outer Taylor at `c = 0` fails) still stays unexpanded.
- **Verified:** `series(log(sin x / x)) = −x²/6 − x⁴/180`,
  `series(log(sinh x / x)) = x²/6 − x⁴/180`, `series(log(tan x / x)) = x²/3 + 7x⁴/90`,
  `series(√(tan x / x)) = 1 + x²/6 + 19x⁴/360`, matching SymPy; `log x` / `√x`
  unexpanded; the `csc²(x)` Laurent series and ordinary/single-function series
  (`exp`, `sin`, `log(cos x)`, `cot`, …) unchanged.
- **Regression test:** `SERIES-COMPOSE-1` in
  `tests/calculus/series_limit_test.cpp`.

### LIMIT-POWFORM-COMPOSITE-1 — `(sin x / x)^(1/x²)` returned `nan` instead of `e^(−1/6)`
- **Problem:** `1^∞` limits whose base tends to 1 through a *composite* expression —
  `(sin x / x)^(1/x²) → e^(−1/6)`, `(tan x / x)^(1/x²) → e^(1/3)` — returned `nan`.
  `try_power_form` resolves `a^b` via `exp(lim b·log a)`, but the inner limit
  `lim log(sin x / x)/x²` failed: the series engine leaves `log(sin x / x)`
  (a log of a quotient) unexpanded, so the `0/0` rate could not be evaluated.
  Single-function bases like `cos(x)^(1/x²)` worked because `log(cos x)` does expand.
- **Fix:** in `src/calculus/limit.cpp`, the `1^∞` branch of `try_power_form` now
  uses the rate `b·(a−1)` instead of `b·log a`. As `a → 1`,
  `log a = (a−1) − (a−1)²/2 + … = (a−1)·(1 + o(1))`, so `lim b·log a = lim b·(a−1)`
  exactly (the correction is `b·(a−1)·(a−1) → 0`). This sidesteps the missing
  log-of-composite series entirely. The `∞^0` and `0^0` forms genuinely need
  `log a` and keep it.
- **Verified:** `(sin x/x)^(1/x²) → e^(−1/6)`, `(tan x/x)^(1/x²) → e^(1/3)`,
  `cos(2x)^(1/x²) → e^(−2)`, `(1+sin x)^(1/x) → e`, matching SymPy; the existing
  `(1+a/x)^x → eᵃ`, `cos(x)^(1/x²)`, `x^x → 1` cases are unchanged.
- **Regression test:** extended `LIMIT-POWFORM-1` in
  `tests/calculus/series_limit_test.cpp`.

### INT-INVTRIG-RECIP-1 — `∫atan(x)/x²` and inverse-trig over a reciprocal power were unevaluated
- **Problem:** `∫atan(x)/x²`, `∫atan(x)/x³`, `∫atanh(x)/x²`, `∫acot(x)/x²` were left
  unevaluated, although they are elementary (by parts the residual is rational).
  The polynomial × by-parts-function block required the `dv` factor to be a
  *polynomial* (`is_polynomial_in(rest, var)`), so a reciprocal power `1/xⁿ` failed
  the gate even though `∫x^(−n)` is elementary.
- **Fix:** in `src/integrals/integrate.cpp`, the block now also admits a bare
  reciprocal power `dv = x^(−n)`, but only for the inverse functions with a
  *rational* derivative — atan/acot/atanh/acoth — where the by-parts residual
  `v·f'` stays rational and `try_rational` closes it exactly. The √-derivative
  functions (asin/acos/asinh/acosh) keep the polynomial-only gate: over a `1/x`
  factor their residual is non-rational and the rational path silently mis-handled
  it (`∫asin(x)/x²` collapsed to a bogus `0`). The marker guard still bails on the
  genuinely non-elementary `n = 1` case (`∫atan(x)/x`, residual `log(x)/(x²+1)`).
- **Verified:** `∫atan(x)/x² = log x − ½log(x²+1) − atan(x)/x`, `∫atan(x)/x³`,
  `∫atanh(x)/x²`, `∫acot(x)/x²` — all diff-back to the integrand (numeric), matching
  SymPy; `∫atan(x)/x` and `∫asin(x)/x²` correctly stay unevaluated.
- **Regression test:** extended `INT-32` in `tests/integrals/integrate_test.cpp`.

### INT-CONSTBASEEXP-1 — `∫2ˣ` and `∫P(x)·aˣ` (constant-base exponential) were unevaluated
- **Problem:** `∫2ˣ`, `∫x·2ˣ`, `∫x²·2ˣ` and every `∫P(x)·a^(b·x+c)` with a constant
  base `a ≠ e` were left unevaluated. SymPP integrated the natural base `eˣ` but had
  no rule for `aˣ`; rewriting `aˣ = exp(x·ln a)` does not help because that form
  canonicalizes straight back to `a^x`.
- **Fix:** added `try_const_base_exp_integral` in `src/integrals/integrate.cpp`
  (dispatched with the other special-exponential rules). It isolates a constant,
  provably-positive base power `a^(affine)` (`a ≠ 1`, exponent affine in var) and a
  polynomial residual, then integrates each monomial by the by-parts reduction
  `∫xⁿ·a^g = xⁿ·a^g/k − (n/k)·∫xⁿ⁻¹·a^g` with `k = b·ln a`, bottoming out at
  `∫a^g = a^g/(b·ln a)`. The natural base `eˣ` (a `Function`, not a `Pow`) is not
  matched, so the existing elementary path for it is untouched.
- **Verified:** `∫2ˣ = 2ˣ/ln 2`, `∫x·2ˣ = 2ˣ(x·ln 2 − 1)/ln²2`, `∫x²·2ˣ`, `∫x·3ˣ`,
  `∫(x+1)·2ˣ`, `∫x·2^(−x)`, `∫2^(3x)` — all diff-back exactly to the integrand,
  matching SymPy; `∫x·eˣ` unchanged.
- **Regression test:** `INT-CONSTBASEEXP-1` in
  `tests/integrals/integrate_test.cpp`.

### SOLVE-ZEROPROD-1 — `solve(x²·eˣ − eˣ)` returned `[]`; `eˣ·(x²−4)` gave a spurious `zoo`
- **Problem:** Equations that factor into a polynomial × transcendental were
  mis-solved. `x²·eˣ − eˣ` returned `[]` (the common `eˣ` is not polynomial, so the
  Poly path could not see `eˣ·(x²−1)`); `eˣ·(x²−4)` returned `[2, −2, zoo]` — the
  spurious `zoo` from solving the never-zero factor `eˣ = 0`; and `x·cos(x)`
  returned only `[0]` because `solve_poly` read it as a degree-1 polynomial whose
  coefficient happened to be `cos(x)`.
- **Fix:** added `solve_zero_product` in `src/solvers/solve.cpp`. A product (or an
  `Add` with a common factor, found by intersecting the per-term factor maps)
  vanishes iff one factor does, so it solves each factor recursively and unions
  the roots — skipping factors that can never be zero (`is_never_zero`: `exp(·)`
  and nonzero constants) and denominator factors (negative powers, whose zeros are
  poles excluded from the surviving roots). It runs ahead of `solve_poly` when a
  function/radical of the variable is present (so the partial polynomial reading
  no longer wins) and again after, for the common-factor `Add` case.
- **Verified:** `x²·eˣ − eˣ → {1,−1}`, `eˣ·(x²−4) → {2,−2}` (no `zoo`),
  `x²·sin x − sin x → {0,π,1,−1}`, `x³·eˣ − x·eˣ → {0,1,−1}`,
  `x·cos x → {0,π/2,3π/2}`, `sin x·(x−1) → {0,1,π}`,
  `eˣ·(x²−1)·(x−3) → {1,−1,3}` — all matching SymPy; the removable-pole case
  `(x²−1)/(x−1) → {−1}` and plain polynomials are unchanged.
- **Regression test:** `SOLVE-ZEROPROD-1` in `tests/solvers/solve_test.cpp`.

### INT-GAUSSSHIFT-1 — `∫exp(−(x−1)²)` and Gaussians with a linear term were unevaluated
- **Problem:** `∫exp(−(x−1)²)`, `∫exp(−x²+x)`, `∫x·exp(−(x−1)²)` and every
  `∫P(x)·exp(a·x²+b·x+c)` with a non-zero linear term (`b ≠ 0`) were left
  unevaluated. The Gaussian rules (`try_gaussian`, `try_poly_times_gaussian`)
  require a *pure*-quadratic exponent (`b = c = 0`); a linear term needs completing
  the square first, which nothing did.
- **Fix:** added `try_shifted_gaussian` in `src/integrals/integrate.cpp`
  (dispatched just before `try_gaussian`). It isolates the `exp(quadratic)` factor
  and a polynomial residual, completes the square
  `a·x²+b·x+c = a·(x−x₀)² + K` with `x₀ = −b/(2a)`, `K = c − b²/(4a)`, substitutes
  `u = x − x₀` (so the exponent becomes the pure Gaussian `e^K·exp(a·u²)`), and
  delegates back to `integrate()` in `u` — reusing the moment/erf rules — before
  back-substituting. The recursion terminates because the shifted exponent has
  `b = 0`, so it never re-enters `try_shifted_gaussian`.
- **Verified:** `∫exp(−(x−1)²) = √π·erf(x−1)/2`, `∫exp(−x²+x) = √π·e^(1/4)·erf(x−1/2)/2`,
  `∫x·exp(−(x−1)²)`, `∫exp(x²+x)` (erfi), `∫exp(−2x²+3x−1)` — all diff-back exactly
  to the integrand, matching SymPy; pure-quadratic cases unchanged.
- **Regression test:** `INT-GAUSSSHIFT-1` in
  `tests/integrals/integrate_test.cpp`.

### INT-GAUSSMOMENT-1 — `∫x²·exp(−x²)` (polynomial × Gaussian) was unevaluated
- **Problem:** `∫x²·exp(−x²)` and every `∫P(x)·exp(c·x²)` with a non-constant
  polynomial `P` were left unevaluated (`Integral(…)` marker). SymPP integrated the
  bare Gaussian `∫exp(−x²) = √π·erf(x)/2` but had no rule for the Gaussian
  *moments*. The improper form was worse: `∫₀^∞ x²·exp(−x²)` evaluated the missing
  antiderivative at the bounds and emitted garbage `−Integral(0,0) + Integral(nan, ∞)`.
- **Fix:** added `try_poly_times_gaussian` in `src/integrals/integrate.cpp`
  (dispatched just before `try_gaussian`). It isolates the `exp(c·x²)` factor
  (pure quadratic exponent, provable-sign `c`) and a polynomial residual, then
  integrates each monomial via the by-parts reduction
  `∫xⁿ·exp(c·x²) = xⁿ⁻¹·exp(c·x²)/(2c) − (n−1)/(2c)·∫xⁿ⁻²·exp(c·x²)`, bottoming out
  at `∫exp(c·x²)` (erf/erfi) for even `n` and `∫x·exp(c·x²) = exp(c·x²)/(2c)` for
  odd `n`. Covers negative `c` (erf) and positive `c` (erfi).
- **Verified:** `∫x²·exp(−x²) = −x·exp(−x²)/2 + √π·erf(x)/4`, `∫x³·exp(−x²)`,
  `∫x⁴·exp(−x²)`, `∫(x²+1)·exp(−x²)`, `∫x²·exp(x²)` (erfi) — all diff-back exactly to
  the integrand; the improper `∫₀^∞ x²·exp(−x²) = √π/4`, matching SymPy.
- **Regression test:** `INT-GAUSSMOMENT-1` in
  `tests/integrals/integrate_test.cpp`.

### LIMIT-EXPPOLY-1 — `lim x²·(2/3)^x` and polynomial × exponential-ratio returned `nan`
- **Problem:** `lim_{x→∞} x²·(2/3)^x` (= 0), `x³·2^x/3^x` (= 0), `x²/2^x` (= 0),
  `x²·3^x/2^x` (= ∞) returned `nan`. The generic product/L'Hôpital path closed a
  degree-1 polynomial against a rational-base exponential (`x·(1/2)^x → 0`) but
  stalled at degree ≥ 2 — each L'Hôpital step lowers the polynomial degree by one
  while reproducing the rational-base exponential, and the recursion did not
  converge (natural-base `x^n·e^(−x)` worked, via the exp-aware reciprocal).
- **Fix:** extended `try_exponential_product` (see `LIMIT-EXPRATIO-1`) to accept a
  residual factor, required to be a polynomial in var. The exponential's growth
  class strictly dominates any polynomial, so once the combined base `B` is known:
  a decaying `B^m` (→ 0) drives the whole product to 0 regardless of polynomial
  degree, and a growing `B^m` (→ +∞) gives ±∞ with the sign of the polynomial
  residual's divergence. A non-polynomial residual is rejected (left to other
  paths) so the dominance argument stays valid.
- **Verified:** `x²·(2/3)^x → 0`, `x³·2^x/3^x → 0`, `x²/2^x → 0`,
  `x²·3^x/2^x → ∞`, `−x²·3^x → −∞`, all matching SymPy; pure exponential ratios
  and `x^n·e^(−x)` unchanged.
- **Regression test:** extended `LIMIT-EXPRATIO-1` in
  `tests/calculus/series_limit_test.cpp`.

### LIMIT-EXPRATIO-1 — `lim 2^x/3^x` and other exponential ratios returned `nan`
- **Problem:** `lim_{x→∞} 2^x/3^x` (= 0), `3^x/2^x` (= ∞), `exp(x)/exp(2x)`,
  `2^x·e^(−3x)` and similar returned `nan`. Each is a product/ratio of distinct
  constant-base exponentials; the limit engine evaluated the factors
  independently (`2^x → ∞`, `3^(−x) → 0`) and saw an `∞·0` indeterminate that
  L'Hôpital cannot crack — differentiating reproduces the same form — so the
  product path stalled and returned `nan`. A single `(2/3)^x` worked, because it
  is one power, not a product.
- **Fix:** added `try_exponential_product` in `src/calculus/limit.cpp`, run before
  the generic product path for `Mul` at `±∞`. When every factor is a constant-base
  exponential `bᵢ^(cᵢ·m)` or `exp(dⱼ·m)` (incl. `exp(g)^k`, the canonical form of
  `1/exp(g)`) sharing one var-monomial `m`, it folds them into a single `B^m` with
  `B = ∏bᵢ^cᵢ·e^(Σdⱼ)` a concrete positive constant, and decides the limit from
  `sign(B−1)` and the direction of `m` (numeric `evalf` fallback signs `B` when the
  base carries an `exp`, e.g. `exp(−1)−1`). A polynomial residual factor is handled
  by growth dominance — see `LIMIT-EXPPOLY-1`.
- **Verified:** `2^x/3^x → 0`, `3^x/2^x → ∞`, `exp(x)/exp(2x) → 0`,
  `2^x·e^(−3x) → 0`, `2^x·2^(−x) → 1`, all matching SymPy at `+∞`. At `−∞` the
  direction flips correctly (`2^x/3^x → ∞`); note SymPy is itself internally
  inconsistent there (`limit((2/3)**x,−∞)=0` vs `limit((2/3)**(−x),∞)=∞`), and the
  numeric values confirm SymPP's `∞` is the correct branch.
- **Regression test:** `LIMIT-EXPRATIO-1` in
  `tests/calculus/series_limit_test.cpp`.

### INT-WEIERSTRASS-NUM-1 — `∫cos(x)/(1+cos x)` and numerator-bearing rational trig unevaluated
- **Problem:** `∫cos(x)/(1+cos x)` (SymPy: `x − tan(x/2)`) was left unevaluated.
  Same root cause as `INT-WEIERSTRASS-DEGEN-1`, but worse: with a non-constant
  numerator the half-angle substitution produces an integrand whose denominator
  is itself a fraction `1 + (1−t²)/(1+t²)` *inside a `Pow` base*. Neither
  `together()` nor `cancel()` descends into a `Pow` base, so the integrand stayed
  a nested fraction and `try_rational` could not integrate it → unevaluated
  marker.
- **Fix:** in `src/integrals/integrate.cpp`, added a file-local `flatten_ratio`
  helper that recursively decomposes a finite rational expression into a single
  numerator/denominator pair, descending into integer-power bases
  (`(p/q)^(−k) = q^k/p^k`). `try_weierstrass` now flattens the substituted
  integrand with it before `cancel()`. The recursion is deliberately *not* added
  to the library `as_numer_denom()` — doing so globally perturbs the limit engine
  when a base carries infinities (e.g. `limit((1+a/x)^x) = e^a`); the
  Weierstrass-substituted integrand is always a finite rational function of `t`,
  so the local helper is both safe and sufficient.
- **Verified:** `∫cos(x)/(1+cos x) = −tan(x/2) + 2·atan(tan x/2) = x − tan(x/2)`
  (diff-back numerically exact; matches SymPy `x − tan(x/2)`);
  `∫(2+cos x)/(1+cos x) = x + tan(x/2)`. All prior `∫1/(a+b·cos x)` cases unchanged.
- **Regression test:** numeric diff-back block added to the Weierstrass test
  (INT-33) in `tests/integrals/integrate_test.cpp`.

### INT-WEIERSTRASS-DEGEN-1 — `∫1/(1+cos x)` returned garbage `zoo·log 2`
- **Problem:** `∫1/(1+cos x)` returned `zoo·log(2)` instead of `tan(x/2)`. The
  half-angle (Weierstrass) substitution `t = tan(x/2)` maps `1/(1+cos x)` to the
  constant integrand `1`, but `try_weierstrass` used `together()` to form the
  substituted integrand, and for this degenerate `a=b` case `together()` left a
  nested, non-reduced denominator `((1−t²)/(1+t²) + 1)·(1+t²)` — which only
  collapses to the constant `2` after full cancellation. Handing that un-reduced
  form to `integrate()` made `try_rational` misparse the denominator and emit
  `zoo`. (`1/(2+cos x)`, `1/(1−cos x)`, `1/(1±sin x)` etc. reduce cleanly under
  `together` and were unaffected.)
- **Fix:** in `src/integrals/integrate.cpp`, `try_weierstrass` now builds the
  integrand by flattening it to a single numerator/denominator with the
  `flatten_ratio` helper (see `INT-WEIERSTRASS-NUM-1`) and then `cancel()`-ing to
  lowest terms, instead of bare `together(...)`. The `has_trig_power_of`
  early-return still backstops the runaway-on-trig-powers case that motivated
  `together`.
- **Verified:** `∫1/(1+cos x) = tan(x/2)` (diff-back is exactly `1/(1+cos x)`,
  matches SymPy); all other `∫1/(a+b·cos x)`, `∫1/(a+b·sin x)` cases unchanged.
- **Regression test:** added the `a=b` cosine case to the Weierstrass oracle
  diff-back set in `tests/integrals/integrate_test.cpp` (INT-33).

### INT-QUADLOG-PARAM-1 — `∫1/(a²−x²)` unevaluated for symbolic coefficients
- **Problem:** `∫1/(a²−x²)` and `∫1/(x²−a²)` (negative-discriminant quadratics,
  the log/atanh case) were unevaluated for symbolic positive coefficients. The
  log branch of `try_arctan_quadratic` carried a `rational_coeffs` gate, even
  though it already requires `is_positive(Δ)` (Δ = b²−4ac).
- **Fix:** in `src/integrals/integrate.cpp`, dropped the `rational_coeffs` gate
  on the log branch; it fires whenever the discriminant is *provably negative*
  (Δ provably positive), e.g. `1/(a²−x²)` with `a > 0` (Δ = 4a²). Completes the
  parametric quadratic-integral family with `INT-ARCTAN-PARAM-1`.
- **Verified:** `∫1/(a²−x²) = log((a+x)/(x−a))/(2a)`,
  `∫1/(x²−a²)` — diff-back verified at concrete positive values; numeric
  quadratics and the arctan branch are unchanged.
- **Regression test:** extended `INT-ARCTAN-PARAM-1` in
  `tests/integrals/integrate_test.cpp`.

### INT-GAUSS-PARAM-1 — parametric Gaussian `∫exp(−a·x²)` unevaluated
- **Problem:** `∫exp(−a·x²)` and `∫exp(a·x²)` were unevaluated for a symbolic
  positive coefficient — `try_gaussian` already branched on `is_negative`/
  `is_positive(c₂)` but a leftover `is_rational(c₂)` gate blocked symbolic ones.
- **Fix:** in `src/integrals/integrate.cpp`, removed the `is_rational(c₂)` gate
  in `try_gaussian`; the sign branches decide erf vs erfi. Same pattern as
  `INT-ARCTAN-PARAM-1` / `INT-SQRTQUAD-PARAM-1`.
- **Verified:** `∫exp(−a·x²) = √π·erf(√a·x)/(2√a)`,
  `∫exp(a·x²) = √π·erfi(√a·x)/(2√a)` (a > 0) — match SymPy exactly; numeric
  Gaussians unchanged, undecidable-sign coefficients left unevaluated.
- **Regression test:** `INT-GAUSS-PARAM-1` in
  `tests/integrals/integrate_test.cpp`.

### INT-SQRTQUAD-PARAM-1 — `∫1/√(x²+a²)` unevaluated for symbolic coefficients
- **Problem:** `∫1/√(x²+a²)`, `∫1/√(a²−x²)`, `∫1/√(x²+a)` came back unevaluated
  for symbolic positive coefficients, even though `try_sqrt_quadratic`'s branches
  already use `is_positive`/`is_negative` (which handle symbolic) — a leftover
  rational-only gate blocked them.
- **Fix:** in `src/integrals/integrate.cpp`, removed the `is_rational(a)/`
  `is_rational(c)` gate in `try_sqrt_quadratic`'s pure-quadratic case; the
  sign-gated branches below decide the asinh / asin / log form. Combined with the
  `MUL-POS-1` fix (`is_positive(a²) = true`), symbolic positive coefficients now
  close.
- **Verified:** `∫1/√(x²+a²) = asinh(x/a)`, `∫1/√(a²−x²) = asin(x/a)`,
  `∫1/√(x²+a) = asinh(x/√a)`, `∫√(a²−x²) = (a²·asin(x/a) + x·√(a²−x²))/2` — the
  reciprocal forms match SymPy exactly. Undecidable-sign coefficients fall
  through unevaluated.
- **Regression test:** `INT-SQRTQUAD-PARAM-1` in
  `tests/integrals/integrate_test.cpp`.

### INT-ARCTAN-PARAM-1 / MUL-POS-1 — `∫1/(x²+a²)` unevaluated; `is_positive(4·a²)` unknown
- **Problem:** `∫1/(x²+a²)` (and `∫1/(x²+a)`, `∫1/(ax²+b)`) came back unevaluated
  for symbolic positive coefficients — `try_arctan_quadratic` required *rational*
  coefficients. Relaxing that exposed a second bug: `is_positive(4·a²)` returned
  *unknown* even for `a > 0`, although `is_positive(4·a)` and `is_positive(a²)`
  were both `true`.
- **Fix:**
  - `src/core/mul.cpp`: the `Positive`/`Negative` handlers now classify each
    factor via its own `Positive`/`Negative` (both imply nonzero), instead of
    requiring `Negative` plus a separate `Nonzero` gate. A factor like `a²`
    (positive, but with unknown `Nonzero`) now counts correctly, so
    `is_positive(4·a²) = true`.
  - `src/integrals/integrate.cpp`: `try_arctan_quadratic` accepts symbolic
    coefficients, firing the arctan branch only when the discriminant is
    *provably positive* (matching SymPy under positivity assumptions). The
    `disc = 0` and log branches stay restricted to rational coefficients.
- **Verified:** `∫1/(x²+a²) = atan(x/a)/a`, `∫1/(x²+a) = atan(x/√a)/√a`,
  `∫1/(ax²+b) = atan(x√(a/b))/√(ab)` (all for positive parameters, diff-back
  verified at concrete values); numeric quadratics are unchanged, and a generic
  (unsigned) parameter is conservatively left unevaluated.
- **Regression test:** `INT-ARCTAN-PARAM-1` in
  `tests/integrals/integrate_test.cpp`.

### INT-DEF-2 / LIMIT-LOG-1 — `∫₀^∞ 1/(1+x⁴) = nan` (log/atan antiderivative at ∞)
- **Problem:** `∫₀^∞ 1/(1+x⁴)` returned `nan` instead of `π√2/4`. Its
  antiderivative has `log(A) − log(B)` and `atan(arg)` terms; at the upper limit
  the logs gave `∞ − ∞` and the `atan` arguments stayed unevaluated. Two root
  causes:
  1. **Infinity arithmetic:** `oo + √2` did not collapse to `oo` — the `Add`
     infinity pre-pass only absorbed numeric *literals*, not closed real
     constants like `√2` or `π`. So `atan(½·(2x+√2)·√2)|_{x=∞}` kept an
     unevaluated `oo + √2` and never reached `atan(∞) = π/2`.
  2. **Limit engine:** no log-continuity or log-combination at `∞`, so
     `limit(log(x+1) − log(x))` was `nan` instead of `0`.
- **Fix:**
  - `src/core/add.cpp`: the `±∞` pre-pass now absorbs any finite *real constant*
    (`is_number` or no free symbols + `is_real`), so `oo + √2 = oo`,
    `oo + π = oo`; `oo + x` (symbolic) is still kept.
  - `src/calculus/limit.cpp`: added `try_log_limit` — log-continuity
    (`limit(log g) = log(lim g)`), `∞ − ∞` log-combination (factor a common κ so
    `Σ cᵢ·log gᵢ = κ·log(∏ gᵢ^(cᵢ/κ))` with a single rational argument), and
    atan-continuity (`limit(atan g) = atan(lim g)`), applied before direct
    substitution.
- **Verified:** `∫₀^∞ 1/(1+x⁴) = π√2/4`, `∫₀^∞ 1/(x⁴+x²+1) = π√3/6`,
  `∫₁^∞ 1/(x(x+1)) = log 2`; `limit(log(x+1) − log x) = 0`,
  `limit(log(x²+x+1) − log(x²−x+1)) = 0`, `limit(atan(2x+1)) = π/2` — all match
  SymPy.
- **Regression tests:** `INT-DEF-2` in `tests/integrals/integrate_test.cpp` and
  `LIMIT-LOG-1` in `tests/calculus/series_limit_test.cpp`.

### SIMP-CXDIV-1 — `simplify((1+I)/(1-I))` left the complex quotient unreduced
- **Problem:** `simplify((1+I)/(1−I))` returned `(1+I)·(1−I)⁻¹` instead of `I`;
  `simplify(1/(1+I))` stayed `(1+I)⁻¹` instead of `1/2 − I/2`. Complex *products*
  expand, but `simplify` never rationalized a complex denominator. (The previous
  iteration fixed `re`/`im`; this fixes the bare `simplify`.)
- **Fix:** exposed `rationalize_complex` (`include/sympp/functions/miscellaneous.hpp`)
  and applied it in `simplify` right after the initial `expand`. Since a
  rationalized quotient can be *larger* than the input (`1/(1+I)` → `1/2 − I/2`),
  the anti-bloat guard now exempts the case where a complex denominator was
  removed — mirroring the existing surd-denominator exemption.
- **Verified:** `(1+I)/(1−I) → I`, `1/(1+I) → 1/2 − I/2`,
  `(2+3I)/(1+I) → 5/2 + I/2`, `I/(2−I) → −1/5 + 2I/5`,
  `(3+4I)/(1+2I) → 11/5 − 2I/5` — all match SymPy; real rational functions
  (`(x²−1)/(x−1) → x+1`) and the anti-bloat guarantee (`(x+1)³` stays factored)
  are unchanged, and a symbolic complex denominator (`1/(x−I)`) is conservatively
  left alone.
- **Regression test:** `SIMP-CXDIV-1` in `tests/simplify/simplify_test.cpp`
  (`[5][simplify][complex][oracle][regression]`).

### REIM-CXDIV-1 — `re`/`im` of an expression with a complex denominator stayed unevaluated
- **Problem:** `re((1+I)/(1−I))` and `im((1+I)/(1−I))` returned an unevaluated
  `re(...)`/`im(...)` instead of `0` and `1`. Complex *products* already expand
  (`(1+I)² = 2I`), but a complex *denominator* `(a+bI)⁻¹` was never rationalized,
  so `re`/`im` couldn't reach the `a+bI` form they already handle.
- **Fix:** added `rationalize_complex` in `src/functions/miscellaneous.cpp`. It
  rewrites every `Pow(d, −m)` whose base `d` carries the imaginary unit and whose
  `|d|² = d·conj(d)` is provably real, as `conj(d)^m/|d|^{2m}` —
  i.e. `1/(a+bI) = (a−bI)/(a²+b²)`. `re`/`im` apply it (then `expand`) to their
  argument and re-enter on the resulting `a+bI` form; the value then folds at
  construction, so `re((1+I)/(1−I))` evaluates to `0` directly.
- **Verified:** `re((1+I)/(1−I)) = 0`, `im = 1`, `1/(1+I) → re 1/2, im −1/2`,
  `(2+3I)/(1+I) → re 5/2, im 1/2`, `I/(2−I) → re −1/5` — all match SymPy;
  symbolic/real arguments (`re(x+Iy) = re(x)−im(y)`) are unchanged.
- **Regression test:** `REIM-CXDIV-1` in
  `tests/functions/miscellaneous_test.cpp` (`[3h][complex][oracle][regression]`).
- **Scope:** numeric (provably-real `|d|²`) complex denominators. A symbolic
  denominator whose `|d|²` stays complex is left untouched.

### DSOLVE-UNIFIED-1 — no single-entry `dsolve(eq, y, x)` (only per-method solvers)
- **Problem:** SymPP exposed `dsolve_first_order`, `dsolve_constant_coeff`,
  `dsolve_cauchy_euler`, … but had no unified `dsolve(eq)` like SymPy's — the
  caller had to know the ODE class and the right signature in advance.
- **Fix:** added `dsolve(eq, y, x)` in `src/ode/dsolve.cpp`. It finds the order
  from the highest derivative of `y` present, delegates a first-order ODE to
  `dsolve_first_order`, and for a linear higher-order ODE linearizes (each
  `y^(k)` → a fresh symbol), extracts the coefficients `aₖ` and rhs `g(x)`, and
  dispatches: constant `aₖ` → `dsolve_constant_coeff` (homogeneous) /
  `dsolve_constant_coeff_nonhomogeneous` (order 2); `aₖ = cₖ·xᵏ` →
  `dsolve_cauchy_euler`. A nonlinear or unrecognized ODE returns an unevaluated
  `Dsolve(...)` marker.
- **Verified:** every general solution substitutes back to an ODE residual of 0
  — `y'=y`, `y'+y=x`, `y''+y=0`, `y''−3y'+2y=0` (distinct roots),
  `y''−2y'+y=0` (repeated root), `y''+y=x` (nonhomogeneous),
  `x²y''−2y=0` (Cauchy-Euler), `y'''−y'=0` (third order).
- **Regression test:** `DSOLVE-UNIFIED-1` in `tests/ode/dsolve_test.cpp`
  (`[11][dsolve][oracle][regression]`).
- **Scope:** linear ODEs with constant or `cₖ·xᵏ` coefficients (any order
  homogeneous; order 2 nonhomogeneous). General variable-coefficient linear and
  nonlinear higher-order ODEs are still per-method / unevaluated.

### INT-EXP-SUB-1 — `∫1/(eˣ+e⁻ˣ)` and other eˣ-rational integrals were unevaluated
- **Problem:** `∫1/(eˣ+e⁻ˣ)`, `∫eˣ/(e²ˣ+1)`, `∫e²ˣ/(1+eˣ)`, `∫1/(eˣ+e²ˣ)` came
  back unevaluated. The heurisch substitution `subs(eˣ → u)` does not catch
  `e²ˣ` or `e⁻ˣ` — those are distinct nodes (`exp(2x)`, `exp(−x)`), not powers
  of `exp(x)` — so the substituted integrand still depended on `x` and bailed.
- **Fix:** added `try_exp_substitution` in `src/integrals/integrate.cpp`. It maps
  every `exp(k·x+d)` (integer `k`) to `e^d·uᵏ` with `u = eˣ`, and `dx = du/u`,
  turning the integrand into a rational function of `u` that `try_rational` /
  `integrate` closes; it back-substitutes `u = eˣ`.
- **Verified:** `∫1/(eˣ+e⁻ˣ) = atan(eˣ)`, `∫eˣ/(e²ˣ+1) = atan(eˣ)`,
  `∫e²ˣ/(1+eˣ) = eˣ−log(1+eˣ)`, `∫1/(eˣ+e²ˣ) = −e⁻ˣ−x+log(eˣ+1)`,
  `∫1/(eˣ+4e⁻ˣ) = atan(eˣ/2)/2` — all differentiate back to the integrand
  (and the headline matches SymPy's `atan(eˣ)`). The previously-working
  `1/(eˣ+1)` family is unchanged.
- **Regression test:** `INT-EXP-SUB-1` in `tests/integrals/integrate_test.cpp`
  (`[7][integrate][oracle][regression]`).
- **Scope:** integrands rational in `eˣ` with integer exponent multiples. A
  fractional rate (`e^(x/2)`) would need `u = e^(x/2)` and is left unhandled.

### SOLVE-EXPBASE-SUM-1 — sums of constant-base exponentials returned `[]`
- **Problem:** `solve(2^x − 3^x)`, `solve(2^(2x) − 5·2^x + 4)`,
  `solve(2^(x+1) − 8)`, `solve(2^x·3^x − 6)` returned `[]`. The existing
  constant-base solver handled only a single `a^x = c`; sums of several
  exponential terms (possibly with different bases) were unhandled.
- **Fix:** added `solve_const_base_exp_sum` in `src/solvers/solve.cpp`. Each
  term reduces to `coeff·exp(rate·x)` with `rate = Σ pⱼ·log(aⱼ)`. After
  combining equal rates: **(A)** when every rate is an integer multiple of a
  common `r₀`, substitute `u = exp(r₀·x)` → a polynomial in `u`
  (`2^(2x)−5·2^x+4 → u²−5u+4`); **(B)** with two incommensurate rates,
  `d₁·exp(r₁x)+d₂·exp(r₂x)=0 ⇒ x = log(−d₂/d₁)/(r₁−r₂)` when `−d₂/d₁ > 0`. Only
  real roots are kept (positive `u`). Pure `exp(…)` equations are deferred to
  `solve_exp_sum` so its complex (period-`2πi`) roots survive.
- **Verified:** `2^x−3^x → 0`, `5^x−2^x → 0`, `2^(2x)−5·2^x+4 → {0,2}`,
  `2^(x+1)−8 → 2`, `2^x·3^x−6 → 1`, `4^x−2^(x+1) → 1`, `9^x−3^(x+1) → 1`,
  `4^x−2 → 1/2`, `2^(2x)−8 → 3/2` — all match SymPy.
- **Regression test:** `SOLVE-EXPBASE-SUM-1` in `tests/solvers/solve_test.cpp`;
  the `SOLVE-EXPBASE-1` "stays unsolved" assertions for `4^x−2` and `2^(2x)−8`
  were updated (they now solve).

### SOLVE-LOGSUM-1 — `solve(log(x)+log(x−1))` returned `[]`
- **Problem:** equations with a *sum* of logarithms — `log(x)+log(x−1)`,
  `log(x)+log(x+1)−log(6)`, `2·log(x)−log(x+2)` — returned `[]`. The existing
  log solver handles only a single log atom; a sum of several is not a
  polynomial in one atom.
- **Fix:** added `solve_log_sum` in `src/solvers/solve.cpp`. It recognizes
  `Σ cᵢ·log(gᵢ(x)) + K` (cᵢ, K var-free), combines via
  `log(∏ gᵢ^cᵢ) = −K ⇒ ∏ gᵢ^cᵢ = exp(−K)`, solves that recursively, and keeps
  only roots in the log domain (every `gᵢ(root) > 0`). The domain filter uses a
  numeric sign from `evalf`, since `is_positive` cannot judge an irrational like
  `(1+√5)/2`.
- **Verified:** `log(x)+log(x−1)=0 → (1+√5)/2` (the negative root dropped),
  `log(x)+log(x+1)=log(6) → 2`, `2log(x)−log(x+2)=0 → 2`,
  `log(x+1)+log(x−1)=0 → √2`, `log(x)−log(x−1)=1 → e/(e−1)` — all match SymPy;
  single-log equations are unchanged.
- **Regression test:** `SOLVE-LOGSUM-1` in `tests/solvers/solve_test.cpp`
  (`[10][solve][oracle][regression]`).
- **Scope:** sums of `cᵢ·log(gᵢ)` with var-free coefficients. A log with the
  variable also outside a log, or symbolic coefficients, is left to other paths.

### SOLVE-ABS-1 — `solve(|x−1|−2)` returned `[]`, and `|g|=c<0` gave spurious roots
- **Problem:** `solve(abs(x−1)−2)` returned `[]` instead of `{3, −1}`.
  `solveset` correctly produced `{3} ∪ {−1}`, but `solve` only extracted roots
  from a single `FiniteSet`, not a **Union** of finite sets. (`abs(x)−3` worked
  only because its solveset is one FiniteSet.) Exposing the Union also revealed
  a soundness bug: `|g| = c` with a negative `c` (e.g. `|x+1|+2 = 0`) returned
  spurious roots, since the inverse never checked `c ≥ 0`.
- **Fix:** in `src/solvers/solve.cpp`, the solveset-extraction step now flattens
  a `FiniteSet`, the empty set, or a `Union` of finite sets into the root list
  (deduplicated); anything with a non-finite component is left empty. The
  solveset `Abs` inverse now returns the empty set when `c` is a concrete
  negative.
- **Verified:** `|x−1|=2 → {3,−1}`, `|2x−1|=5 → {3,−2}`, `|x²−1|=3 → {2,−2}`,
  `|x|=0 → {0}`, and `|x+1|+2`, `|x|+5 → ∅` — all match SymPy (real domain).
- **Regression test:** `SOLVE-ABS-1` in `tests/solvers/solve_test.cpp`
  (`[10][solve][oracle][regression]`).
- **Scope:** `|affine or polynomial| = const`. An absolute value with the
  variable also outside (`|x−1| = x`) or a coefficient on the abs (`2|x| = 6`)
  is still unhandled.

### INT-RECIP-SUB-1 — `∫1/(xⁿ√(a x²+c))` was unevaluated
- **Problem:** `∫1/(x√(x²+1))`, `∫1/(x²√(x²+1))`, `∫1/(x√(x²+4))` came back
  unevaluated. These need the reciprocal substitution `x = 1/u`, which the
  engine lacked.
- **Fix:** added `try_reciprocal_substitution` in
  `src/integrals/integrate.cpp`. It gates on an integrand with a negative power
  of the variable AND a half-integer power of a degree-2 polynomial, substitutes
  `x = 1/u` (`dx = −u⁻² du`), and — since SymPP can't pull a power out of a
  radical on its own — does the targeted rewrite `(a·u⁻²+c)^e =
  u^(−2e)·(a+c·u²)^e`, leaving an ordinary `√(quadratic)` integral that the
  existing machinery closes. Back-substitutes `u = 1/x`. The integrand is
  `expand`ed first so `(x·√(…))⁻¹` flattens to `x⁻¹·(…)^(−1/2)` for the gate.
- **Verified:** `∫1/(x√(x²+1)) = −asinh(1/x)`, `∫1/(x√(x²+4)) = −asinh(2/x)/2`,
  `∫1/(x√(1+9x²)) = −asinh(1/(3x))/3`, and the `x²`/`x³` denominator cases — all
  match SymPy (diff-back verified on the `x>0` principal branch, the same
  convention SymPy's answers use).
- **Regression test:** `INT-RECIP-SUB-1` in
  `tests/integrals/integrate_test.cpp` (`[7][integrate][oracle][regression]`).
- **Scope:** `√(a·x²+c)` (no linear term). The `√(x²−1)`/`√(1−x²)` variants give
  branch-dependent Piecewise answers in SymPy and are left to the cleaner paths.

### SUM-POLYEXPAND-1 — `Σ k·(k+1)` and other product summands stayed unevaluated
- **Problem:** `summation(k·(k+1))`, `(k+1)²`, `(2k+1)(k−1)` returned an
  unevaluated `Sum(...)`, even though the expanded `Σ(k²+k)` summed fine via
  Faulhaber. A product or power summand isn't matched by the closed-form
  branches, and the constant-pull only fires when there's a var-free factor.
- **Fix:** in `src/calculus/summation.cpp`, before the closed-form dispatch,
  expand a `Mul`/`Pow` summand and — when expansion produces an `Add` — recurse,
  so linearity splits it into individually-summable terms (monomials `kᵖ`, or
  poly·geometric). This also picks up mixed forms like `(k+1)·2ᵏ`.
- **Verified:** `Σ k(k+1) = n(n+1)(n+2)/3`, `Σ k(k−1) = n(n−1)(n+1)/3`,
  `Σ (k+1)² = n(2n²+9n+13)/6`, `Σ (2k+1)(k−1) = n(n−1)(4n+7)/6` — all match
  SymPy; pure geometric/exponential summands (`2ᵏ`, `k·2ᵏ`) are unaffected
  (they don't expand to an `Add`).
- **Regression test:** `SUM-POLYEXPAND-1` in
  `tests/calculus/series_limit_test.cpp` (`[6][summation][oracle][regression]`).

### POLYOP-2 — `resultant` and `discriminant` parsed to unevaluated nodes
- **Problem:** `resultant(x²−1, x−1)` and `discriminant(x²+1)` came back as
  opaque function nodes, even though `resultant(p, q, var)` and
  `discriminant(p, var)` already existed and were tested — they just required an
  explicit variable and weren't registered with the parser.
- **Fix:** added parser-facing `resultant(p, q)` (two-arg) and
  `discriminant(p)` (one-arg) wrappers in `src/polys/poly.cpp` that infer the
  variable from the single free symbol (reusing `inferred_var`), and registered
  them. Same pattern as `POLYOP-1`.
- **Verified:** `discriminant(x²+2x+1) = 0`, `discriminant(x²−5x+6) = 1`,
  `discriminant(x²+1) = −4`, `discriminant(x³−3x+1) = 81`,
  `resultant(x²−1, x−1) = 0`, `resultant(x²+1, x−2) = 5`, and the sign
  convention `resultant(x−1, x−2) = −1` vs `resultant(x−2, x−1) = 1` — all match
  SymPy.
- **Regression test:** `POLYOP-2` in `tests/polys/poly_test.cpp`
  (`[4][poly][regression]`).

### POLYOP-1 — `degree`, `quo`, `rem`, `cancel` parsed to unevaluated nodes
- **Problem:** `degree(x³+2x)`, `quo(x²−1, x−1)`, `rem(x², x−1)` and the
  one-argument `cancel((x²−1)/(x−1))` came back as opaque function nodes. The
  `cancel(expr, var)` C++ function existed but needed an explicit variable, and
  `degree`/`quo`/`rem` were not implemented or registered with the parser.
- **Fix:** added parser-facing wrappers in `src/polys/poly.cpp` that infer the
  polynomial variable from the single free symbol (`inferred_var`), then call
  the `Poly` primitives: `degree → Poly::degree`, `quo`/`rem →
  `Poly::divmod`, and a 1-argument `cancel` over the existing `cancel(expr,
  var)`. Each falls back to an unevaluated node when the argument is not a
  univariate polynomial expression. Registered `cancel`, `degree` (one-arg) and
  `quo`, `rem` (two-arg) in the parser.
- **Verified:** `degree(x³+2x) = 3`, `degree(5) = 0`, `quo(x²−1, x−1) = x+1`,
  `quo(x³−1, x−1) = x²+x+1`, `rem(x², x−1) = 1`,
  `cancel((x²−1)/(x−1)) = x+1` — all match SymPy.
- **Regression test:** `POLYOP-1` in `tests/polys/poly_test.cpp`
  (`[4][poly][oracle][regression]`).
- **Scope:** univariate. `degree(0) = −∞` and `degree(c≠0) = 0` for constants,
  matching SymPy.

### LCM-POLY-1 — `lcm` of polynomials stayed unevaluated
- **Problem:** `lcm(x²−1, x−1)` returned an unevaluated `lcm(...)` node instead
  of `x²−1`. Like `gcd`, the `lcm` function only handled two integers.
- **Fix:** in `src/functions/combinatorial.cpp`, `lcm(a, b)` now computes the
  univariate polynomial LCM as `a·b / gcd(a, b)` (reusing the polynomial gcd
  from `GCD-POLY-1`) via exact `Poly` division. The division restores the right
  content automatically.
- **Verified:** `lcm(x²−1, x−1) = x²−1`, `lcm(x−1, x+1) = x²−1`,
  `lcm(2x−2, 3x−3) = 6x−6`, `lcm(x, x²) = x²`,
  `lcm(x²−1, x²+2x+1) = x³+x²−x−1` — all match SymPy.
- **Regression test:** `LCM-POLY-1` in `tests/functions/combinatorial_test.cpp`
  (`[3i][lcm][oracle][regression]`).
- **Note:** `lcm(x, n)` now eagerly evaluates to `n·x` (matching SymPy), so the
  two integer-lcm tests that relied on the old lazy node were updated. As with
  gcd, multivariate LCM stays an unevaluated node (the `Poly` class is
  univariate).

### GCD-POLY-1 — `gcd` of polynomials stayed unevaluated
- **Problem:** `gcd(x²−1, x−1)` returned an unevaluated `gcd(...)` node instead
  of `x−1`. The `gcd` function only handled two integers, even though the `Poly`
  class already provides a Euclidean polynomial GCD.
- **Fix:** in `src/functions/combinatorial.cpp`, `gcd(a, b)` now detects a common
  single variable (via `free_symbols`), builds `Poly`s, and computes the GCD.
  SymPy's convention is the **primitive integer** gcd (integer coefficients,
  content 1, positive leading) scaled by the gcd of the integer contents, so the
  monic `Poly` GCD is re-primitivized (`gcd_to_primitive`): clear denominators,
  divide by the integer content, then multiply by `gcd(content a, content b)`.
- **Verified:** `gcd(x²−1, x−1) = x−1`, `gcd(2x²−2, 2x−2) = 2x−2`,
  `gcd(6x²+11x+3, 2x²−x−6) = 2x+3` (primitive, not the monic `x+3/2`),
  `gcd(x²+1, x−1) = 1`, `gcd(x²−1, 2) = 1`, `gcd(x, 18) = 1` — all match SymPy.
- **Regression test:** `GCD-POLY-1` in `tests/functions/combinatorial_test.cpp`
  (`[3i][gcd][oracle][regression]`).
- **Note:** `gcd(x, n)` now eagerly evaluates to `1` (x and a constant are
  coprime over ℚ[x]), matching SymPy; the parse-round-trip test that relied on
  the old lazy node was updated. Multivariate GCD (`gcd(x²−y², x−y)`) remains an
  unevaluated node — the `Poly` class is univariate.

### LIMIT-CONJUGATE-1 — `x − √(x²+1)` and radical ∞−∞ limits returned nan
- **Problem:** `limit(x − √(x²+1), ∞)` returned `nan` instead of `0`; likewise
  `x − √(x²−1)`, `√(x+1) − √x`. Direct substitution gives the indeterminate
  `∞ − ∞`, and the existing polynomial / L'Hôpital paths don't handle radicals.
- **Fix:** added `try_conjugate_difference` in `src/calculus/limit.cpp`. For a
  two-term sum `t₁ + t₂` containing a radical whose limit is `∞ − ∞`, it
  rationalizes via the conjugate: `t₁ + t₂ = (t₁² − t₂²)/(t₁ − t₂)`. Squaring
  clears the radical from the numerator, and the resulting ratio resolves. A
  key subtlety: the ratio is passed to `limit` **unsimplified**, because
  `simplify` would rationalize the denominator straight back to the original
  `∞ − ∞` form and loop (`limit` substitutes before simplifying, so the pole
  collapses first).
- **Verified:** `x − √(x²+1) → 0`, `x − √(x²−1) → 0`, `√(x+1) − √x → 0`;
  the non-indeterminate `x + √(x²+1) → ∞` is unaffected.
- **Regression test:** `LIMIT-CONJUGATE-1` in
  `tests/calculus/series_limit_test.cpp` (`[6][limit][infinity][oracle][regression]`).
- **Scope:** the conjugate resolves cases where squaring leaves a *constant* (or
  lower-degree) numerator. `√(x²+x) − x → 1/2` is still open — its conjugate
  leaves an `∞/∞`-with-radical ratio that needs a leading-term asymptotic
  expansion (factoring the dominant power out of the radical). The log-ratio
  `log x / log(2x) → 1` is also still open (different root cause).

### LIMIT-POWFORM-1 — `(1+x)^(1/x)` and other 1^∞ limits returned 1 instead of e
- **Problem:** `limit((1+x)^(1/x), x, 0)` returned `1` instead of `e` — the
  textbook definition of e. Likewise `(1+2x)^(1/x) → 1` (should be `e²`),
  `cos(x)^(1/x²) → 1` (should be `e^(−1/2)`), `(1−x)^(1/x) → 1` (`e⁻¹`). At a
  finite target, direct substitution evaluates the exponent `1/x` to `zoo` and
  collapses `pow(1, zoo)` to `1` *before* the `1^∞` indeterminate handler runs,
  so the indeterminacy was lost. (The same forms at `∞` already worked, because
  `pow(1, ∞)` surfaced as `nan` there.)
- **Fix:** in `src/calculus/limit.cpp`, call `try_power_form` for a `Pow`
  expression *before* the direct-substitution step. It resolves the genuine
  indeterminate forms `1^∞`, `0^0`, `∞^0` via `exp(lim exponent·log base)` and
  returns `nullopt` for any determinate power, so ordinary powers
  (`(1+x)²`, `2^x`, `x^x`) are unaffected.
- **Verified:** `(1+x)^(1/x) → e`, `(1+2x)^(1/x) → e²`, `(1+x)^(2/x) → e²`,
  `(1−x)^(1/x) → e⁻¹`, `(1+3x)^(2/x) → e⁶`, `cos(x)^(1/x²) → e^(−1/2)` — all
  match SymPy; determinate powers and the `∞`-target cases are unchanged.
- **Regression test:** `LIMIT-POWFORM-1` in
  `tests/calculus/series_limit_test.cpp` (`[6][limit][oracle][regression]`).
- **Note:** correctness bug (confidently wrong answers). Other limit gaps
  surfaced in the same survey — `x − √(x²+1) → 0` and `log x / log(2x) → 1`
  still return `nan` — remain open (different root causes).

### SERIES-LAURENT-1 — functions with a pole at 0 had no series expansion
- **Problem:** `series(cot(x))`, `csc(x)`, `coth(x)`, `csch(x)`, `csc(x)²`,
  `1/(eˣ−1)` all returned the input unexpanded. The series engine was a pure
  Taylor expansion: at a pole the leading coefficient is non-finite, so it gave
  up. (Even `x·cot(x)`, which is analytic, failed — the Taylor path's
  higher derivatives hit ∞−∞ forms the limit engine could not resolve.)
- **Fix:** rewrote `src/calculus/series.cpp` around a **power-series division**
  Laurent path. When the ordinary Taylor expansion fails, the engine rewrites
  reciprocal trig/hyperbolic functions to sin/cos ratios
  (`cot→cos/sin`, `csc→1/sin`, …), splits the result into numerator `N` and
  denominator `D`, Taylor-expands both (analytic), and divides the power series:
  `f = x^(v_N − v_D)·(Ñ/D̃)` with `Ñ(0), D̃(0) ≠ 0`. This yields the Laurent
  series directly, including negative powers, without differentiating the
  singular function. Genuine singularities (`log x`) still return unexpanded.
- **Verified:** `cot(x) = 1/x − x/3 − x³/45 − …`,
  `csc(x) = 1/x + x/6 + 7x³/360 + …`, `coth`, `csch`, `csc²(x) = 1/x² + 1/3 + …`,
  `1/(eˣ−1) = 1/x − 1/2 + x/12 − …`, and `x·cot(x) = 1 − x²/3 − …` — all match
  SymPy; analytic functions (`exp`, `sin`, `1/(1−x)`) and `log x`, `1/x`, `1/x²`
  are unchanged.
- **Regression test:** `SERIES-LAURENT-1` in
  `tests/calculus/series_limit_test.cpp` (`[6][series][oracle][regression]`).
- **Scope:** Laurent expansion at `x0 = 0`. A pole at a non-zero point would
  need the same division after shifting the expansion variable.

### LIMIT-RECIPTRIG-1 — limits of cot/csc/sec (and hyperbolic) returned nan
- **Problem:** `limit(x·cot(x), 0)` returned `nan` instead of `1`; likewise
  `cot(x)·sin(x)`, `x·csc(x)`, `x·coth(x)`, `x²·csc²(x)`. The limit machinery
  (direct substitution, L'Hôpital) understands sin/cos but treats the
  reciprocal functions cot/csc/sec/coth/csch/sech as opaque, so any `0·∞` form
  built from them failed.
- **Fix:** added `rewrite_reciprocal_trig` in `src/calculus/limit.cpp`, applied
  at the top of `limit_impl`: it rewrites `cot→cos/sin`, `csc→1/sin`,
  `sec→1/cos`, `coth→cosh/sinh`, `csch→1/sinh`, `sech→1/cosh` and retries. The
  rewrite is exact, so the limit is unchanged; the sin/cos form is one the
  L'Hôpital path resolves.
- **Verified:** `x·cot(x) → 1`, `cot(x)·sin(x) → 1`, `x·csc(x) → 1`,
  `x·coth(x) → 1`, `x²·csc²(x) → 1`, `tan(x)·cot(x) → 1`,
  `(cos x − 1)·csc(x) → 0` — all match SymPy. (`limit(cot(x), 0)` is `zoo`, the
  correct two-sided value; SymPy's default one-sided gives `oo`.)
- **Regression test:** `LIMIT-RECIPTRIG-1` in
  `tests/calculus/series_limit_test.cpp` (`[6][limit][oracle][regression]`).
- **Note:** this also unblocks part of the still-open Laurent-series gap
  (`series(cot(x)) = 1/x − x/3 − …`), which additionally needs pole handling in
  the series engine.

### SOLVE-EQ-1 — `solve(Eq(lhs, rhs))` and relational parsing returned `[]`
- **Problem:** `solve(Eq(x**2, 4))` returned `[]` instead of `{2, −2}`. Two
  causes: (1) the parser built `Eq(a, b)` (and `Ne`/`Lt`/`Le`/`Gt`/`Ge`) as an
  opaque user-function node rather than a `Relational`, and (2) `solve` had no
  branch to reduce an equation to `lhs − rhs = 0`.
- **Fix:**
  - registered `Eq`, `Ne`, `Lt`, `Le`, `Gt`, `Ge` in the parser's two-argument
    table (`src/parsing/parser.cpp`), so they build proper `Relational` nodes;
  - in `src/solvers/solve.cpp`, `solve` now reduces a `Relational` of kind `Eq`
    to `solve(lhs − rhs, var)` (matching SymPy's `solve(Eq(...))`). Inequalities
    describe a region, not a discrete root list, so they are not forced into the
    vector API.
- **Verified:** `Eq(x², 4) → {2, −2}`, `Eq(x³, x) → {0, 1, −1}`,
  `Eq(sin x, 1/2) → {π/6, 5π/6}`, `Eq(eˣ, 3) → {log 3}`, `Eq(2x+1, 5) → {2}`,
  and the parsed-string forms — all match SymPy; `Eq(x, x)` still evaluates to
  `True`.
- **Regression test:** `SOLVE-EQ-1` in `tests/solvers/solve_test.cpp`
  (`[10][solve][oracle][regression]`).

### SUM-EXP-2 — polynomial × exponential series Σ P(k)·rᵏ/k! stayed unevaluated
- **Problem:** `Σ k/k!`, `Σ k²/k!`, `Σ (2k+3)/k!`, `Σ k·xᵏ/k!` came back
  unevaluated. `SUM-EXP-1` closed only a bare `rᵏ/k!`; a polynomial numerator
  `P(k)` was an unrecognized factor and bailed.
- **Fix:** generalized `sum_exponential_series` in
  `src/calculus/summation.cpp` to collect a polynomial numerator `P(var)` and
  fold it through the **falling-factorial basis**: since
  `Σ_{k≥0} k^{(m)}·rᵏ/k! = rᵐ·eʳ`, writing `P = Σ_m c_m·k^{(m)}` gives
  `Σ P(k)·rᵏ/k! = (Σ_m c_m·rᵐ)·eʳ = Q(r)·eʳ`. The transform
  (`exp_series_poly_transform`) extracts the monic falling factorials top-down
  (a triangular solve, no Stirling-number table). Head terms for `lo > 0` use
  the full `P(k)·rᵏ/k!`.
- **Verified:** `Σ k/k! = e`, `Σ k²/k! = 2e`, `Σ k³/k! = 5e`,
  `Σ (k+1)/k! = 2e`, `Σ (2k+3)/k! = 5e`, `Σ k·xᵏ/k! = x·eˣ`,
  `Σ_{k≥1} k/k! = e` — all match SymPy; the bare-`rᵏ/k!` cases are unchanged.
- **Regression test:** extended `SUM-EXP-1` in
  `tests/calculus/series_limit_test.cpp`; the `SUM-3` unrecognised-sum stand-in
  moved to `Σ 1/(k!+1)` (no elementary closed form), since `Σ k/k!` now closes.

### SUM-EXP-1 — exponential series Σ rᵏ/k! stayed unevaluated
- **Problem:** `Σ_{k=0}^∞ 1/k!`, `Σ x^k/k!`, `Σ 2^k/k!`, `Σ (−1)^k/k!` all came
  back as an unevaluated `Sum(...)`. SymPP already closed the convergent
  p-series (`Σ1/k²=π²/6`) and geometric sums, but not the factorial/exponential
  family. (Note: these were reachable only through the 4-argument `summation`
  with an `∞` bound — the CLI probe hides the bounds, which made it look like
  even `Σ1/k²` failed when it did not.)
- **Fix:** added `sum_exponential_series` in `src/calculus/summation.cpp`.
  It recognizes `c · (∏ baseᵢ^(aᵢ·k + bᵢ)) · k!^(−1)` at an `∞` upper bound:
  each base contributes `baseᵢ^{aᵢ}` to the rate `r` and `baseᵢ^{bᵢ}` to the
  constant `c`, giving `Σ_{k=0}^∞ c·rᵏ/k! = c·eʳ`. For a lower bound `lo > 0`
  the omitted head `Σ_{k=0}^{lo−1} rᵏ/k!` is subtracted. The series is entire,
  so no convergence test is required.
- **Verified:** `Σ1/k! = e`, `Σ_{k≥1}1/k! = e−1`, `Σx^k/k! = e^x`,
  `Σ2^k/k! = e²`, `Σ(−1)^k/k! = e⁻¹`, `Σ1/(2^k·k!) = e^(1/2)`, `Σ3/k! = 3e`,
  `Σ_{k≥2}x^k/k! = e^x − x − 1` — all match SymPy.
- **Regression test:** `SUM-EXP-1` in `tests/calculus/series_limit_test.cpp`
  (`[6][summation][oracle][regression]`); the `SUM-3` unrecognised-sum test now
  uses `Σ k/k!` as its stand-in since `Σ1/k!` closes.
- **Scope:** pure `rᵏ/k!`. A polynomial-times-`1/k!` numerator (`Σ k/k! = e`)
  needs an index-shift reduction and is still left unevaluated.

### LIMIT-GAMMA-1 — limits of gamma/factorial at ∞ returned wrong answers
- **Problem:** `limit(gamma(x+1)/gamma(x), ∞)` returned **`1`** (should be `∞` —
  the ratio *is* `x`), `exp(x)/gamma(x) → ∞` (should be `0`) and
  `gamma(x)/exp(x) → 0` (should be `∞`) were **inverted**, and bare `gamma(x)` /
  `factorial(x)` stayed unevaluated. Root cause: `limit` substitutes `x → ∞`
  *before* simplifying, producing `gamma(∞)/gamma(∞)`, which `simplify` then
  cancels to `1`; the engine had no model of gamma's super-exponential growth.
- **Fix:** three coordinated changes —
  - **(A)** `gamma(+∞) = +∞` and `factorial(+∞) = +∞` at the factories
    (`src/functions/combinatorial.cpp`);
  - **(B)** in `src/calculus/limit.cpp`, when the target is `+∞` and a
    gamma/factorial is present, normalize `factorial(u) → gamma(u+1)` and
    `simplify` first — this collapses `gamma(x+1)/gamma(x) → x`, after which the
    existing rational-at-∞ machinery gives `∞`;
  - **(C)** a conservative growth-rank rule (`gamma ≫ exp ≫ polynomial ≫ log`):
    when one factor strictly dominates, the limit is `∞` or `0`
    (`exp(x)/gamma(x) → 0`, `gamma(x)/exp(x) → ∞`, `2^x/x! → 0`). A genuine tie
    at the top rank (`gamma(2x)/gamma(x)`) is left unevaluated rather than
    guessed. A fallback rewrites any leftover factorial to gamma so the
    L'Hôpital path never differentiates `factorial` into a `Derivative` node.
- **Verified:** 15 cases match SymPy — bare gamma/factorial, integer-shift
  ratios, and the cross-class quotients. The only unresolved case,
  `gamma(2x)/gamma(x)`, now returns `nan` (honest) instead of a wrong value;
  it needs a Stirling comparison to fold to `∞`.
- **Regression test:** `LIMIT-GAMMA-1` in
  `tests/calculus/series_limit_test.cpp`
  (`[6][limit][infinity][gamma][oracle][regression]`).
- **Note:** this was a correctness bug (confidently wrong answers), not just a
  missing feature.

### LOG-BASE-1 — two-argument `log(x, b)` was an opaque node with no derivative
- **Problem:** `log(x, 2)` parsed to a generic user-function node, so it never
  simplified and `diff(log(x, 2), x)` came back as an unevaluated
  `Derivative(log(x, 2), x)` instead of `1/(x·log 2)`. SymPP had no two-argument
  `log` factory; the parser fell through to `function_symbol("log")`.
- **Fix:** added `log(arg, base)` in `src/functions/exponential.cpp` returning
  `log(arg)/log(base)` (matching SymPy), and registered `log` in the parser's
  two-argument table. Exact integer powers fold to the exponent
  (`log(8, 2) = 3`, `log(1024, 2) = 10`) via a divisibility loop; base `E`
  collapses to the natural log. Because the result is built from standard
  one-argument logs, `diff`, `simplify`, `series`, etc. handle it for free.
- **Verified:** `log(x, 2) → log(x)/log(2)`, `log(x, E) → log(x)`,
  `log(8, 2) → 3`, `log(100, 10) → 2`, `log(1024, 2) → 10`,
  `diff(log(x, 2)) → 1/(x·log 2)` — all match SymPy.
- **Regression test:** `LOG-BASE-1` in `tests/functions/exponential_test.cpp`
  (`[3c][log][oracle][regression]`).
- **Scope:** exact integer powers fold; a non-power integer ratio
  (`log(8, 4)`) stays the correct-but-unfolded `log(8)/log(4)` rather than
  SymPy's fully-reduced `3/2`.

### EXP-LOG-INVERSE-1 — `exp(log(x))` stayed unevaluated for a generic argument
- **Problem:** `exp(log(x))` returned `exp(log(x))` instead of `x`, and
  `exp(2·log(x)) → exp(2·log(x))` instead of `x²`, for a generic (unknown-sign)
  symbol. The factory gated `exp∘log` folding on `is_positive(arg) == true`,
  which is unknown for a bare symbol, so the simplification never fired —
  diverging from SymPy, which folds unconditionally. Sums such as
  `exp(log x + 1)` and `exp(log x + log y)` were likewise left intact.
- **Fix:** in `src/functions/exponential.cpp`,
  - `exp(log(w)) → w` is now unconditional: `exp(log(z)) = z` for every `z ≠ 0`
    on the principal branch, and a concrete negative argument never reaches this
    case (`log(−3)` is already `iπ + log(3)`, an Add);
  - `exp(c·log(w)) → w^c` folds unconditionally when `c` is numeric (then `w^c`
    *is* `exp(c·log w)` by definition); a symbolic `c` still requires `w > 0`;
  - `exp(Σ tᵢ)` pulls out every numeric-coefficient log term:
    `exp(log w₁ + c·log w₂ + r) = w₁ · w₂^c · exp(r)`.
- **Verified:** `exp(log x) → x`, `exp(2 log x) → x²`, `exp(log x / 2) → √x`,
  `exp(−log x) → 1/x`, `exp(log x + 1) → E·x`, `exp(log x + log y) → x·y`,
  `exp(log x − log y) → x/y`, `exp(log x + y) → x·exp(y)` — all match SymPy;
  `exp(y·log x)` (symbolic coefficient) and `exp(x + 1)` (no log term) are left
  unchanged, also matching SymPy.
- **Regression tests:** updated the `exp`/`log` inverse-pair cases and added
  `EXP-LOGSUM-1` in `tests/functions/exponential_test.cpp`. The two earlier
  tests that asserted the over-conservative "stays unevaluated" behavior were
  updated to the principal-branch result.
- **Note:** this intentionally replaces a previously deliberate branch-cut
  conservatism; the folding is mathematically exact (principal branch), so it
  introduces no incorrect results.

### EXPAND-POWBASE-1 — `expand` left a power of a product unflattened
- **Problem:** `expand((2x)²)` returned `(2x)²` rather than `4x²`; likewise
  `(xy)² → (xy)²`, `(3xy)² → (3xy)²`. The core `expand` only multinomial-expanded
  a power of an **Add** base — a power of a **Mul** base was never distributed,
  diverging from SymPy (and the convolution workaround in `solve_trig_reduce`
  existed precisely because of this).
- **Fix:** extended `expand_pow` in `src/core/expand.cpp` to distribute
  `(a·b·…)^n → a^n·b^n·…` when the exponent is any integer, or — for a
  non-integer exponent — when every base factor is provably positive (matching
  SymPy's `expand(power_base=True, force=False)`). The distributed product is
  re-expanded so numeric factors fold (`2² → 4`) and any `(Add)^n` factor
  multinomial-expands.
- **Verified:** `(2x)² → 4x²`, `(2x)³ → 8x³`, `(xy)² → x²y²`, `(3xy)² → 9x²y²`,
  `(2xy²)³ → 8x³y⁶`, `(2x)⁻² → 1/(4x²)`, `(−x)² → x²`,
  `((x+1)y)² → x²y²+2xy²+y²` — all match SymPy; a non-integer power over a
  possibly-negative factor (`(2x)^(1/2)`) is correctly left intact.
- **Regression test:** `EXPAND-POWBASE-1` in `tests/core/expand_test.cpp`
  (`[1i][expand][oracle][regression]`).
- **Scope:** integer exponents unconditionally; non-integer only under provable
  positivity. This is the conservative, sign-safe subset SymPy applies by
  default.

### SOLVE-TRIG-MULTIANGLE-1 — `solve` returned `[]` for mixed multiple-angle trig equations
- **Problem:** equations combining different integer multiples of one angle —
  `sin(x) − cos(2x)`, `cos(2x) + cos(x)`, `sin(2x) − sin(x)`,
  `cos(2x) + 3·sin(x) − 2` — all returned `[]`. The existing trig solvers only
  handle a single trig atom; a `cos(2x)` term alongside `sin(x)` is two distinct
  atoms with different arguments.
- **Fix:** added `solve_trig_reduce` in `src/solvers/solve.cpp`. It applies
  `expand_trig` to rewrite every multiple angle in terms of `sin θ`, `cos θ` for
  the common base angle θ, substitutes `s = sin θ`, `c = cos θ`, then reduces
  `s² → 1 − c²` to the canonical form `P(c) + s·Q(c) = 0`:
  - `Q ≡ 0` → `P(cos θ) = 0` (polynomial in cosine, both ± angles per root);
  - `P ≡ 0` → `sin θ = 0` together with `Q(cos θ) = 0`;
  - mixed → `s = −P/Q` with `s² = 1 − c²` gives `P² − (1−c²)Q² = 0`, a polynomial
    in `cos θ`; each in-range real root fixes `sin θ`'s sign (via `evalf`), hence
    a unique representative angle.
  The resultant is built by convolving coefficient vectors rather than squaring
  symbolically, since `expand` does not distribute a power over a product
  (`(2c)²` stays unflattened). Complex cosine roots (e.g. `sin x + cos x − 5`,
  which has no real solution) are rejected by a strict real-in-[−1,1] test.
- **Verified:** `sin(x)−cos(2x)`, `cos(2x)±cos(x)`, `cos(2x)±sin(x)`,
  `sin(2x)−sin(x)`, `sin(2x)−cos(x)`, `cos(2x)+3sin(x)−2`, `cos(3x)−cos(x)` and
  more — every returned root substitutes back to satisfy the equation and the
  real solution set is complete (checked against SymPy); `sin(x)+cos(x)−5`
  stays empty.
- **Regression test:** `SOLVE-TRIG-MULTIANGLE-1` in
  `tests/solvers/solve_test.cpp` (`[10][solve][trig][oracle][regression]`).
- **Scope:** a single base angle affine in the variable. Genuinely
  transcendental mixed equations (`tan x = x`, `sin x + x`) remain separate gaps.

### ASIN-ACOS-NEGSPECIAL-1 — `asin`/`acos` left negative √-special values unevaluated
- **Problem:** `acos(−√3/2)`, `acos(−√2/2)`, `asin(−√3/2)` stayed unevaluated,
  while `acos(−1/2)` folded correctly. The odd-identity reduction
  `asin(−x) = −asin(x)` is driven by `strip_neg`, which only recognized a
  leading integer `−1`; a Mul with a negative *rational* coefficient
  (`−½·√3`) was not stripped, so the positive-argument special-value table never
  ran.
- **Fix:** generalized `strip_neg` in `src/functions/trigonometric.cpp` to pull
  any negative numeric leading coefficient (`−½·√3 → ½·√3`, `−3·g → 3·g`), not
  just `−1·g`.
- **Verified:** `asin(−√3/2) = −π/3`, `asin(−√2/2) = −π/4`, `acos(−√3/2) = 5π/6`,
  `acos(−√2/2) = 3π/4`. This also cleans up `solve_trig_reduce`'s output (clean
  `5π/6` instead of `acos(−√3/2)`).
- **Regression test:** extended `asin`/`acos` exact-special-argument cases in
  `tests/functions/inverse_trig_test.cpp`.

### SOLVE-TRIG-PHASE-1 — `solve` returned `[]` for trig arguments with an additive phase
- **Problem:** `solve(sin(x+1)-1/2)`, `solve(cos(2x+π/3))`, `solve(tan(x+1)-1)`
  and similar all returned `[]`. The trig solvers (`solve_trig`,
  `solve_trig_poly`) only accepted a bare `B·x` argument — any additive phase
  inside the inner function (`x+1`, `2x+π/3`) failed the `arg == B·var` test.
- **Fix:** added an `affine_arg` helper in `src/solvers/solve.cpp` that
  decomposes a trig argument as `B·var + P` (B, P var-free, B ≠ 0) via a
  degree-1 polynomial read. `append_trig_roots` now takes the phase `P` and
  inverts each principal angle θ as `var = (θ − P)/B`. Both `solve_trig` and
  `solve_trig_poly` route their argument through `affine_arg`, so a bare `B·x`
  (P = 0) is just the special case.
- **Verified:** `sin(x+1)-1/2`, `cos(2x+π/3)`, `sin(2x+1)`, `tan(x+1)-1`,
  `sin(x+1)²-1/4` each return representative roots that substitute back to
  satisfy the equation (checked against SymPy); the bare-argument and
  scaled-argument cases (`2sin(x)-1`, `cos(2x)-1/2`) are unchanged.
- **Regression test:** `SOLVE-TRIG-PHASE-1` in
  `tests/solvers/solve_test.cpp` (`[10][solve][trig][oracle][regression]`,
  6 assertions).
- **Scope:** single trig atom with an argument affine in the variable. A
  nonlinear inner argument (`sin(x²)`) or genuinely transcendental mixed
  equations (`tan(x) = x`, `sin(x) + x`) remain separate gaps.

### INT-BIQUAD-NUM-1 — `∫x²/(x⁴+1)` (even numerator over a biquadratic) was unevaluated
- **Problem:** `INT-BIQUADRATIC-1` only handled `1/(biquadratic)`. A polynomial
  numerator over the same irreducible biquadratic — `∫x²/(x⁴+1)`,
  `∫(x²+1)/(x⁴+1)`, `∫(3x²+2)/(x⁴+4)` — still came back unevaluated, since the
  rational path factors only over ℚ.
- **Fix:** generalized `try_biquadratic` in `src/integrals/integrate.cpp` to
  accept `N·base^(-1)` where `N` is an even polynomial of degree ≤ 2
  (`n₀+n₂·x²`; odd numerators are left to the cleaner substitution paths). The
  numerator is distributed across the two real-quadratic partial fractions
  (`s_q=n₀/b`, `d_p=(n₂−s_q)/a`), each `(α·x+β)/(x²±a·x+b)` piece integrating to
  `log + arctan` as before. To keep clean closed forms, `try_rational` now
  defers (returns `nullopt`) whenever a monomial substitution applies — so
  `x²/(x⁶+1)` still yields `⅓·atan(x³)` rather than a partial-fraction expansion,
  while non-candidate forms like `x²/(x⁴+1)` reach the biquadratic path.
- **Verified:** `x²/(x⁴+1)`, `(x²+1)/(x⁴+1)`, `x²/(x⁴+x²+1)`, `1/(x⁴+1)`,
  `(3x²+2)/(x⁴+4)` each differentiate back to the integrand; `x²/(x⁶+1)`,
  `x/(x²+1)` and the prior biquadratic/rational cases are unchanged.
- **Regression test:** `INT-BIQUAD-NUM-1` in
  `tests/integrals/integrate_test.cpp` (`[7][integrate][oracle][regression]`,
  8 assertions).
- **Scope:** even numerators (degree ≤ 2) over an irreducible biquadratic. Odd
  numerators and higher-degree numerators remain on the substitution/rational
  paths.

### LIMIT-RAT-SYM-1 — limit of a rational function at ∞ broke with symbolic coefficients
- **Problem:** `limit((x+a)/x, x, ∞)` returned `0` (should be `1`), and
  `a·x/(x+1)` stayed unevaluated. Numeric-coefficient rationals (`(2x+1)/(x+1)→2`)
  worked, but a symbolic (var-free) coefficient broke direct ∞ substitution (it
  collapsed `∞·0` to `0`) and the L'Hôpital fallback. As a knock-on,
  `(1+a/x)^x` (the `1^∞` form) returned `nan` instead of `eᵃ`, because its base
  limit `1+a/x → 1` was being computed as `0`.
- **Fix:** added `rational_limit_at_infinity` in `src/calculus/limit.cpp`, run at
  the top of `limit_impl` for an infinite target (before the unreliable
  substitution). It splits `N/D`, requires polynomial var-free coefficients, and
  returns the limit by degree comparison and the leading-coefficient ratio
  (`deg N < deg D → 0`, `= → lead_N/lead_D`, `> → ±∞`).
- **Verified:** `(x+a)/x → 1`, `a·x/(x+1) → a`, `(a·x²+1)/(x²+x) → a`,
  `(x²+a)/(x+1) → ∞`, and `(1+a/x)^x → eᵃ` all match `sympy.limit`; numeric cases
  unchanged.
- **Regression test:** `LIMIT-RAT-SYM-1` in
  `tests/calculus/series_limit_test.cpp` (`[6][limit][infinity][oracle][regression]`,
  7 assertions).
- **Scope:** rational functions with var-free coefficients. A symbolic exponent
  rate `(1+a/x)^(b·x)` still needs `lim b·x = sign(b)·∞`, which depends on `b`'s
  unknown sign — left unevaluated.

### INT-BIQUADRATIC-1 — `∫1/(x⁴+1)` (irreducible biquadratic) was unevaluated
- **Problem:** `∫1/(x⁴+1)` — the canonical biquadratic that is irreducible over ℚ
  — and its scaled relatives (`1/(x⁴+4)`, `1/(2x⁴+2)`) came back unevaluated. The
  rational path only factors over ℚ; `x⁴+1` factors only over ℚ(√2). This was the
  long-standing headline integration gap (previously deferred to a full
  Lazard–Rioboo–Trager implementation).
- **Fix:** added `try_biquadratic` in `src/integrals/integrate.cpp` (after
  `try_arctan_quadratic`). For a denominator `a₄x⁴+a₂x²+a₀` it normalizes to
  `x⁴+p·x²+q` and, when `q>0` and `|p|<2√q` (real irreducible quadratic factors),
  factors `x⁴+p·x²+q = (x²+a·x+b)(x²−a·x+b)` with `b=√q`, `a=√(2√q−p)`. Partial
  fractions give two `(P·x+Q)/(irreducible quadratic)` pieces, each integrated
  directly to a `log + arctan` (handling the irrational `a`,`b` that the generic
  `try_linear_over_quadratic` can't). Runs only when `try_rational` (which handles
  the ℚ-factorable biquadratics like `x⁴+x²+1`) has already failed.
- **Verified:** `1/(x⁴+1)`, `1/(x⁴+4)`, `1/(2x⁴+2)`, `1/(x⁴+9)`, `1/(3x⁴+12)` each
  differentiate back to the integrand; the ℚ-factorable (`x⁴+x²+1`, `x⁴−1`) and
  quadratic cases are unchanged.
- **Regression test:** `INT-BIQUADRATIC-1` in
  `tests/integrals/integrate_test.cpp` (`[7][integrate][oracle][regression]`,
  9 assertions).
- **Scope:** biquadratics (`x⁴+p·x²+q`, no `x`/`x³` term) with real irreducible
  factors. A general irreducible quartic with an `x³` or `x` term, or one needing
  a higher algebraic extension, is still a separate (LRT-scale) gap.

### INT-QUAD-IRRATIONAL-1 — `∫1/(quadratic)` failed for positive-discriminant irrational roots
- **Problem:** `∫1/(x²−3)`, `∫1/(2x²−3)`, `∫1/(x²+x−1)` came back unevaluated. The
  arctan handler only covered `Δ<0` (no real roots), and `try_rational` only
  factors over ℚ — so a quadratic with a positive discriminant but irrational
  roots (no rational factorization) fell through both.
- **Fix:** in `src/integrals/integrate.cpp`, `try_arctan_quadratic` now handles
  the `disc < 0` (i.e. `Δ = b²−4ac > 0`) branch too:
  `∫1/(a·x²+b·x+c) = [log(2a·x+b−√Δ) − log(2a·x+b+√Δ)]/√Δ`. It only reaches this
  branch for irrational roots, since rational-root quadratics are split into
  clean logs by `try_rational`, which runs first.
- **Verified:** `1/(x²−3)`, `1/(3−x²)`, `1/(2x²−3)`, `1/(x²−2)`, `1/(x²+x−1)`,
  `1/(5x²−7)` each differentiate back to the integrand; the rational-root
  (`1/(x²−1)`) and `Δ<0` (`1/(x²+1)`) cases are unchanged. This also retroactively
  closes `∫1/(x·√(2x+3))` (INT-LINRADICAL-SUB-1's documented limitation), whose
  reduced `∫2/(u²−3)` now resolves.
- **Regression test:** `INT-QUAD-IRRATIONAL-1` in
  `tests/integrals/integrate_test.cpp` (`[7][integrate][oracle][regression]`,
  10 assertions).
- **Scope:** quadratics with rational coefficients; the irreducible-over-ℚ
  *quartic* `1/(x⁴+1)` (algebraic-extension factorization) is still a separate,
  larger gap.

### INT-LINRADICAL-SUB-1 — `integrate` missed the √(a·x+b) substitution
- **Problem:** integrands containing a radical of a non-trivial linear inner —
  `∫1/(x·√(x+1))`, `∫√(x+1)/x` — came back unevaluated. `try_radical_substitution`
  only handles `√x` (inner = var), not `√(a·x+b)`.
- **Fix:** added `try_linear_radical_substitution` in
  `src/integrals/integrate.cpp` (run after `try_radical_substitution`). It finds a
  `√(a·x+b)` factor (linear base, `b ≠ 0`), substitutes `x = (u²−b)/a`,
  `dx = (2u/a) du` (mapping `√(a·x+b) → u` and `x → (u²−b)/a` via `xreplace`),
  integrates in `u`, and back-substitutes `u = √(a·x+b)`.
- **Verified:** `∫1/(x·√(x+1))`, `∫√(x+1)/x`, `∫1/(√(x+1)+1)`, `∫x·√(x+1)` each
  differentiate back to the integrand; a bare `√(x+1)` still goes through the
  power rule.
- **Regression test:** `INT-LINRADICAL-SUB-1` in
  `tests/integrals/integrate_test.cpp` (`[7][integrate][oracle][regression]`,
  9 assertions).
- **Scope:** the reduced `u`-integral must be solvable downstream — e.g.
  `∫1/(x·√(2x+3))` reduces to `∫2/(u²−3)`, which needs a √3 factorization the
  rational integrator doesn't do, so it stays unevaluated (graceful, not wrong).

### INT-RADICAL-SUB-1 — `integrate` missed the radical substitution u = x^(1/d)
- **Problem:** integrands that are functions of a radical `x^(1/d)` came back
  unevaluated — `∫exp(√x)`, `∫sin(√x)`, `∫1/(1+√x)`, `∫1/(1+x^(1/3))` all returned
  `Integral(…)`, where SymPy gives elementary closed forms.
- **Fix:** added `try_radical_substitution` in `src/integrals/integrate.cpp`
  (run after `try_heurisch`). It takes `d = lcm` of the denominators of every
  `var` exponent; if `d > 1` it substitutes `x = uᵈ` (`dx = d·u^(d-1) du`,
  rewriting each `x^e → u^(e·d)` via `xreplace`), integrates the now-rational/
  elementary integrand in `u`, and back-substitutes `u = x^(1/d)`. This parallels
  the `solve_radical_poly` substitution.
- **Verified:** `∫exp(√x)`, `∫sin(√x)`, `∫1/(1+√x)`, `∫1/(√x+x)`, `∫1/(1+x^(1/3))`
  each differentiate back to the integrand (oracle), with the explicit
  `∫exp(√x) = 2√x·exp(√x) − 2exp(√x)`; the power rule still handles plain `√x`.
- **Regression test:** `INT-RADICAL-SUB-1` in
  `tests/integrals/integrate_test.cpp` (`[7][integrate][oracle][regression]`,
  12 assertions).
- **Scope:** a single radical generator `x^(1/d)`; nested or mixed radicals
  (`√(x+√x)`) need a different substitution.

### SUM-POLYGEOM-1 — `summation` of P(k)·rᵏ only worked for degree 1
- **Problem:** `Σ k²·2ᵏ`, `Σ (2k+1)·3ᵏ`, `Σ k²/2ᵏ` came back unevaluated; only the
  hardcoded degree-1 `Σ k·rᵏ` formula existed.
- **Fix:** added `sum_poly_geometric` in `src/calculus/summation.cpp`. For a
  summand `P(k)·r^(c·k+d)` (P a polynomial of degree ≥1, ratio `r = base^c` a
  concrete number ≠1) it finds the antidifference `S(k) = Q(k)·rᵏ`, where `Q` has
  the same degree and satisfies `r·Q(k+1) − Q(k) = P(k)`. The coefficients solve
  top-down: `q_j = [P_j − r·Σ_{i>j} C(i,j) q_i]/(r−1)`. The sum is then
  `S(hi+1) − S(lo)`.
- **Verified:** `k²·2ᵏ`, `k³·2ᵏ`, `k⁴·2ᵏ`, `(2k+1)·3ᵏ`, `(k²−k+1)·3ᵏ`, `k²/2ᵏ`
  all checked equal to `sympy.summation`, plus a concrete numeric range; the
  degree-1 path is unchanged.
- **Regression test:** `SUM-POLYGEOM-1` in
  `tests/calculus/series_limit_test.cpp` (`[6][summation][oracle][regression]`,
  8 assertions).
- **Scope:** concrete numeric ratio `r ≠ 1` (a symbolic ratio would need a
  `Piecewise` on `r = 1`, matching SymPy's restriction).

### SUM-FAULHABER-1 — `summation` of kᵖ only worked for p ∈ {2, 3}
- **Problem:** `Σ_{k=1}^n k⁴` (and any higher power) came back as an unevaluated
  `Sum(…)`; only the hardcoded `k²` and `k³` formulas existed.
- **Fix:** in `src/calculus/summation.cpp`, replaced the special cases with
  general Faulhaber's formula `Σ_{k=1}^n kᵖ = 1/(p+1)·Σ_{j=0}^p C(p+1,j) B_j
  n^(p+1−j)`, using the `bernoulli` numbers (added in BERNOULLI-EULER-1) and
  binomial coefficients. Evaluated for any positive integer `p` (capped at 200
  for cost) and over an arbitrary range via `F(hi) − F(lo−1)`.
- **Verified:** `Σ kᵖ` for `p = 2…15` checked equal to `sympy.summation`,
  including a partial range `Σ_{k=2}^n k⁴`; no leftover `Sum()` marker.
- **Regression test:** `SUM-FAULHABER-1` in
  `tests/calculus/series_limit_test.cpp` (`[6][summation][oracle][regression]`,
  13 assertions). The existing `k²`/`k³` oracle tests still pass on the new
  (unfactored but equivalent) form.
- **Scope:** integer powers of the bare summation variable; `Σ P(k)` for a
  general polynomial `P` already works through the existing Add-linearity path.

### SOLVE-EXPSUM-1 — `solve` returned `[]` for a Laurent polynomial in exp(x)
- **Problem:** `exp(x)+exp(-x)-2` and `exp(2x)-3·exp(x)+2` came back empty (or
  incomplete). They mix several `exp(m·x)` atoms, so `solve_exp_log_poly`
  (single atom, unit rate) bailed; SymPy returns `[0]` and `[0, log(2)]`.
- **Fix:** added `solve_exp_sum` in `src/solvers/solve.cpp`. It collects every
  `exp(m·x)` (integer m), substitutes `u = exp(x)` so `exp(m·x) → uᵐ`, solves the
  resulting equation in `u` recursively (the rational/poly paths clear the
  negative powers from `exp(-x)`), and inverts each root via `x = log(u)`. Because
  the multiplicity lives in the polynomial in `u`, scaled exponents now yield
  every complex representative, matching SymPy: `exp(2x)=1 → {0, iπ}`,
  `exp(3x)=1 → {0, ±2iπ/3}`. This also closes the previously-deferred composite
  case `exp(2x)-3·exp(x)+2`.
- **Verified:** `exp(x)±exp(-x)[-2]`, `exp(2x)-{3,5}exp(x)+{2,6}`, `exp(2x)-1`,
  `exp(3x)-1`, `exp(x)+exp(-x)-5/2` all checked equal to `sympy.solve`.
- **Regression test:** `SOLVE-EXPSUM-1` in `tests/solvers/solve_test.cpp`
  (`[10][solve][transcendental][oracle][regression]`, 6 assertions).
- **Scope:** exponents must be integer multiples of `x` (`exp(x/2)` would need a
  finer base). A bare `x` outside the exponentials (`x·eˣ`) leaves it to the
  LambertW path.

### SOLVE-RADPOLY-1 — `solve` returned `[]` for a polynomial in a radical x^(1/d)
- **Problem:** `x − √x − 2` came back empty, where SymPy gives `[4]` (a quadratic
  in `√x`: `u²−u−2=0 → u=2`, dropping `u=−1`). The polynomial path can't build a
  `Poly` through the fractional power, and SOLVE-RAD-1 only inverts the single
  form `gᵖ=c`.
- **Fix:** added `solve_radical_poly` in `src/solvers/solve.cpp`. It collects
  every `x^e` (e rational) and the bare `var`, takes `d = lcm` of the exponent
  denominators, substitutes `t = x^(1/d)` (rewriting each `x^e → t^(e·d)` via
  `xreplace`), solves the polynomial in `t`, and back-substitutes `x = t^d`. Each
  candidate is verified against the original equation, so extraneous roots
  (`√x = −1 ⇒ x = 1`) are dropped automatically.
- **Verified:** `x−√x−2 → [4]`, `x−3√x+2 → [1,4]`, `x+√x−6 → [4]` (extraneous 9
  dropped), `x−√x → [0,1]`, `x^(1/3)−2 → [8]` all equal to `sympy.solve`; plain
  polynomials untouched.
- **Regression test:** `SOLVE-RADPOLY-1` in `tests/solvers/solve_test.cpp`
  (`[10][solve][oracle][regression]`, 6 assertions).
- **Scope:** a single radical generator `x^(1/d)`. Mixed independent radicals
  (`√x + √(x+1) − 3`) need rationalization first and remain a follow-up.

### BERNOULLI-EULER-1 — `bernoulli` and `euler` numbers were missing
- **Problem:** `bernoulli(4)` and `euler(4)` parsed only as undefined functions,
  where SymPy gives `−1/30` and `5`.
- **Fix:** added `Bernoulli` and `Euler` functions (FunctionId, classes, builders,
  parser entries) in `src/functions/combinatorial.cpp`, each from its binomial
  recurrence. `bernoulli(n)` builds `Bₖ` exactly as rationals via
  `B_m = −1/(m+1)·Σ_{k<m} C(m+1,k) Bₖ`, returning SymPy's convention `B₁ = +1/2`;
  `euler(n)` builds the secant numbers via `E_{2m} = −Σ_{k<m} C(2m,2k) E_{2k}`.
  Odd `Bₙ>1` and odd `Eₙ` vanish; symbolic/negative arguments stay unevaluated.
- **Verified:** `bernoulli(0…20)` and `euler(0…16)` checked against SymPy,
  including `bernoulli(12)=−691/2730` and `euler(10)=−50521`.
- **Regression test:** `BERNOULLI-EULER-1` in
  `tests/functions/combinatorial_test.cpp` (`[3i][bernoulli][euler][oracle]`,
  16 assertions).
- **Scope:** integer-index Bernoulli/Euler numbers only; the polynomial forms
  `bernoulli(n, x)` / `euler(n, x)` remain. This completes the elementary
  special-number family (harmonic, bernoulli, euler).

### HARMONIC-FACT2-1 — `harmonic` and `factorial2` were missing
- **Problem:** `harmonic(5)` and `factorial2(5)` parsed only as undefined
  functions, where SymPy gives `137/60` and `15`.
- **Fix:** added `Harmonic` and `Factorial2` functions (FunctionId, classes,
  builders, parser entries) in `src/functions/combinatorial.cpp`. `harmonic(n)`
  accumulates `Σ_{k=1}^n 1/k` exactly as an `mpq_class` and returns a Rational;
  `factorial2(n)` multiplies `n(n−2)(n−4)…` down to 1 or 2, with the empty-product
  conventions `factorial2(0)=factorial2(−1)=1`. Symbolic and out-of-domain
  arguments stay unevaluated, matching SymPy.
- **Verified:** `harmonic` on `{0,1,2,5,10,20,50,100}` and `factorial2` on
  `{0,1,2,5,6,7,10,15,20,−1}` checked against SymPy.
- **Regression test:** `HARMONIC-FACT2-1` in
  `tests/functions/combinatorial_test.cpp` (`[3i][harmonic][factorial2][oracle]`,
  13 assertions).
- **Scope:** single-argument `harmonic` (Hₙ) only; SymPy's generalized
  `harmonic(n, m) = Σ 1/kᵐ`, plus the `bernoulli` and `euler` numbers (which need
  a recurrence), remain.

### SOLVE-EXPBASE-1 — `solve` returned `[]` for constant-base exponentials a^x = c
- **Problem:** `2^x−8`, `3^x−9`, `5^x−3` all came back empty. `a^x` is a `Pow`
  with a numeric base (not the `exp` function), so it never reached the
  transcendental solve path. SymPy returns `[3]`, `[2]`, `[log(3)/log(5)]`.
- **Fix:** added `solve_const_base_exp` in `src/solvers/solve.cpp` (tried right
  after `solve_rational`, since these skip the function-of-var branch). For
  `A·a^x + C = 0` it returns `x = log(−C/A)/log(a)`, collapsing to the exact
  integer exponent when `−C/A` is a power of `a` (`2^x=8 → 3`). `a^x=0` correctly
  yields no solution, and a negative RHS gives the complex log form matching
  SymPy (`2^x=−7 → (log7+iπ)/log2`).
- **Verified:** `2^x−{8,7,1}`, `3^x−9`, `5^x−3`, `10^x−100`, `6^x−36`, `2^x`,
  `2^x+7` all checked equal to `sympy.solve`.
- **Regression test:** `SOLVE-EXPBASE-1` in `tests/solvers/solve_test.cpp`
  (`[10][solve][transcendental][oracle][regression]`, 11 assertions).
- **Scope:** restricted to a non-perfect-power integer base and a bare exponent
  (`b=1`). For a perfect-power base (`4^x`) or scaled exponent (`2^(2x)`), SymPy
  enumerates extra complex representatives (the smaller complex period) that a
  single-root inversion would miss, so those stay unsolved — the same
  real-vs-complex completeness boundary as the `exp(B·x)` case (SOLVE-EXPLOG-POLY-1).

### SERIES-SINGULAR-1 — `series` emitted zoo/nan garbage at a singular point
- **Problem:** `series` built Taylor coefficients by substituting `x=x0` into each
  derivative. At a point where the function is singular this gave non-finite
  values that poisoned the whole sum: `series(log(x), x, 0)` →
  `zoo + x·zoo + …`, and removable singularities like `sin(x)/x` → `nan`. SymPy
  returns `log(x)`/`1/x`/`sqrt(x)` as-is and gives the proper Taylor series for
  `sin(x)/x`.
- **Fix in `src/calculus/series.cpp`:** when a coefficient `subs(dᵏf, x, x0)` is
  non-finite, fall back to `limit(dᵏf, x, x0)`. A removable singularity has a
  finite limit, so the correct Taylor coefficient is recovered
  (`sin(x)/x → 1 − x²/6 + …`); a genuine singularity stays non-finite, so the
  series can't be formed and the input is returned unchanged
  (`series(log(x),x,0)=log(x)`), matching SymPy.
- **Verified:** `log(x)`, `1/x`, `√x`, `1/x²` return as-is; `sin(x)/x`,
  `(eˣ−1)/x` give the SymPy Taylor series; ordinary analytic functions unchanged.
- **Regression test:** `SERIES-SINGULAR-1` in
  `tests/calculus/series_limit_test.cpp` (`[6][series][oracle][regression]`,
  7 assertions).
- **Scope:** genuine poles return the function unexpanded rather than a Laurent
  series (`cot(x)`, `eˣ/x` — SymPy gives `1/x − x/3 + …`), and a few removable
  cases whose high-order derivative limits don't resolve (`tan(x)/x`) also return
  unexpanded — graceful, never garbage. A true Laurent/Puiseux `series` is a
  larger effort.

### SOLVE-LAMBERT-1 — `solve` returned `[]` for Lambert-W equations
- **Problem:** `x·eˣ−1`, `x·e^(2x)−1`, `x·eˣ+1` all came back empty, where SymPy
  inverts them with the Lambert W function: `[W(1)]`, `[W(2)/2]`, `[W(−1)]`.
- **Fix:** added `solve_lambert` in `src/solvers/solve.cpp`. It detects the
  canonical form `a·x·exp(b·x)+c=0` (a, b var-free, b≠0, x to the first power,
  exp argument linear) and returns `x = W(−c·b/a)/b` using the defining identity
  `W(z)·e^(W(z))=z`. SymPP already had the `lambertw` function; this wires `solve`
  to produce it.
- **Verified:** `x·eˣ−{1,2}`, `x·eˣ+1`, `x·e^(2x)−1`, `2x·eˣ−1`, `x·e^(3x)−5`,
  and the homogeneous `x·eˣ → 0` all checked equal to `sympy.solve`.
- **Regression test:** `SOLVE-LAMBERT-1` in `tests/solvers/solve_test.cpp`
  (`[10][solve][transcendental][oracle][regression]`, 7 assertions). SOLVE-VAR-1's
  comment was updated — it no longer claims SymPP lacks a LambertW solver.
- **Update (rearrangement forms now covered):** `solve_lambert` was extended to
  the product-log `a·x·log(x)+c → exp(W(−c/a))` and the additive
  `x+eˣ+c → −c−W(e^(−c))` / `x+log(x)+c → W(e^(−c))` forms (additive forms
  require unit coefficients and argument `var`). Verified against SymPy for
  `x·log(x)∓{1,2}`, `x+log(x)+{−1,0,1}`, `x+eˣ+{−1,0,1}`, incl. the
  auto-evaluating `x+log(x)−1 → 1` and `x+eˣ−1 → 0`. Remaining: scaled arguments
  / non-unit coefficients in the additive forms, and `log(b·x)` in the product
  form.

### LIMIT-SIGN-1 — `limit` of a discontinuous `sign`/`abs` returned the point value
- **Problem:** `limit(sign(x), x, 0)` returned `0` — the point value `sign(0)=0`
  — instead of recognising the discontinuity. `sign(x²)` gave `0` (should be 1),
  and `|x|/x` gave `0` (via L'Hôpital → `sign(x)` → `0`). The two-sided limit of
  these is the sign of the inner just around the target, not at it.
- **Fix:** added `resolve_sign_limit` in `src/calculus/limit.cpp`, run at the top
  of `limit_impl`. It collects every `sign(g)`/`abs(g)` subexpression whose
  argument tends to 0, then evaluates the limit in the right- and
  left-neighborhoods with each replaced by its sampled one-sided value
  (`sign(g)→±1`, `abs(g)=g·sign(g)→g·(±1)`). Agreeing sides give the limit;
  disagreeing sides return `nan` (the two-sided limit does not exist). The
  side-sampling `side_sign_at` helper is shared with `signed_pole`.
- **Verified:** matches SymPy's two-sided `limit(…, '+-')` for `sign(x)`,
  `sign(-x)`, `sign(x³)` → DNE; `sign(x²)` → 1; `|x|/x`, `x/|x|`,
  `sin(x)/|x|` → DNE; `sign(x)·x` → 0; `|x|`, `|x−3|` continuous. Note `0·zoo`
  was already `nan` in core arithmetic — the bug was solely the point-value
  substitution.
- **Regression test:** `LIMIT-SIGN-1` in `tests/calculus/series_limit_test.cpp`
  (`[6][limit][regression]`, 9 assertions).
- **Scope:** when a side substitution leaves its own singularity (`|x|/x²`), the
  still-two-sided sub-limit yields `zoo` rather than the directional `+oo`; and
  `floor`/`ceiling`/`Heaviside` at their jumps still return the point value (a
  separate discontinuity class). SymPP keeps its two-sided convention throughout
  (SymPy defaults to the right-hand limit).

### ARITH-FN-1 — `mobius`, `divisor_count`, `divisor_sigma` were missing
- **Problem:** the multiplicative arithmetic functions stayed symbolic —
  `mobius(30)`, `divisor_count(12)`, `divisor_sigma(12)` parsed only as undefined
  functions, where SymPy gives `−1`, `6`, `28`.
- **Fix:** added `Mobius`, `DivisorCount`, `DivisorSigma` functions sharing one
  trial-division `factorize` helper in `src/functions/combinatorial.cpp`. From the
  `(prime, exponent)` list: `μ(n)` is 0 on any squared factor else `(−1)^#primes`;
  `σ₀(n)=∏(eᵢ+1)`; `σ₁(n)=∏(p^(eᵢ+1)−1)/(p−1)`. Symbolic and non-positive
  arguments stay unevaluated, matching SymPy.
- **Verified:** all three checked against SymPy for `{1,2,7,12,30,36,60,100,210,
  720,9973}`, including perfect numbers (`σ₁(6)=12`, `σ₁(28)=56`) and a large
  composite (`σ₁(720)=2418`).
- **Regression test:** `ARITH-FN-1` in `tests/functions/combinatorial_test.cpp`
  (`[3i][mobius][divisor][oracle]`, 17 assertions).
- **Scope:** single-argument `divisor_sigma` (σ₁) only; SymPy's two-argument
  `divisor_sigma(n, k)` (σ_k) and `isprime`/`factorint`/`divisors` remain.

### PRIME-PRIMEPI-1 — `prime(n)` and `primepi(n)` were missing
- **Problem:** following TOTIENT-1, the prime-indexing/counting functions stayed
  symbolic — `prime(5)`, `primepi(10)` parsed only as undefined functions, where
  SymPy gives `11` and `4`.
- **Fix:** added `Prime` and `PrimePi` functions (FunctionId, classes, builders,
  parser entries) in `src/functions/combinatorial.cpp` /
  `include/sympp/functions/combinatorial.hpp`. `prime(n)` walks `mpz_nextprime`
  `n` times for the n-th prime; `primepi(n)` counts primes ≤ n the same way
  (clamping to 0 below 2). Symbolic args and a non-positive `prime` index stay
  unevaluated, matching SymPy.
- **Verified:** `prime` on `{1,5,10,100,1000}` and `primepi` on
  `{1,2,10,100,10000}` checked against SymPy, plus the round-trips
  `prime(primepi(13))=13` and `primepi(prime(7))=7`.
- **Regression test:** `PRIME-PRIMEPI-1` in
  `tests/functions/combinatorial_test.cpp` (`[3i][prime][primepi][oracle]`,
  16 assertions).
- **Scope:** the `mpz_nextprime` walk is linear in the index/count (safety-bounded
  at `prime(10⁶)` / `primepi(10⁸)`); a very large `primepi` would want a
  Meissel–Mertens sieve. `isprime`, `factorint`, `divisors` remain.

### TOTIENT-1 — Euler's totient `totient(n)` was missing
- **Problem:** `totient(n)` parsed only as an undefined function and never
  evaluated, where SymPy's `totient` computes Euler's φ for positive integers
  (`totient(12)=4`, `totient(7)=6`).
- **Fix:** added a `Totient` function (FunctionId, class, builder, parser entry)
  in `src/functions/combinatorial.cpp` / `include/sympp/functions/combinatorial.hpp`.
  For a positive Integer it computes `φ(n) = n·∏_{p|n}(1−1/p)` by trial-dividing
  out each distinct prime; symbolic and non-positive arguments stay unevaluated,
  matching SymPy.
- **Verified:** values for `{1,2,7,12,17,36,100,1000000}` and a large composite
  (`totient(360360)=69120`) checked against SymPy.
- **Regression test:** `TOTIENT-1` in `tests/functions/combinatorial_test.cpp`
  (`[3i][totient][oracle]`, 11 assertions).
- **Scope:** trial-division factorization is fine for everyday inputs; a
  cryptographically large `n` with two huge prime factors would be slow (SymPy
  has the same characteristic). Related number-theory functions (`isprime`,
  `primepi`, `factorint`) remain unimplemented.

### SOLVE-INVFN-1 — `solve` returned `[]` for inverse trig/hyperbolic equations
- **Problem:** `asin(x)−1`, `atan(x)−1`, `asinh(x)−2`, … all came back empty,
  where SymPy returns `[sin(1)]`, `[tan(1)]`, `[sinh(2)]` — the forward-function
  inverse. Range-violating equations (`asin(x)−2`, with `2 > π/2`) should give
  `[]`.
- **Fix:** added `solve_inverse_func_poly` in `src/solvers/solve.cpp`. It detects
  a polynomial in one inverse atom `g⁻¹(B·x)` (`g⁻¹ ∈
  {asin,acos,atan,asinh,acosh,atanh}`), substitutes/solves for the inner value
  `c`, and inverts with the forward function: `g⁻¹(B·x)=c → x = g(c)/B`. Each `c`
  is range-checked against the inverse function's codomain (asin `[−π/2,π/2]`,
  acos `[0,π]`, atan `(−π/2,π/2)`, acosh `[0,∞)`; asinh/atanh unbounded) via a
  numeric `evalf`, so out-of-range roots are dropped. Inverse functions are
  single-valued, so any `B` is handled.
- **Verified:** fourteen equations checked equal to `sympy.solve`, including an
  auto-evaluating RHS (`asin(x)=π/6 → 1/2`), a scaled argument
  (`atan(2x)=1 → tan(1)/2`), a quadratic (`asin(x)²=1 → ±sin(1)`), and three
  range rejections.
- **Regression test:** `SOLVE-INVFN-1` in `tests/solvers/solve_test.cpp`
  (`[10][solve][transcendental][oracle][regression]`, 12 assertions).

### SOLVE-RATIONAL-1 — `solve` returned `[]` for rational equations
- **Problem:** any equation with a var-dependent denominator came back empty —
  `x+1/x−2`, `1/x−1/2`, `1/(x−1)+1/(x+1)`, `(x²−1)/(x−1)` all yielded `[]`, where
  SymPy returns `[1]`, `[2]`, `[0]`, `[−1]`. The polynomial path can't build a
  `Poly` from a `1/x` term, and rational equations carry no Function/radical so
  they skipped the transcendental branch too.
- **Fix:** added `solve_rational` in `src/solvers/solve.cpp` (tried after the
  polynomial path fails). It `together()`s the equation into `N/D`, solves
  `N(var)=0` recursively, and discards any root that vanishes the denominator
  (`subs(D, var, r) = 0`) — so a removable pole like `x=1` in `(x²−1)/(x−1)` is
  dropped rather than returned.
- **Verified:** solution sets checked set-equal to `sympy.solve` for ten
  equations, including pole removal, an irrational two-root case
  (`1/x+1/(x−1)−2`), and no-solution constant numerators (`1/(x+1)−1/(x−1)`).
- **Regression test:** `SOLVE-RATIONAL-1` in `tests/solvers/solve_test.cpp`
  (`[10][solve][oracle][regression]`, 6 assertions).
- **Scope:** denominators of integer powers (`1/(x−1)²`); a radical denominator
  (`1/√x`) stays the radical path's job.

### SOLVE-EXPLOG-POLY-1 — `solve` returned `[]` for a polynomial in exp(x) or log(x)
- **Problem:** `solve` handled a single `exp(x)−c` / `log(x)−c` but came back empty
  for any higher-degree polynomial in one transcendental atom — `exp(x)²−3exp(x)+2`,
  `log(x)²−1`, `log(x)²−3log(x)+2` all yielded `[]`, where SymPy returns
  `[0, log(2)]`, `[e, 1/e]`, `[e, e²]`. (This is the exp/log analogue of the trig
  case closed in SOLVE-TRIG-POLY-1.)
- **Fix:** added `solve_exp_log_poly` in `src/solvers/solve.cpp`. It detects a
  polynomial in exactly one `exp`/`log` atom `g(B·x)`, substitutes `u = g(B·x)`,
  solves the polynomial in `u` with the `Poly` root finder, and inverts each
  root: `exp(B·x)=c → log(c)/B` (skipping `c=0`, which has no solution),
  `log(B·x)=c → exp(c)/B`. Complex roots are included where SymPy includes them
  (`exp(x)=−1 → iπ`).
- **Verified:** solution sets checked set-equal to `sympy.solve` for nine
  equations (quadratic/cubic in exp and log, repeated roots, a scaled log
  argument, complex roots).
- **Regression test:** `SOLVE-EXPLOG-POLY-1` in `tests/solvers/solve_test.cpp`
  (`[10][solve][transcendental][oracle][regression]`, 5 assertions).
- **Scope:** the `exp` inversion is taken only for `B=1` — `exp(B·x)=c` with
  `B≠1` has `B` complex representatives (period `2πi/B`) whose enumeration is out
  of scope; `log` is single-valued, so any `B` is handled. A composite like
  `exp(2x)` alongside `exp(x)` is two distinct atoms and needs an
  `exp(2x)=exp(x)²` rewrite — a follow-up.

### RADSIMP-SIMPLIFY-2 — `simplify` didn't rationalize two-surd denominators
- **Problem:** following RADSIMP-SIMPLIFY-1, denominators that are a sum of two
  surds with no rational part (`√3−√2`, `√5+√2`) were still left as reciprocals,
  where SymPy returns `√2+√3`, `(√5−√2)/3`, etc. Two causes: `radsimp` only
  handled a single `a+b√c` binomial; and even when extended, a result with a
  non-unit rational denominator (`(√5−√2)/3`) has a larger node count than the
  reciprocal, so `simplify`'s anti-bloat guard reverted it.
- **Fix in `src/simplify/simplify.cpp`:** (1) `radsimp` now also rationalizes
  `b₁√c₁ ± b₂√c₂` via the conjugate `b₁√c₁ ∓ b₂√c₂` (product `b₁²c₁ − b₂²c₂`,
  rational); (2) the anti-bloat guard in `simplify` is relaxed via a new
  `has_surd_denominator` check — when the pipeline removes a surd denominator
  that the input still carries, the (possibly larger) rationalized form is kept,
  while ordinary expansion bloat is still rejected.
- **Verified:** `1/(√3−√2)`, `1/(√5+√2)`, `1/(√7−√3)`, `2/(√3+√2)`,
  `x/(√5−√3)` all checked equal to SymPy and free of a surd reciprocal; the
  single-binomial cases (RADSIMP-SIMPLIFY-1) and unrelated forms are unchanged.
- **Regression test:** `RADSIMP-SIMPLIFY-2` in
  `tests/simplify/simplify_test.cpp` (`[5][simplify][radsimp][oracle][regression]`,
  6 assertions).
- **Scope:** a rational part plus two surds (`1/(1+√2+√3)`) needs a two-step
  conjugate and remains a follow-up.

### RADSIMP-SIMPLIFY-1 — `simplify` left a surd in a binomial denominator
- **Problem:** `simplify(1/(1+√2))` returned the reciprocal unchanged instead of
  `√2−1`. Two issues compounded: (1) `radsimp` only looked for the denominator
  inside a `Mul`, so a *bare* reciprocal `Pow(1+√2, −1)` was returned untouched;
  (2) even when it did rationalize, it produced a non-distributed `−(−√2+1)`
  whose node count exceeded the reciprocal, so `simplify`'s anti-bloat guard
  reverted to `1/(1+√2)`.
- **Fix in `src/simplify/simplify.cpp`:** `radsimp` now handles a bare
  `Pow(den, −1)` (num = 1) in addition to the `Mul` case, and `expand`s its
  result so the rationalized form is the compact `√2−1` — small enough to pass
  the guard.
- **Verified:** `simplify` of `1/(1+√2)`, `1/(2+√3)`, `1/(√5−1)`, `1/(3−2√2)`,
  and `x/(1+√2)` all checked equal to SymPy and free of a surd-binomial
  reciprocal.
- **Regression test:** `RADSIMP-SIMPLIFY-1` in
  `tests/simplify/simplify_test.cpp` (`[5][simplify][radsimp][oracle][regression]`,
  5 assertions).
- **Scope:** single binomial surd `a + b√c`. A two-surd denominator
  (`1/(√7−√3)`) still needs a multi-term conjugate and remains a follow-up.

### COMB-RATIO-1 — `combsimp`/`gammasimp` only cancelled ratios when the numerator was larger
- **Problem:** `simplify_func_ratio` cancelled `factorial(a)/factorial(b)` (and
  the gamma analogue) only when `a − b` was a *non-negative* integer. When the
  denominator held the larger argument it gave up: `factorial(k−1)/factorial(k)`
  stayed instead of `1/k`, `gamma(z)/gamma(z+1)` instead of `1/z`. This also
  blocked binomial-ratio simplification entirely — `binomial(n,k)/binomial(n,k−1)`
  stayed where SymPy returns `(n−k+1)/k`.
- **Fix (two parts) in `src/simplify/simplify.cpp`:**
  1. Added the negative-difference branch to `simplify_func_ratio`: for
     `a − b = −m` it emits `1/falling_factorial(b, m)` (factorial) or
     `1/falling_factorial(b−1, m)` (gamma).
  2. Added `combsimp_binomial_ratio`, which rewrites `binomial(a,b) =
     a!/(b!·(a−b)!)` inside a `Mul` so the now-bidirectional factorial
     cancellation closes binomial shifts — but only adopts the result when it
     fully resolves (no factorial/binomial left), so a lone `binomial(2n,n)`
     keeps its compact form rather than expanding into factorials.
- **Verified:** `factorial(k−1)/factorial(k) → 1/k`, `gamma(z)/gamma(z+1) → 1/z`,
  and `binomial(n,k)/binomial(n,k−1)`, `binomial(n+1,k)/binomial(n,k)`,
  `binomial(n,k)/binomial(n−1,k)`, `binomial(n,k)/binomial(n,k−2)` all checked
  equal to SymPy; non-reducible binomials/factorials are unchanged.
- **Regression test:** `COMB-RATIO-1` in `tests/simplify/simplify_test.cpp`
  (`[5][combsimp][oracle][regression]`, 7 assertions).

### TRIG-PI5-1 — `sin/cos/tan` of the pentagon angles (π/5, π/10) stayed unevaluated
- **Problem:** following TRIG-PI8-1, the remaining common special angles — the
  pentagon family π/5 (36°) and π/10 (18°) — were still symbolic, where SymPy
  gives `cos(π/5) = (1+√5)/4`, `cos(2π/5) = (√5−1)/4`, `tan(π/5) = √(5−2√5)`, and
  the `√(10±2√5)` nested radicals for the π/10 cosines.
- **Fix:** added den-5 (num 1,2) and den-10 (num 1,3) reference angles to
  `base_cos_pi`, and the four matching `tan` values to `base_tan_pi`, in
  `src/functions/trigonometric.cpp`. `sin` derives from the co-function
  reflection (`sin(π/10)=cos(2π/5)`, `sin(3π/10)=cos(π/5)`, …) and every multiple
  reduces through the existing period/reflection folds.
- **Verified:** all `sin/cos/tan` of `{1,2,3,4,6}·π/5` and `{1,3,7,9}·π/10` plus
  negatives checked equal to SymPy via the oracle.
- **Regression test:** `TRIG-PI5-1` in `tests/functions/trigonometric_test.cpp`
  (`[3b][trig][oracle][regression]`, 9 assertions).
- **Scope:** the special-angle table now covers denominators
  {1,2,3,4,5,6,8,10,12} — the standard constructible angles. Denominators like 7,
  9, 11 (non-constructible / `cos` not expressible in real radicals) stay
  symbolic, matching SymPy.

### TRIG-PI8-1 — `sin/cos/tan(π/8)` stayed unevaluated
- **Problem:** the special-angle table covered denominators {1,2,3,4,6,12} but
  not 8, so the π/8 family (22.5°, the half-angles of π/4) came back symbolic —
  `cos(π/8)`, `tan(π/8)`, … — where SymPy returns `√(2+√2)/2`, `√2−1`, etc. A
  prior test deliberately kept π/8 unevaluated, but that was a presentational
  preference (the π/12 entries already emit radicals), not a correctness
  constraint, so it diverged from SymPy for no reason.
- **Fix:** added the den-8 reference angles to `base_cos_pi` (`cos(π/8)=√(2+√2)/2`,
  `cos(3π/8)=√(2−√2)/2`) and `base_tan_pi` (`tan(π/8)=√2−1`, `tan(3π/8)=√2+1`) in
  `src/functions/trigonometric.cpp`. `sin` derives from the existing co-function
  reflection (`sin_pi(r)=cos_pi(½−r)`), and all multiples (5π/8, 7π/8, …) reduce
  through the period/reflection folds.
- **Verified:** every `sin/cos/tan` of `{±1,3,5,7,9,11}·π/8` checked equal to
  SymPy via the oracle.
- **Regression test:** `TRIG-PI8-1` in `tests/functions/trigonometric_test.cpp`
  (`[3b][trig][oracle][regression]`, 6 assertions) — replaces the old
  "stays unevaluated" assertion.
- **Scope:** the π/5 / π/10 pentagon family (`cos(π/5)=(1+√5)/4`, …) is still
  symbolic and remains a follow-up.

### INT-RATIONAL-NOPARTIAL-1 — `try_rational` leaked partial results with an `Integral` marker
- **Problem:** when `apart()` couldn't fully split a denominator, `try_rational`
  integrated the partial-fraction sum term by term and returned the half-answer —
  e.g. `∫x²/(x⁶+1)` gave `−⅓atan(x) + Integral((⅓x²+⅓)/(x⁴−x²+1))`. The leaked
  `Integral(…)` sat inside an `Add`, so the existing top-level `is_integral_marker`
  guard missed it, and the bogus partial shadowed cleaner strategies.
- **Fix:** added a recursive `contains_integral_marker` in
  `src/integrals/integrate.cpp` and used it for both intermediate guards and a
  final check in `try_rational`: if the assembled antiderivative still hides an
  `Integral` anywhere, `try_rational` returns `nullopt`. The integral then either
  falls through to a strategy that closes it — `∫x²/(x⁶+1) = ⅓atan(x³)` via the
  monomial substitution (INT-MONOMIAL-SUB-1) — or is returned honestly
  unevaluated (`∫1/(x⁶+1)`), never as a leaked partial.
- **Verified:** `∫x²/(x⁶+1)` closes to `⅓atan(x³)` (oracle); `∫1/(x⁶+1)` and
  `∫1/(x⁵+1)` carry no partial sum; fully-solvable rationals (`1/(x⁴−1)`,
  `1/(x³+1)`, `x/(x⁶+1)`) are unchanged.
- **Regression test:** `INT-RATIONAL-NOPARTIAL-1` in
  `tests/integrals/integrate_test.cpp` (`[7][integrate][oracle][regression]`,
  5 assertions).
- **Scope:** denominators that are irreducible over ℚ but elementary-integrable
  (`1/(x⁶+1)`, `1/(x⁵+1)`) still return unevaluated — closing those needs
  integration over an algebraic extension (Lazard–Rioboo–Trager), a larger
  effort. The fix removes the misleading partials in the meantime.

### INT-MONOMIAL-SUB-1 — `integrate` missed the monomial substitution u = x^d
- **Problem:** integrands of the form `x^(d-1)·f(x^d)` whose `f(x^d)` is hidden
  inside a rational or radical expression came back unevaluated — `∫x/(x⁴+1)`,
  `∫x³/(x⁸+1)`, `∫x/√(1−x⁴)` all returned `Integral(…)`, while SymPy gives
  `½atan(x²)`, `¼atan(x⁴)`, `½asin(x²)`. `try_heurisch` couldn't help: its
  substitution is structural and `x⁴` does not contain `x²` as a subexpression,
  so `u = x²` never linearised the denominator.
- **Fix:** added `try_monomial_substitution` in `src/integrals/integrate.cpp`
  (run before `try_heurisch`). For `d = 2…6` it forms `t = expr/(d·x^(d-1))`,
  rewrites every `x^k` with `d | k` to `u^(k/d)` (via `xreplace`), and — if no
  bare `x` survives — integrates the resulting `f(u)` and back-substitutes
  `u = x^d`. A `x^k` with `d ∤ k`, or a leftover `var`, aborts that `d`.
- **Verified:** each antiderivative differentiates back to its integrand
  (oracle), with the explicit `∫x/(x⁴+1) = ½atan(x²)`; unrelated integrands
  (`x/(x²+1)`, `1/(x²+1)`) are unchanged.
- **Regression test:** `INT-MONOMIAL-SUB-1` in
  `tests/integrals/integrate_test.cpp` (`[7][integrate][oracle][regression]`,
  8 assertions).
- **Scope:** numerator must be the exact `x^(d-1)` the substitution needs.
  Cases such as `∫x²/(x⁶+1)` (clean `⅓atan(x³)`) are still intercepted earlier by
  `try_rational`, which returns a partial result with a leftover `Integral` —
  a separate issue in the rational-integration path.

### TRIG-ANGLE-ADD-1 — `simplify` didn't fold the angle-addition identities
- **Problem:** `simplify` collapsed same-argument trig combinations (Pythagorean,
  power-reduction, double-angle) but left the two-argument angle-addition forms
  untouched: `sin(x)cos(y)+cos(x)sin(y)`, `cos(x)cos(y)−sin(x)sin(y)`, etc.
  stayed as products where SymPy returns `sin(x+y)`, `cos(x+y)`, …
- **Fix:** added `trigsimp_angle_addition` (with a `as_two_trig_term` classifier)
  in `src/simplify/simplify.cpp`, run inside `trigsimp_node`. On a two-term `Add`
  whose terms are products of two first-power sin/cos factors it recognises
  `sin(a)cos(b) ± cos(a)sin(b) → sin(a±b)` and
  `cos(a)cos(b) ∓ sin(a)sin(b) → cos(a±b)`, carrying a shared coefficient
  through. The classifier bails on any non-clean shape (squared trig, a third
  trig factor, a leftover function in the coefficient) so nothing is mis-folded.
- **Verified:** each result checked equal to `sympy.simplify` for all four
  identities, a coefficient-carrying case, and a composing case
  (`sin(2x)cos(x)+cos(2x)sin(x) → sin(3x)`).
- **Regression test:** `TRIG-ANGLE-ADD-1` in `tests/simplify/simplify_test.cpp`
  (`[5][trigsimp][oracle][regression]`, 7 assertions).
- **Scope:** matches a two-term `Add` with the exact product structure; it does
  not yet expand `sin(3x)` into single-argument powers (the reverse direction) or
  handle products of more than two angles.

### POLY-FACTOR-ROOTS-1 — `solve`/`Poly::roots` returned nested radicals for factorable quartics
- **Problem:** a quartic with no rational root but which factors over ℚ into two
  quadratics (e.g. `x⁴+x²+1 = (x²+x+1)(x²−x+1)`) went straight to Ferrari's
  resolvent, producing nested radicals like `sqrt((I·sqrt(3)−1)/2)` instead of
  the clean `±1/2 ± √3·i/2`. `x⁶−1` inherited the same mess for its degree-4
  cofactor, and a degree-≥5 polynomial with no rational root but a non-trivial
  factorization (e.g. `(x²+x+1)(x³+2)`) returned no roots at all.
- **Fix:** added `roots_by_factoring` in `src/polys/poly.cpp`, used by
  `Poly::roots()` for degree 4 (before the Ferrari fallback) and degree ≥5
  (before giving up). It calls the existing `kronecker_find_factor` to split the
  polynomial over ℚ and recurses on each factor, so each piece is solved by the
  cleanest applicable path (quadratic formula, Cardano, …). SymPy factors before
  solving for the same reason.
- **Verified:** root sets checked set-equal to `sympy.solve` for `x⁴+x²+1`,
  `x⁶−1`, `x⁸−1`, `(x²+x+1)(x³+2)` and the existing rational/biquadratic cases;
  the quartic roots carry no nested-radical wrapping.
- **Regression test:** `POLY-FACTOR-ROOTS-1` in `tests/polys/poly_test.cpp`
  (`[4][poly][roots][oracle][regression]`, 7 assertions).
- **Scope:** helps only polynomials that are *reducible* over ℚ. Genuinely
  irreducible quartics such as `x⁴+1` still go through Ferrari and keep a radical
  form (`(-I)^(1/2)` etc.) — correct but not simplified to `±√2/2 ± √2·i/2`;
  cleaning those is a separate radical-denesting gap.

### SUM-TELESCOPE-1 — `summation` returned unevaluated for telescoping rational sums
- **Problem:** `summation` handled polynomial (Faulhaber), geometric and
  arithmetic-geometric summands, but every rational summand `c/D(k)` came back as
  an unevaluated `Sum(…)`: `Σ 1/(k(k+1))`, `Σ 1/(k(k+2))`, `Σ 1/(4k²−1)` all had
  closed forms in SymPy (`n/(n+1)`, …) and none in SymPP.
- **Fix:** added `telescope_rational` in `src/calculus/summation.cpp`. For a
  summand `c/D(k)` where `D` is a quadratic with two distinct rational roots
  `r₁,r₂` whose difference `d = r₁−r₂` is a nonzero integer, partial fractions
  give `c/(lead·d)·[u(k) − u(k+d)]` with `u(k)=1/(k−r₁)`, which telescopes to
  `c/(lead·d)·[Σ_{j=0}^{d−1} u(lo+j) − Σ_{j=1}^{d} u(hi+j)]`. A pole inside the
  summation range (an integer root ≥ `lo`) is detected and the sum is left
  unevaluated rather than producing a bogus closed form (`Σ 1/(k(k−1))` stays).
- **Verified:** closed forms checked equal to `sympy.summation` for 9 summands
  (unit and non-unit leading coefficients, pole gaps `d∈{1,2,3}`, scaled
  numerators, shifted factors like `(3k−1)(3k+2)`).
- **Regression test:** `SUM-TELESCOPE-1` in
  `tests/calculus/series_limit_test.cpp` (`[6][summation][oracle][regression]`,
  6 assertions).
- **Scope:** limited to a denominator that is one quadratic with two distinct
  integer-spaced rational roots and a var-free numerator. Higher-degree
  denominators (≥3 linear factors) and non-integer-spaced roots — which need
  full partial-fraction grouping or Gosper's algorithm — remain open.

### SOLVE-TRIG-LINEAR-1 — `solve` returned `[]` for `a·sin(x)+b·cos(x)+c` (R-method)
- **Problem:** `solve` had no path for a linear combination of sin and cos of the
  same argument. `sin(x)+cos(x)`, `√3·sin(x)+cos(x)`, `sin(x)+cos(x)−1`,
  `3·sin(x)+4·cos(x)`, `sin(2x)+cos(2x)` all came back `[]` (two distinct trig
  atoms, so neither the single-term nor the polynomial-in-one-function path
  applied).
- **Fix:** added `solve_trig_linear` in `src/solvers/solve.cpp`. It recognises
  `a·sin(B·x)+b·cos(B·x)+c` (var-free `a,b,c`, shared argument `B·x`). The
  homogeneous case (`c=0`) reduces to `tan(B·x)=−b/a`, a single representative
  `atan(−b/a)/B`. Otherwise it applies the R-method: `a·sin+b·cos = R·sin(B·x+φ)`
  with `R=√(a²+b²)`, `φ=atan2(b,a)`, so `sin(B·x+φ)=−c/R` yields two
  representatives `(θ−φ)/B`. The `trig_value_in_range` guard makes
  `|c|>R` return no real solution (`sin(x)+cos(x)−5 → []`).
- **Verified:** every solution set checked set-equal to `sympy.solve` for 8
  equations (homogeneous, non-homogeneous, scaled argument, irrational
  coefficient, out-of-range); transcendental forms that don't share a closed
  shape were confirmed numerically equal.
- **Regression test:** `SOLVE-TRIG-LINEAR-1` in `tests/solvers/solve_test.cpp`
  (`[10][solve][trig][oracle][regression]`, 8 assertions).
- **Scope:** sin and cos must share the same argument. Multiple-angle mixes
  (`sin(2x)−sin(x)`, `cos(2x)−cos(x)`) and products (`sin(x)·cos(x)−1/2`, needing
  the double-angle identity) remain open.

### SOLVE-TRIG-POLY-1 — `solve` returned `[]` for a polynomial in one trig function
- **Problem:** `solve` handled a *single* trig term `A·f(B·x)+C=0` (SOLVE-TRIG-1)
  but came back empty for any higher-degree polynomial in one trig function:
  `sin(x)²−1`, `2·sin(x)²−1`, `sin(x)²−sin(x)`, `cos(x)²−1/4`, `tan(x)²−1`,
  `2·cos(x)²−cos(x)−1` all yielded `[]`. The single-term path bailed because the
  power form was only solved homogeneously (`f^n = 0`), so any constant term
  killed it, and the polynomial path can't see through `sin`.
- **Fix:** added `solve_trig_poly` in `src/solvers/solve.cpp`. It detects a
  polynomial in exactly one trig atom `f(B·x)` (`f ∈ {sin,cos,tan}`, `B`
  var-free), substitutes `u = f(B·x)`, solves the polynomial in `u` via the
  existing `Poly` root finder, then inverts each in-range root to representative
  angles over one principal period — matching SymPy's `solve` as a set. The
  per-function inversion is now a shared `append_trig_roots` helper reused by
  both `solve_trig` and `solve_trig_poly`. A numeric root with `|c|>1` for
  sin/cos contributes no real solution (`sin(x)²=4 → []`), via a `trig_value_in_range`
  guard; an in-range irrational `c` (e.g. `asin(1/3)`) is preserved unevaluated,
  exactly as SymPy reports it.
- **Verified:** every solution set checked set-equal to `sympy.solve` for 9
  equations (sin/cos/tan, scaled argument, quadratic-with-endpoints, out-of-range).
- **Regression test:** `SOLVE-TRIG-POLY-1` in `tests/solvers/solve_test.cpp`
  (`[10][solve][trig][oracle][regression]`, 8 assertions).
- **Scope:** still one trig function per equation. Mixed-function equations
  (`sin·cos`, `sin(2x)−sin(x)`) and the `a·sin+b·cos` R-method remain open; the
  structural zero-product path (SOLVE-TRIG-1) continues to cover products of
  distinct trig factors.

### ORACLE-XCHECK-1/2/3 — oracle tests only compared against hand-written literals
- **Problem (test-rig integrity, not a math gap):** all ~794 oracle assertions
  were `equivalent(sympp_output, "literal_I_typed")` — SymPy adjudicated the
  equivalence, but the *expected answer* was hand-authored. A wrong literal that
  happened to match a wrong SymPP result would pass undetected. The
  `oracle.diff` / `oracle.integrate` / `oracle.simplify` helpers (which let
  SymPy compute the reference answer independently) existed but were used zero
  times.
- **Fix:** added `tests/oracle/oracle_crosscheck_test.cpp`
  (`[0][oracle][crosscheck]`). For a spread of inputs, SymPP parses and computes
  `diff`/`integrate`/`simplify` while SymPy computes the *same* operation on the
  *same* input; the two independently-produced results must be `equivalent`.
  Neither side is hand-authored, so the engines genuinely check each other
  (12 derivative, 8 antiderivative, 6 simplification inputs).
- **Verified:** all 26 cross-check assertions pass on local `1.13.3`.
- **Regression test:** `tests/oracle/oracle_crosscheck_test.cpp`.
- **Scope:** indefinite integrals here share the `+C = 0` convention so the
  antiderivatives compare equal (not merely up to a constant). Inputs are
  limited to operations both engines support; expanding the cross-checked
  surface is incremental future work.

### VERSION-GUARD-1 — the oracle accepted any SymPy version silently
- **Problem (test-rig integrity, not a math gap):** the whole suite trusts the
  SymPy oracle to adjudicate `equivalent(...)`, but nothing pinned *which* SymPy
  it ran against. CI installed `sympy==1.14` while local development ran
  `1.13.3`, and the smoke test only checked the version string was non-empty. A
  silent upgrade to an unvetted release could change a canonical form or add an
  auto-simplification and quietly alter the "equivalent" verdict with no test
  noticing.
- **Fix:** added a `[0][oracle]` smoke test (`VERSION-GUARD-1`,
  `tests/oracle/oracle_smoke_test.cpp`) that parses the oracle's reported SymPy
  version and `FAIL`s loudly unless its `major.minor` is on a validated
  allowlist (`1.13` / `1.14`). Bumping SymPy now requires re-validating the
  suite and extending the allowlist in the same commit. The `.github/workflows/
  ci.yml` pin and the `docs/09-known-issues.md` header were annotated to match.
- **Verified:** passes on local `1.13.3` and CI `1.14`; a standalone check
  confirms the guard fires on `1.15`, `2.0`, and `1.12`.
- **Regression test:** `tests/oracle/oracle_smoke_test.cpp` — `[0][oracle]`.
- **Scope:** guards the major.minor only; it intentionally does not detect
  behavioural drift *within* an allowed minor release.

### SOLVE-TRIG-1 — `solve` returned `[]` for transcendental trig equations
- **Input:** `solve(sin(x), x)`, `solve(cos(x), x)`, `solve(2*sin(x)-1, x)`,
  `solve(sin(2*x), x)`, `solve(cos(3*x), x)`, `solve(tan(x)-1, x)`.
- **Was:** empty `[]` — these have infinite (periodic) solution sets, so
  `solveset` returns an `ImageSet` and the vector-returning `solve` surfaced
  nothing finite.
- **Expected (SymPy `solve`):** representative roots over one period:
  `[0, pi]`, `[pi/2, 3*pi/2]`, `[pi/6, 5*pi/6]`, `[0, pi/2]`, `[pi/6, pi/2]`,
  `[pi/4]`.
- **Fix (`src/solvers/solve.cpp`):** added `solve_trig`, tried in `solve`
  before the `solveset` fallback. It matches a single trig term
  `A*f(B*x) + C = 0` (`f ∈ {sin,cos,tan}`, `A`,`C` var-free, `B*x` linear with
  no additive phase), forms `c = -C/A`, and inverts the *inner* function to its
  principal solutions — `sin`: {asin c, π−asin c}; `cos`: {acos c, 2π−acos c};
  `tan`: {atan c} — then divides each by `B` and dedupes. This mirrors SymPy,
  which inverts and divides by `B` rather than enumerating every `x ∈ [0,2π)`.
- **Verified against SymPy:** all ten forms match `solve` as a set (root order
  is presentation-dependent). Pre-existing paths unchanged: `sinh(x)→[0]`,
  `log(x)-1→[E]`, polynomials, radicals.
- **Regression test:** `tests/solvers/solve_test.cpp`
  — `[10][solve][trig][oracle][regression]` (order-independent set comparison).
- **Scope:** a single `sin`/`cos`/`tan` term, linear argument, var-free
  coefficient — plus a homogeneous positive-integer power `f(B*x)^n = 0`, which
  reduces to `f(B*x) = 0` (so `sin(x)²→[0,π]`, `cos(x)²→[π/2,3π/2]`,
  `tan(x)²→[0]`, `2·sin(x)²→[0,π]`) — plus a **zero-product** of var-dependent
  factors, solved factor-by-factor (recursively) and unioned:
  `sin x·cos x→[0,π/2,π,3π/2]`, `sin x·(cos x−1)→[0,π,2π]`,
  `(sin x−1)(cos x+1)→[π/2,π]`. Trig *combinations* (`sin x + cos x`),
  *non-homogeneous* powers (`sin(x)²−1`), and phase-shifted arguments remain out
  of scope (decline cleanly). A mixed polynomial·trig product such as `x·sin x`
  still yields only the algebraic root `[0]` — the polynomial path resolves it
  before the trig handler runs.

### GAMMA-REC-1 — gamma recurrence `z*gamma(z) → gamma(z+1)` wasn't applied
- **Input:** `combsimp(x*gamma(x))`, `combsimp(x*(x+1)*gamma(x))`,
  `combsimp(gamma(x+1)/x)`.
- **Was:** left unevaluated — `combsimp`/`gammasimp` only handled gamma ratios,
  reflection, and binomial collapse, never absorbed polynomial factors.
- **Expected (SymPy):** `gamma(x+1)`, `gamma(x+2)`, `gamma(x)`.
- **Fix (`src/simplify/simplify.cpp`):** added `gamma_recurrence`, applied by
  both `combsimp_node` and `gammasimp_node`. Within a `Mul` it iterates to a
  fixpoint, absorbing a numerator factor equal to a gamma argument `z` (raising
  it to `gamma(z+1)` per Γ(z+1)=z·Γ(z)) or a denominator factor equal to `z-1`
  (lowering `gamma(z)` to `gamma(z-1)`). The fixpoint loop lets chains collapse:
  `x*(x+1)*gamma(x)→gamma(x+2)`. Spectator factors (`y`, `2`) are preserved, and
  `x*gamma(x+1)` is correctly left alone.
- **Verified against SymPy:** all six forms match `combsimp`; the reflection and
  ratio passes still hold.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[5][combsimp][gamma][oracle][regression]`.
- **Scope:** integer-step recurrence via exact factor matching. Non-integer
  argument shifts and rational-function denominators beyond a single `z-1` term
  are out of scope.

### BINOM-COMB-1 — `combsimp` didn't collapse `binomial(n,k)` to polynomial form
- **Input:** `combsimp(binomial(n,n))`, `combsimp(binomial(n,n-1))`,
  `combsimp(binomial(n+1,n))`, `combsimp(binomial(n,2))`, `combsimp(binomial(n,3))`.
- **Was:** all left unevaluated — `combsimp` only cancelled factorial *ratios*,
  and the `binomial` factory collapses only `k ∈ {0,1}` or literal-integer args.
- **Expected (SymPy `combsimp`):** `1`, `n`, `n+1`, `n*(n-1)/2`, `n*(n-1)*(n-2)/6`.
- **Fix (`src/simplify/simplify.cpp`):** added `combsimp_binomial`, which folds
  `binomial(n,k)` whenever `k` or `n-k` is a small non-negative integer `m`, via
  the gamma identity `binomial(n,k) = falling_factorial(n,m)/m!` (valid for
  symbolic `n`). `m = n-k` is tried first (the symmetric tail: `n`, `n-1`, …),
  then `m = k` (the small head: `0,1,2,…`). `combsimp_node` runs this after the
  factorial-ratio pass. Fully symbolic `binomial(n,k)` is left untouched.
- **Verified against SymPy:** all five forms match `combsimp`; the `binomial`
  factory's integer fast-paths (`binomial(5,2)=10`, `binomial(n,0)=1`) are intact.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[5][combsimp][binomial][oracle][regression]`.
- **Scope:** one of `k` / `n-k` a non-negative integer ≤ 50. Genuinely symbolic
  binomials and the gamma/factorial recurrence collapse (`x*gamma(x)→gamma(x+1)`)
  remain out of scope.

### GAMMA-REFL-1 — `gammasimp` missed the Euler reflection formula
- **Input:** `gammasimp(gamma(x)*gamma(1-x))`, `gammasimp(gamma(2*x)*gamma(1-2*x))`.
- **Was:** left as `gamma(-x + 1)*gamma(x)` — `gammasimp` only cancelled gamma
  *ratios* (`gamma(a)/gamma(b)` with integer `a-b`), never products.
- **Expected (SymPy `gammasimp`):** `pi/sin(pi*x)`, `pi/sin(2*pi*x)`.
- **Fix (`src/simplify/simplify.cpp`):** added `gamma_reflection`, which scans a
  `Mul` for two numerator gamma factors whose arguments sum to 1 and folds each
  pair via Γ(z)·Γ(1−z) = π/sin(πz). The surviving argument is chosen free of a
  leading additive constant so the output reads `pi/sin(pi*z)` rather than the
  equivalent `pi/sin(pi*(1-z))`. `gammasimp_node` now runs the reflection pass
  after the existing ratio pass, so both compose (`gamma(n+1)/gamma(n)` → `n`
  still holds, and a spectator factor `y` is preserved).
- **Verified against SymPy:** all three forms match `gammasimp`.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[5][gammasimp][reflection][oracle][regression]`.
- **Scope:** numerator gamma pairs with arguments summing to 1 (any shared
  scalar multiple of the argument works). Reflection of gamma *ratios* or of
  more general argument relations is out of scope.

### CONJ-FN-1 / ARG-CX-1 — conjugate over analytic functions; arg of a complex value
- **Input:** `conjugate(exp(I*x))` (x real), `conjugate(cosh(z))`;
  `arg(I)`, `arg(1+I)`, `arg(-1+I)`, `arg(-I)`.
- **Was:** `conjugate(exp(x·I))` left unevaluated; `arg(I)`, `arg(1+I)` left
  unevaluated (and `"arg"` wasn't even recognised by the parser).
- **Expected (SymPy):** `exp(-I*x)`, `cosh(conjugate(z))`; `pi/2`, `pi/4`,
  `3*pi/4`, `-pi/2`.
- **Fix (`src/functions/miscellaneous.cpp`, `src/parsing/parser.cpp`):**
  - `conjugate(f(g)) = f(conjugate(g))` for an entire function with real Taylor
    coefficients — `exp`, `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh` (`log`
    excluded for its branch cut, matching SymPy). With CONJ-DIST-1 this gives
    `conjugate(exp(I·x)) = exp(−I·x)` for real `x`.
  - `arg(z) = atan2(im z, re z)` when the real/imaginary parts resolve (free of
    unevaluated `Re`/`Im`) and the imaginary part is nonzero — `atan2` already
    evaluates the quadrant values, so `arg(I) = π/2`, `arg(1+I) = π/4`, etc.
  - the parser now maps `"arg"` to the `arg_` factory.
- **Verified against SymPy:** the conjugate-over-function family and the four
  `arg` values match; `conjugate(log(z))` correctly stays unevaluated.
- **Regression tests:** `tests/functions/miscellaneous_test.cpp`
  — `[3h][conjugate][oracle][regression]` (CONJ-FN-1) and
  `[3h][arg][oracle][regression]` (ARG-CX-1).
- **Scope:** the listed analytic functions and complex-value `arg`. `arg` of a
  symbolic `a+b·I` reduces only when `atan2(b,a)` itself has a closed form.

### ABS-MOD-1 — symbolic complex modulus `|a + b·I|` wasn't computed
- **Input:** `Abs(x + I*y)`, `Abs(2 + I*y)` (x, y real).
- **Was:** the unevaluated `Abs(x + y·I)`. The `abs` factory computed the modulus
  only for a *numeric* `a + b·I` (`rational_complex`); a symbolic one fell
  through.
- **Expected (SymPy):** `sqrt(x**2 + y**2)`, `sqrt(y**2 + 4)`.
- **Fix (`src/functions/miscellaneous.cpp`):** when `re(arg)` and `im(arg)`
  resolve to expressions free of unevaluated `Re`/`Im` nodes (now possible after
  REIM-DIST-1) and the imaginary part is nonzero, return `√(re² + im²)`. A
  generic `Abs(z)` keeps its `Re(z)`/`Im(z)` split and so stays unevaluated,
  matching SymPy.
- **Verified against SymPy:** `Abs(x+I·y) → √(x²+y²)`, `Abs(2+I·y) → √(y²+4)`;
  generic `Abs(z)` and real `Abs(x)` are unchanged.
- **Regression test:** `tests/functions/miscellaneous_test.cpp`
  — `[3d][abs][oracle][regression]` (ABS-MOD-1).
- **Scope:** complex numbers whose parts are real-determinable. `|exp(I·x)| = 1`
  (modulus of a transcendental imaginary) and `arg(...)` of a complex value are
  follow-ups.

### REIM-DIST-1 — `re`/`im` didn't distribute or handle `I`
- **Input:** `re(I*x)`, `im(I*x)`, `re(x+I*y)`, `im(x+I*y)` (x, y real);
  `re(I*z)` (z generic).
- **Was:** the unevaluated `re(x*I)`, `im(x*I)`, `re(x + y*I)`, … . The `re`/`im`
  factories handled only a real argument, a numeric `a+b·I`, and a leading
  negative factor — they didn't distribute over a sum or recognise the imaginary
  unit.
- **Expected (SymPy):** `0`, `x`, `x`, `y`; `re(I*z) = -im(z)`.
- **Fix (`src/functions/miscellaneous.cpp`):**
  - `re`/`im` are linear over `Add` (`re(Σ aᵢ) = Σ re(aᵢ)`);
  - a Mul is decomposed as `c · Iᵏ · w` (`decompose_mul_complex`: real factors →
    `c`, I-count mod 4 → `k`, the rest → `w`), and the real coefficient is pulled
    out with the `Iᵏ` rotation: `re(c·w)=c·re(w)`, `re(c·I·w)=−c·im(w)`,
    `im(c·I·w)=c·re(w)`, … . The decomposition returns nothing unless a real
    factor or an `I` was peeled off, so the recursion terminates.
- **Verified against SymPy:** the real-symbol cases collapse exactly
  (`re(I·x)=0`, `im(I·x)=x`, `re(x+I·y)=x`, `im(x+I·y)=y`, `re(2x)=2x`), and the
  generic I-rotation `re(I·z)=−im(z)`, `im(I·z)=re(z)` matches.
- **Regression test:** `tests/functions/miscellaneous_test.cpp`
  — `[3h][complex][oracle][regression]` (REIM-DIST-1).
- **Scope:** linearity + the imaginary-unit rotation. `re`/`im` of a generic
  *product* of two unknown-reality factors stays unevaluated.

### CONJ-DIST-1 / ABS-I-1 — `conjugate` didn't distribute; `Abs(I·x)` not reduced
- **Input:** `conjugate(I*x)`, `conjugate(x*y)`, `conjugate(x+y)`,
  `conjugate(x**2)`; `Abs(I*x)`, `Abs(2*I*x)`.
- **Was:** `conjugate(x*I)`, `conjugate(I*x)` left unevaluated; `Abs(x*I)`
  unreduced. The `conjugate` factory handled only real / numeric-complex /
  involution; the `Abs` Mul-pull-out pulled out numeric and known-sign factors
  but not the imaginary unit.
- **Expected (SymPy):** `-I*conjugate(x)`, `conjugate(x)*conjugate(y)`,
  `conjugate(x)+conjugate(y)`, `conjugate(x)**2`; `Abs(x)`, `2*Abs(x)`.
- **Fix (`src/functions/miscellaneous.cpp`):**
  - `conjugate` now distributes over `Mul`, `Add`, and integer `Pow` (it is a
    ring homomorphism); recursion reduces each part (`conjugate(I) = −I`,
    `conjugate(real) = real`), so `conjugate(I·x) = −I·conjugate(x)`.
  - `abs` pulls the imaginary unit out of a product (`|I| = 1 ⇒ |I·x| = |x|`)
    alongside the numeric / known-sign factors.
- **Verified against SymPy:** all the distribution cases and `Abs(I·x)`,
  `Abs(2·I·x)`, `Abs(I·x·y)` match; existing `Abs` reductions are unchanged.
- **Regression tests:** `tests/functions/miscellaneous_test.cpp`
  — `[3h][conjugate][oracle][regression]` (CONJ-DIST-1) and
  `[3d][abs][oracle][regression]` (ABS-I-1).
- **Scope:** conjugate distribution and the imaginary-unit `Abs` reduction.

### ASSUME-IMAG-1 — no Imaginary / Complex assumption predicates
- **Was:** the assumption vocabulary had no `Imaginary` or `Complex` key, so
  `I.is_imaginary`, `is_real(I·x)`, `is_complex(x)` had no answer — SymPP modelled
  complex values structurally (`a + b·I`) but couldn't *reason* about
  imaginariness.
- **Fix:** added `Complex` and `Imaginary` to `AssumptionKey`, mask fields +
  builders (`set_complex`/`set_imaginary`), hash, and deductive closure
  (`imaginary ⇒ complex ∧ ¬real ∧ finite ∧ nonzero ∧ ¬rational/integer/sign/
  parity`; `real ⇒ complex ∧ ¬imaginary`; `zero ⇒ ¬imaginary` since 0 is real).
  Wired:
  - `ImaginaryUnit::ask` — `I` is imaginary, complex, finite, ¬real;
  - generic `ask()` derivations — `real ∨ imaginary ⇒ complex`,
    `real ⇒ ¬imaginary`, `imaginary ⇒ ¬real`;
  - `Mul::ask` — an odd number of imaginary factors (rest real, all nonzero) is
    imaginary (`I·real = imaginary`), an even number is real (`I·I = −1`);
  - `Add::ask` — a sum of imaginaries is imaginary, a real + imaginary mix is
    complex but neither;
  - `Pow::ask` — `imaginary^(odd integer)` is imaginary, `^(even)` is real.
  - the infinities answer `¬complex` (∞/zoo aren't finite complex numbers).
- **Verified against SymPy:** `I`, `2·I`, `I·x` (x real ≠ 0), `xi` (declared
  imaginary), `xi²`, `xi³`, `I·I`, `x·y` (reals), `x_r + I·y_r`, plain reals and
  `0` — the `is_imaginary` / `is_real` / `is_complex` triples match
  `sympy`'s on 9/10 (the 10th, `x_r + I·x_r` → `is_real`, is conservatively
  `Unknown` in SymPP vs `False` in SymPy — proving ¬real needs imaginary-part
  cancellation analysis).
- **Regression test:** `tests/core/assumptions_test.cpp`
  — `[2][assumptions][imaginary][regression]` (ASSUME-IMAG-1).
- **Scope:** the imaginary/complex ontology + arithmetic propagation. Still
  deferred: irrational/prime/algebraic/hermitian/commutative predicates and the
  SAT-based `ask(query, assumptions)` reasoner. ~14 of SymPy's ~30+ predicates.

### REWRITE-EXP-1 — no `rewrite(target)` API (exp ↔ trig)
- **Was:** SymPP had no analogue of SymPy's `expr.rewrite(target)` — a common
  cross-cutting operation (Euler / hyperbolic identities, used in solving and
  simplification).
- **Fix (`src/simplify/simplify.cpp`, `include/sympp/simplify/simplify.hpp`):**
  new `rewrite(expr, "exp")` re-expresses `sin`/`cos`/`tan` and
  `sinh`/`cosh`/`tanh` as exponentials (`sin(x) → −i·(e^{ix}−e^{−ix})/2`,
  `cosh(x) → (e^x+e^{−x})/2`, …), applied recursively so combinations and
  composite arguments (`sin(2x)`) are handled. An unknown target is a no-op.
- **Verified against SymPy:** the six trig/hyperbolic forms plus `sin(x)+cos(x)`
  and `sin(2x)` all equal `expr.rewrite(exp)` symbolically.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[5][rewrite][oracle][regression]` (REWRITE-EXP-1).
- **Scope:** target `"exp"`. Other targets (`rewrite(exp, sin/cos)` Euler
  direction, `rewrite(tan, …)`, gamma/factorial cross-rewrites) are follow-ups.

### SOLVE-DEDUP-1 — `solve` returned duplicate roots for repeated factors
- **Input:** `solve((x+2)**2)`, `solve(x**2*(x-1))`,
  `solve((x-1)**2*(x+1))`, `solve((x-1)**3)`.
- **Was:** `[-2, -2]`, `[0, 0, 1]`, `[-1, 1, 1]`, `[1, 1, 1]` — `solve_poly`
  emits a root once per (square-free) factor, so a repeated factor produced
  duplicates. (Surfaced after SOLVE/INEQ's `expand` made factored polynomials
  reachable.)
- **Expected (SymPy):** `[-2]`, `[0, 1]`, `[-1, 1]`, `[1]` — SymPy's `solve`
  returns the distinct solution set.
- **Fix (`src/solvers/solve.cpp`):** collapse structurally-equal roots in
  `solve`, preserving order.
- **Verified against SymPy:** the repeated-factor cases now return the distinct
  set; genuinely distinct roots (`(x-1)(x-2)(x-3) → {1,2,3}`) are unchanged.
- **Regression test:** `tests/solvers/solve_test.cpp`
  — `[10][solve][regression]` (SOLVE-DEDUP-1).
- **Scope:** root-set deduplication. Multiplicity is not reported (SymPy's
  `solve` default also drops it; `roots()` keeps it — not implemented).

### INEQ-EXACT-1 — inequalities used float endpoints + a 1e30 ∞ proxy; `solve` ignored factored polynomials
- **Input:** `solve_univariate_inequality(x²−4 < 0)`, `x²−4 > 0`, `x²+1 > 0`,
  `(x−1)(x−2) < 0`; and `solve((x−1)*(x−2))`.
- **Was:**
  - endpoints came back as `Float`s (`−2.0000…`, `2.0000…`) instead of exact
    integers, because every root was round-tripped through a `double`;
  - the unbounded ends used a literal `1e30` as an `∞` proxy
    (`x²+1>0 → (−1e30, 1e30)` instead of `Reals`), with a code comment noting
    *"we don't have an Infinity singleton"* — stale since the infinity atoms
    shipped;
  - `solve((x−1)*(x−2))` returned `[]` (the `Poly` machinery couldn't build from
    the un-expanded `Mul`), so `(x−1)(x−2) < 0` wrongly gave `EmptySet`.
- **Expected (SymPy):** `(−2, 2)`, `(−∞,−2) ∪ (2,∞)`, `Reals`, `(1, 2)`,
  `[1, 2]`.
- **Fix (`src/solvers/solve.cpp`):**
  - `solve` now `expand`s the input before the polynomial path, so a factored
    polynomial is solved (`(x−1)(x−2) → [1, 2]`);
  - `solve_univariate_inequality` keeps each root as its **exact** `Expr` (paired
    with a `double` only for ordering / sign-sampling), emits the real
    `S::Infinity()` / `S::NegativeInfinity()` at the unbounded ends, and returns
    `reals()` when there are no roots and the sign matches.
- **Verified against SymPy:** the family above matches exactly, including the
  closed-endpoint `≤` case and the ray `Union`.
- **Regression test:** `tests/solvers/solve_test.cpp`
  — `[10][inequality][regression]` (INEQ-EXACT-1).
- **Scope:** real univariate inequalities with numeric roots. Irrational roots
  still order via their numeric value but appear exactly in the endpoints.

### SET-COMPL-1 — `ℝ \ interval` wasn't computed, and ray membership was Unknown
- **Input:** `set_complement(Reals, Interval(1,3))`,
  `set_complement(Reals, Interval.open(1,3))`,
  `set_complement(Reals, Interval(1, oo))`; and membership queries on the result.
- **Was:** an unevaluated `Complement(Reals, [1,3])` node. Even when a ray was
  built by hand, `Interval::contains` bailed because a ±∞ endpoint is not a
  number (`is_number(oo) == false`), so every membership test on `(−∞,1)`
  returned Unknown.
- **Expected (SymPy):** `(−∞,1) ∪ (3,∞)`, `(−∞,1] ∪ [3,∞)`, `(−∞,1)`; and
  `0 ∈`, `2 ∉` for the first.
- **Fix (`src/sets/sets.cpp`):**
  - `set_complement(Reals, [a,b])` builds `(−∞,a) ∪ (b,∞)` with each boundary's
    open/closed flag flipped (a point removed from ℝ is excluded from the
    complement); a ±∞ endpoint drops that ray, and `ℝ \ ℝ → EmptySet`.
  - `Interval::contains` now treats a ±∞ endpoint as an always-satisfied
    unbounded side, so membership on a ray (and hence on the complement) is
    decided.
- **Verified against SymPy:** `ℝ\[1,3]`, `ℝ\(1,3)` (endpoints flip),
  `ℝ\[1,∞)` and `ℝ\(−∞,3]` (single ray) all match; membership `0∈`, `1∉`, `2∉`,
  `4∈` is now decided.
- **Regression test:** `tests/sets/sets_test.cpp`
  — `[10][sets][complement][interval][regression]` (SET-COMPL-1).
- **Scope:** `ℝ \ interval`. Complement of a `FiniteSet` or within a bounded
  universal set is a follow-up.

### SET-INTERVAL-1 — interval `∪` / `∩` weren't computed
- **Input:** `set_union(Interval(1,3), Interval(2,4))`,
  `set_intersection(Interval(1,3), Interval(2,4))`,
  `set_intersection(Interval(1,2), Interval(3,4))`.
- **Was:** the operands wrapped in an unevaluated `Union` / `Intersection` node
  (`[1,3] ∪ [2,4]`, `[1,3] ∩ [2,4]`). `set_union`/`set_intersection` only folded
  the empty-set cases.
- **Expected (SymPy):** `Interval(1, 4)`, `Interval(2, 3)`, `EmptySet`.
- **Fix (`src/sets/sets.cpp`):** for two real intervals,
  - **intersection** = `[max(los), min(his)]` with the open flags carried from
    the winning endpoint (OR'd on a tie); `lo > hi → EmptySet`, `lo == hi →` a
    single-point `FiniteSet` (or `EmptySet` if either endpoint there is open);
  - **union** merges when the intervals overlap or touch
    (`ib.lo ≤ ia.hi ∧ ia.lo ≤ ib.hi`) into `[min(los), max(his)]`, otherwise
    stays a `Union`.
  Endpoint ordering uses a sign comparison (`endpoint_cmp`), so symbolic bounds
  that can't be ordered fall back to the unevaluated node.
- **Verified against SymPy:** overlap/adjacent merge (`[1,3]∪[2,4]=[1,4]`,
  `[1,2]∪[2,3]=[1,3]`), disjoint union stays a `Union`, intersection
  (`[1,3]∩[2,4]=[2,3]`), containment (`[1,5]∩[2,3]=[2,3]`), disjoint → `EmptySet`,
  closed-touch → `{3}`, open-touch → `EmptySet` — all match.
- **Regression test:** `tests/sets/sets_test.cpp`
  — `[10][sets][interval][regression]` (SET-INTERVAL-1).
- **Scope:** pairs of real intervals with orderable endpoints. Multi-set unions,
  interval-vs-FiniteSet, and `Complement(Reals, …) → ray ∪ ray` are follow-ups.

### ILAPLACE-QUAD-2 — inverse Laplace of a LINEAR numerator over a quadratic
- **Input:** `inverse_laplace_transform(s/(s**2+2*s+2))`,
  `(s+1)/(s**2+2*s+2)`, `s/((s-2)**2+1)`, `(2*s+1)/(s**2+2*s+5)`.
- **Was:** the unevaluated marker. ILAPLACE-QUAD-1 closed the *constant*-numerator
  case, but a numerator linear in `s` (the damped-cosine family) still fell
  through — `inverse_laplace_term` bails as soon as the numerator contains `s`.
- **Expected (SymPy):** `exp(-t)·cos t − exp(-t)·sin t`, `exp(-t)·cos t`,
  `exp(2t)·cos t + 2·exp(2t)·sin t`, `2·exp(-t)·cos 2t − exp(-t)·sin 2t/2`.
- **Fix (`src/integrals/transforms.cpp`):** new `inverse_laplace_linear_quad` —
  split `F = num·den^(-1)`, require `num` linear and `den` an irreducible
  quadratic (`Poly` degrees 1 and 2, `β' ≠ 0`); complete the square and
  decompose the numerator about `(s − a)`:
  `(α·s+β)/α' = A·(s−a) + B` with `A = α/α'`, `B = (β + α·a)/α'`, giving
  `A·exp(a·t)·cos(b·t) + (B/b)·exp(a·t)·sin(b·t)`.
- **Verified against SymPy:** all four inputs match; the constant-numerator
  (ILAPLACE-QUAD-1) and pure `s/(s²+a²)` paths are unaffected (the new handler
  requires a genuine linear numerator and `β' ≠ 0`).
- **Regression test:** `tests/integrals/transforms_test.cpp`
  — `[8][inverse_laplace][oracle][regression]` (ILAPLACE-QUAD-2).
- **Scope:** numerators up to degree 1 over an irreducible quadratic. Higher-
  degree rational functions still rely on `apart` to split first.

### ILAPLACE-QUAD-1 — inverse Laplace couldn't invert a completed-square quadratic
- **Input:** `inverse_laplace_transform(1/(s**2+2*s+2))`,
  `2/(s**2+4*s+5)`, `1/(s**2+2*s+10)`.
- **Was:** the unevaluated `InverseLaplaceTransform(...)` marker. The inverse
  table handled `(s−a)^n` and `c/(s²+a²)` (no linear term), but a quadratic with
  a **linear** term — `s²+ps+q` with discriminant `< 0` — matched neither.
- **Expected (SymPy):** `exp(-t)·sin(t)`, `2·exp(-2t)·sin(t)`,
  `exp(-t)·sin(3t)/3` (the inverse s-shift, symmetric to LAPLACE-SHIFT-1).
- **Fix (`src/integrals/transforms.cpp`):** a Case 3 in `inverse_laplace_term` —
  build a `Poly` in `s`; for a degree-2 denominator `α·s²+β·s+γ` with `β ≠ 0` and
  `b² = γ/α − (β/2α)² > 0`, complete the square (`a = −β/(2α)`) and return
  `(num/α)·exp(a·t)·sin(b·t)/b`.
- **Verified against SymPy:** the completed-square family matches (modulo the
  `Heaviside(t)` factor SymPP omits); the existing `(s−a)^n` and `s²+a²` paths
  are unchanged (β = 0 still routes to the plain `sin` case).
- **Regression test:** `tests/integrals/transforms_test.cpp`
  — `[8][inverse_laplace][oracle][regression]` (ILAPLACE-QUAD-1).
- **Scope:** constant numerator over an irreducible quadratic. A linear
  numerator (`(s−a)/((s−a)²+b²) → exp(a·t)·cos(b·t)`) over the shifted quadratic
  is the symmetric follow-up; and the `1/(s²−a²) → sinh/cosh` inverse still
  prints via the complex `−I·sin(i·t)` form rather than `sinh`.

### LAPLACE-SHIFT-1 — Laplace transform missed `sinh`/`cosh` and the s-shift theorem
- **Input:** `laplace_transform(sinh(t))`, `cosh(t)`, `exp(-t)·sin(t)`,
  `t·exp(t)`, `t²·exp(t)`, `exp(2t)·cos(t)`.
- **Was:** the unevaluated `LaplaceTransform(...)` marker. The table covered
  `t^n`, `exp`, `sin`, `cos`, but not the hyperbolics, and the `Mul` case only
  pulled out constant factors — so any `exp(a·t)·g(t)` product (damped
  oscillations, `t^n·exp`) fell through.
- **Expected (SymPy):** `1/(s²−1)`, `s/(s²−1)`, `1/((s+1)²+1)`, `1/(s−1)²`,
  `2/(s−1)³`, `(s−2)/((s−2)²+1)`.
- **Fix (`src/integrals/transforms.cpp`):**
  - `sinh`/`cosh` table entries: `L{sinh(a·t)} = a/(s²−a²)`,
    `L{cosh(a·t)} = s/(s²−a²)`.
  - the **s-shift theorem** in the `Mul` case: every `exp(a·t)` factor is pulled
    out, the `a`'s summed, and the rest's transform `G(s)` is shifted to
    `G(s − a)` — closing the damped-oscillation and `t^n·exp` families.
- **Verified against SymPy:** all six inputs match, including the scaled
  `3·exp(−2t)·sin(3t) → 9/((s+2)²+9)`; the existing `t`, `sin`, `cos`, `exp`,
  linearity entries are unchanged.
- **Regression test:** `tests/integrals/transforms_test.cpp`
  — `[8][laplace][oracle][regression]` (LAPLACE-SHIFT-1).
- **Scope:** the table + s-shift. The general Meijer-G-driven transform of
  arbitrary inputs stays deferred (it depends on the hypergeometric machinery).

### SIMP-EXP-POW-1 — `simplify((exp(x))**2)` didn't fold to `exp(2x)`
- **Input:** `simplify(exp(x)**2)`, `exp(x)**3`, `exp(x)**(-1)`, `exp(x+1)**2`.
- **Was:** unchanged (`exp(x)**2`, …). `combine_exp` merged `exp` factors inside
  a `Mul`, but a standalone `Pow(exp(g), k)` was never folded.
- **Expected (SymPy):** `exp(2*x)`, `exp(3*x)`, `exp(-x)`, `exp(2*x + 2)`.
- **Fix (`src/simplify/simplify.cpp`):** `combine_exp_node` now folds a
  `Pow(exp(g), k)` with an **integer** `k` to `exp(expand(k·g))`. A fractional or
  symbolic exponent is left as a `Pow` — matching SymPy, which keeps
  `sqrt(exp(x))` and `exp(x)**n` for branch-cut safety.
- **Verified against SymPy:** `exp(x)**{2,3,-1}` and `exp(x+1)**2` fold exactly;
  `exp(x)**(1/2)` (≡ SymPy's `sqrt(exp(x))`) and `exp(x)**n` are left unfolded.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[5][simplify][oracle][regression]` (SIMP-EXP-POW-1).
- **Scope:** integer power of a single `exp`. A power of a *product* of exps
  (`(exp(x)·exp(y))**2`) needs a second combine pass and is left as-is (still
  correct, just not maximally combined).

### LIMIT-HANG-1 — `limit` hung on a radical `∞/∞` form
- **Input:** `limit(sqrt(x**2+x) - x, x, oo)`,
  `limit(x/(sqrt(x**2+x)+x), x, oo)`.
- **Was:** an effectively-infinite hang (CPU spin, no result). `lhopital_nd`
  differentiates num/den each step; for a radical integrand the nested radicals
  grow every iteration (the ratio never stabilises), so the per-iteration
  `simplify`/`together`/`diff` over an ever-larger expression dominated the
  runtime — a CAS that locks up on a finite input.
- **Fix (`src/calculus/limit.cpp`):** a node-count budget in the `lhopital_nd`
  loop — when `node_count(num) + node_count(den)` exceeds 400, bail to the
  unevaluated `nan`. Legitimate limits resolve in a handful of iterations far
  under the budget, so none are affected.
- **Verified against SymPy:** the radical inputs now **terminate** (returning
  `nan`) instead of hanging, and the neighbouring limits still resolve
  (`sin(x)/x → 1`, `x·e^(-x) → 0`, `(1+1/x)^x → E`, `(x²-1)/(x-1) → 2`).
- **Regression test:** `tests/calculus/series_limit_test.cpp`
  — `[6][limit][infinity][regression]` (LIMIT-HANG-1): the radical forms return
  *a* value (no hang).
- **Scope:** this is a **robustness** fix — it stops the hang but does not
  compute the limit. `sqrt(x²+x) − x = 1/2` needs asymptotic-series / Gruntz
  machinery (the `x = 1/t` substitution that turns the ∞−∞ into a 0/0 at the
  origin), which stays deferred-deep.

### SIMP-EXP-HYP-1 — `simplify` didn't fold `e^x + e^(−x)` into `2·cosh(x)`
- **Input:** `simplify(exp(x) + exp(-x))`, `simplify(exp(x) - exp(-x))`,
  `simplify(exp(2x) + exp(-2x))`, `simplify(3·exp(x) + 3·exp(-x))`.
- **Was:** unchanged (`exp(-x) + exp(x)`, …). SymPP had the cosh/sinh → exp
  direction (TRIG-HYP-2) but not the reverse, so an exponential sum never
  collapsed to a hyperbolic function.
- **Expected (SymPy):** `2·cosh(x)`, `2·sinh(x)`, `2·cosh(2x)`, `6·cosh(x)`.
- **Fix (`src/simplify/simplify.cpp`):** a new `exp_to_hyp_add` pass (mirror of
  `hyp_to_exp_add`) collects, per argument `g`, the coefficients of `e^g` and
  `e^(−g)`; equal coefficients fold to `2a·cosh(g)`, opposite to `2a·sinh(g)`.
  The argument is normalised to its positive representative (`cosh` even, `sinh`
  odd) so the output matches SymPy's `2·cosh(2x)` rather than `2·cosh(−2x)`.
  Wired into the `simplify` pipeline after `combine_exp`. No oscillation with
  TRIG-HYP-2: a pure `2·cosh(x)` has no `sinh` partner to convert back.
- **Verified against SymPy:** the cosh/sinh folds for arguments `x` and `2x`
  with integer coefficients all match exactly; an unequal-coefficient sum
  (`e^x + 2·e^(−x)`) is correctly left alone.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[5][simplify][oracle][regression]` (SIMP-EXP-HYP-1).
- **Scope:** real exponentials with equal/opposite coefficients per argument.

### SOLVESET-TRIG-SCALE-1 — `solveset(cos(2x)=1)` returned EmptySet; redundant cos union
- **Input:** `solveset(cos(2*x) - 1, x)`, `solveset(cos(x) - 1, x)`,
  `solveset(cos(x) + 1, x)`.
- **Was:** `cos(2x) − 1` → `EmptySet` (wrong — it has solutions `{nπ}`);
  `cos(x) − 1` → `ImageSet ∪ ImageSet` of two **identical** `{2nπ}` sets. The
  `invert_solveset` trig branches only emitted an ImageSet when the argument was
  exactly `var`, so a scaled argument `a·var` fell through; and the cos branch
  always emitted a two-branch `±acos(c)` union even when the branches coincide.
- **Expected (SymPy):** `{nπ}`, `{2nπ}`, `{2nπ + π}`.
- **Fix (`src/solvers/solve.cpp`):** the Sin/Cos/Tan branches now accept a linear
  argument `g = a·var` (a `linear_coeff` helper) and divide the periodic image
  through by `a`. The cos branch emits a **single** ImageSet when
  `acos(c) ∈ {0, π}` (`c = ±1`, where the `±` branches coincide), and the
  two-branch union otherwise.
- **Verified against SymPy:** `cos(2x)=1 → {nπ}`, `cos(x)=1 → {2nπ}`,
  `cos(x)=-1 → {2nπ+π}` match exactly; `sin(2x)`, `tan(2x)=1`, `cos(3x)=1/2`,
  `cos(x)=1/2` equal SymPy's solution sets (SymPP often emits the cleaner single
  ImageSet where SymPy emits an equivalent union); a generic RHS keeps the
  two-branch union.
- **Regression test:** `tests/solvers/solve_test.cpp`
  — `[10][solveset][regression]` (SOLVESET-TRIG-SCALE-1).
- **Scope:** linear (a·var) trig arguments. Genuinely nested arguments
  (`cos(x²)`) still need parametric back-substitution.

### SOLVESET-POW0-1 — `solveset(f(x)**n)` returned EmptySet for a transcendental base
- **Input:** `solveset(sin(x)**2, x)`, `solveset(sin(x)**4, x)`,
  `solveset(tan(x)**2, x)`.
- **Was:** `EmptySet` — clearly wrong (`sin(x)² = 0` has the solutions of
  `sin(x) = 0`). `invert_solveset` only peeled a `Pow` with a *non-integer*
  exponent (SOLVE-RAD-1); an integer power fell through to the polynomial path,
  which can't build a `Poly` in `x` from `sin(x)` and so returned no solutions.
- **Expected (SymPy):** the solution set of the base — `{n·π}` for `sin(x)²`,
  `tan(x)²`; `{(2n+1)π/2}` for `cos(x)²`.
- **Fix (`src/solvers/solve.cpp`):** in `invert_solveset`'s `Pow` branch, when the
  exponent is a positive integer and the RHS is `0`, recurse with
  `solveset(base)` — `g^n = 0 ⟺ g = 0`.
- **Verified against SymPy:** `sin(x)²`, `sin(x)⁴`, `cos(x)²`, `tan(x)²` all now
  return the (periodic) solution set instead of EmptySet, and equal SymPy's set.
  (SymPP emits the cleaner single ImageSet `{n·π}` where SymPy emits the
  equivalent union `{2nπ} ∪ {2nπ+π}`.) A polynomial base (`(x-1)² → {1}`) is
  unaffected.
- **Regression test:** `tests/solvers/solve_test.cpp`
  — `[10][solveset][regression]` (SOLVESET-POW0-1).
- **Scope:** `g^n = 0`. A non-zero RHS with a scaled trig argument
  (`cos(2x) = 1`) still needs the scaled-argument trig inversion and remains a
  follow-up; `rewrite(target)` (exp↔trig, etc.) is not implemented at all.

### EXPAND-TRIG-MULTI-1 — `expand_trig` didn't expand multiple angles `sin(n·x)`
- **Input:** `expand_trig(sin(2*x))`, `expand_trig(cos(2*x))`,
  `expand_trig(sin(3*x))`, `expand_trig(cos(3*x))`, `expand_trig(sin(4*x))`.
- **Was:** the argument unchanged (`sin(2*x)`, …). `expand_trig_node` applied the
  angle-addition formula only when the argument was an `Add` (`sin(x+y)`); a
  multiple angle `n·x` is a `Mul`, so it fell through.
- **Expected (SymPy):** `2·sin(x)·cos(x)`, `2·cos²x − 1`, `3·sin x − 4·sin³x`,
  `4·cos³x − 3·cos x`, …
- **Fix (`src/simplify/simplify.cpp`):** `expand_trig_node` now also splits a
  `Mul` argument with an integer factor `n ≥ 2` as `n·g = g + (n−1)·g` and
  applies the same angle-addition formula; the fixpoint loop (raised to 32
  passes) reduces `(n−1)·g` recursively down to `sin(x)`/`cos(x)`. Works for
  `sin`/`cos`/`tan`, and composes with the `Add` case (`cos(2x+y)`).
- **Verified against SymPy:** `sin/cos/tan(n·x)` for n = 2…4 and the combined
  `cos(2x+y)` all match `sympy.expand_trig` **up to trig-identity equivalence**
  (SymPP keeps the `cos²−sin²`/nested-product form; SymPy applies a final
  Chebyshev normalization to the minimal all-sin / all-cos form — the two are
  equal, just shaped differently).
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[5][expand_trig][oracle][regression]` (EXPAND-TRIG-MULTI-1).
- **Scope:** the expansion is correct but not minimal; collapsing to SymPy's
  Chebyshev-reduced form (all-sin for `sin`, all-cos for `cos`) is a further
  normalization. Downstream `fu`/`simplify` re-`trigsimp` the result, so the
  verbose intermediate doesn't leak into their output.

### LIMIT-POLY-INF-1 — polynomial `∞−∞` limits returned `nan`
- **Input:** `limit(x**2 - x, x, oo)`, `limit(x - x**2, x, oo)`,
  `limit(2*x**2 - 5*x, x, oo)`, `limit(-x**3 + x, x, -oo)`.
- **Was:** `nan`. Direct substitution gave `∞ − ∞`, and (after LIMIT-EXP-1's
  Add-linearity, which bails when a term diverges) nothing recovered the
  dominant term.
- **Expected (SymPy):** `oo`, `-oo`, `oo`, `oo` — a polynomial at ±∞ is governed
  by its leading term.
- **Fix (`src/calculus/limit.cpp`):** when direct substitution at an infinite
  target is `nan` and the expression is a polynomial in `var` (all `Poly`
  coefficients free of `var`), take the limit of the leading term
  `c·var^deg` — `subs` then folds it through the infinity arithmetic with the
  correct even/odd-degree sign at `−∞`.
- **Verified against SymPy:** the polynomial family at both `+∞` and `−∞`
  (signs correct for even and odd leading degree). Non-polynomial `∞−∞`
  (`e^x − x`, `x − log x`) correctly **stays `nan`** — it needs the dominant-term
  / Gruntz asymptotics that remain deferred.
- **Regression test:** `tests/calculus/series_limit_test.cpp`
  — `[6][limit][infinity][regression]` (LIMIT-POLY-INF-1).
- **Scope:** polynomials. Mixed exponential/logarithmic `∞−∞` still needs Gruntz.

### FACTOR-HOM-1 — `factor` left common multivariate (homogeneous bivariate) polynomials unfactored
- **Input:** `factor(x**2 - y**2, x)`, `factor(x**2 + 2*x*y + y**2, x)`,
  `factor(x**3 - y**3, x)`, `factor(x**2*y - y**3, x)`.
- **Was:** the input unchanged. `factor` builds a `Poly` in `var`; a genuinely
  multivariate polynomial has symbolic (polynomial-in-the-other-variable)
  coefficients, which the ℚ-only `factor_list` (square-free + rational-root +
  Kronecker) can't take, so it returned the expression unfactored.
- **Expected (SymPy):** `(x - y)*(x + y)`, `(x + y)**2`,
  `(x - y)*(x**2 + x*y + y**2)`, `y*(x - y)*(x + y)`.
- **Fix (`src/polys/poly.cpp`):** a `factor_homogeneous_bivariate` pass for the
  two-symbol case. When every monomial shares the same total degree, the
  polynomial is **dehomogenized** (other variable → 1), factored over ℚ with the
  existing machinery, and each factor **re-homogenized** to its own degree
  (`Σ aₖ·xᵏ ↦ Σ aₖ·xᵏ·yⁿᵈᵉᵍ⁻ᵏ`), with a `y^(n−deg)` factor for any pure-`y`
  part. The product is **verified to expand back to the input**, so a
  non-homogeneous or irreducible polynomial is rejected rather than mis-factored.
- **Verified against SymPy:** difference of squares/cubes, sum of cubes,
  perfect-square trinomials, `9x²−4y² → (3x−2y)(3x+2y)`, `x⁴−y⁴`, the pure-`y`
  pull-out; `x²+y²` correctly stays irreducible, and univariate factoring is
  unchanged.
- **Regression test:** `tests/polys/poly_test.cpp`
  — `[4][poly][factor][oracle][regression]` (FACTOR-HOM-1).
- **Scope:** homogeneous **bivariate** polynomials (the common textbook cases).
  Non-homogeneous multivariate (`x²−y²+x`, three or more variables) still needs
  the full Wang / multivariate-GCD port and is rejected by the self-check.

### LIMIT-EXP-1 / INT-DEF-1 — `0·∞` limits with a decaying exponential, and improper definite integrals
- **Input:** `limit(x*exp(-x), x, oo)` (and `x²·e^(-x)`, …); the definite
  integrals `∫₀^∞ x^n·e^(-x) dx`.
- **Was:** `nan`. Two compounding defects:
  1. `limit` recast `x·e^(-x)` (an `∞·0` product) as the **0/0** form
     `e^(-x)/(1/x)`, where each L'Hôpital step only *raises* the polynomial
     degree, so it never converged → `nan`. It also had no linearity over `Add`
     or `Mul`, so a sum/product of such terms (the shape of these
     antiderivatives) stayed `nan`.
  2. Definite integration was literal Newton–Leibniz (`subs(F, var, oo)`), so an
     infinite bound substituted `oo` into `-(x+1)·e^(-x)` and got `∞·0 = nan`
     instead of the boundary *limit*.
- **Expected (SymPy):** `limit(x^n·e^(-x), oo) = 0`; `∫₀^∞ x^n·e^(-x) dx = n!`.
- **Fix:**
  - `src/calculus/limit.cpp`: `try_product_form` now tries **both** the 0/0 and
    `∞/∞` arrangements (the latter, `x^n / e^x`, is the one L'Hôpital cracks),
    with an **exp-aware reciprocal** (`1/exp(g) = exp(−g)`) so the exponential
    stays in the denominator across iterations instead of flipping back into the
    numerator. `limit_impl` gained **linearity over `Add` and `Mul`**: when every
    term/factor has a determinate limit (and there is no `∞−∞` / `0·∞` conflict)
    the result is their sum/product; a genuinely divergent term makes it fall
    through rather than guess.
  - `src/integrals/integrate.cpp`: the definite integral evaluates each boundary
    with `limit(antider, var, bound)` for an infinite bound (or when `subs` lands
    on `nan` / an infinity), and plain substitution otherwise.
- **Verified against SymPy:** `x^n·e^(-x) → 0` for n up to 5, `x·e^(-2x) → 0`,
  the Gamma integrals `∫₀^∞ x^n·e^(-x) = {1,2,6,24}`, `∫₀^∞ x·e^(-2x) = 1/4`;
  finite-bound integrals and convergent sums (`e^(-x) − e^(-2x) → 0`,
  `x + 1/x → oo`) unchanged.
- **Regression tests:** `tests/calculus/series_limit_test.cpp`
  (`[6][limit][infinity][regression]`, LIMIT-EXP-1) and
  `tests/integrals/integrate_test.cpp` (`[7][integrate][definite][regression]`,
  INT-DEF-1).
- **Scope:** `0·∞` where an exponential dominates a polynomial. True `∞−∞`
  forms (`x² − x`, `e^x − x`) still return `nan` — they need dominant-term /
  Gruntz asymptotics and are deliberately left rather than mis-evaluated.

### TOGETHER-LCM-1 — `together` combined fractions over the product, not the LCM, of denominators
- **Input:** `together(a/b + c/b)`, `together(x/(x+1) + 1/(x+1))`,
  `together(1/(x-1) + 1/(x-1)**2)`.
- **Was:** `(a*b + b*c)*b**(-2)`, `(x + x*(x+1) + 1)*(x+1)**(-2)` (which is
  actually `1`), and `(...)*(x-1)**(-3)`. `as_numer_denom`'s `Add` branch used
  the **product** of the part denominators as the common denominator and
  cross-multiplied, so a shared factor was squared (`b·b`, `(x+1)²`) and the
  result was left inflated and unreduced.
- **Expected (SymPy):** `(a + c)/b`, `1`, `x/(x-1)**2`.
- **Fix (`src/polys/poly.cpp`):** `as_numer_denom` now combines a sum of
  fractions over the **LCM** of the denominators. Each denominator is decomposed
  into `base^power` factors (`accumulate_denom_factors`, peeling `Pow(base,+int)`
  and `Mul`, treating anything else — Symbol, `(x+1)`, non-integer power — as an
  opaque base); the common denominator takes the max power per base, and each
  numerator is scaled by its per-base power deficit. Distinct denominators still
  cross-multiply (`1/x + 1/y → (x+y)/(x·y)`); a shared factor is no longer
  duplicated, and an exact cancellation (`x/(x+1)+1/(x+1)`) collapses through the
  canonical `Mul` to `1`.
- **Verified against SymPy:** `a/b+c/b`, the 3-term `a/b+c/b+d/b`, the `(x+1)`
  and `(x-1)²` shared-factor cases, `1/x+1/x**2`, and the distinct-denominator
  baseline all match. The fix flows through to `simplify` (the SIMP-4 follow-up
  `simplify(a/b+c/b) → (a+c)/b` now works) and leaves `cancel` / `apart`
  unchanged.
- **Regression test:** `tests/polys/poly_test.cpp`
  — `[4][together][oracle][regression]` (TOGETHER-LCM-1).
- **Scope:** structural factor sharing (identical bases, power relationships).
  Denominators sharing a *non-obvious* polynomial factor (`x²−1` and `x+1`)
  still combine over their product — that needs the multivariate-GCD work
  (CANCEL-1) to detect; the result stays correct, just not minimal.

### MAT-CHARPOLY-1 — `Matrix::charpoly` returned an unexpanded cofactor form
- **Input:** `Matrix{{1,2},{3,4}}.charpoly(λ)` and other square matrices.
- **Was:** `(λ - 1)*(λ - 4) - 6` — the raw cofactor-expansion shape produced by
  `det(λI − A)`. Mathematically a characteristic polynomial, but not the form a
  caller expects.
- **Expected (SymPy):** the expanded, like-terms-collected polynomial
  `λ² − 5λ − 2` (SymPy's `charpoly().as_expr()`).
- **Fix (`src/matrices/matrix.cpp`):** `charpoly` now returns
  `expand(det(λI − A))`. `eigenvals` is unaffected (it rebuilds a `Poly`, which
  expands regardless).
- **Verified against SymPy:** the expanded polynomial matches for 2×2, 3×3,
  singular, and symbolic matrices (`λ² − 2λx + x² − 1` for `[[x,1],[1,x]]`).
  (Term *ordering* still differs — SymPP's canonical `Add` order vs SymPy's
  descending-degree — but the polynomials are identical; ordering is a separate
  printer concern.)
- **Regression test:** `tests/matrices/matrix_test.cpp`
  — `[9][matrix][charpoly][oracle][regression]` (MAT-CHARPOLY-1): no surviving
  `)*(` factor, and oracle-equivalence for 2×2 / 3×3 / symbolic.
- **Scope:** the rest of the dense-matrix surface (det, inverse, rank, rref,
  nullspace, eigenvals, eigenvects) was cross-checked against SymPy in this pass
  and already matches.

### INT-RECIP-2 — `∫1/cosh(x)` (reciprocal hyperbolic as a `Pow`) wasn't integrated
- **Input:** `integrate(1/cosh(x))`, `integrate(1/sinh(x))`, and affine variants
  (`1/cosh(2x)`, `1/sinh(3x+1)`).
- **Was:** the unevaluated `Integral(cosh(x)**(-1), x)` — the hyperbolic analogue
  of INT-RECIP-1. `integrate(sech(x))` / `integrate(csch(x))` worked, but the
  `Pow(cosh(x), -1)` form fell through.
- **Expected (SymPy):** `∫1/cosh(x) = atan(sinh(x))`,
  `∫1/sinh(x) = log(tanh(x/2))`.
- **Fix (`src/integrals/integrate.cpp`):** extended the INT-RECIP-1
  reciprocal-first-power `Pow` branch with `Pow(cosh(u), -1) → sech(u)` and
  `Pow(sinh(u), -1) → csch(u)`, reusing the Sech/Csch antiderivatives
  (`atan(sinh)/a`, `log(tanh(u/2))/a`).
- **Verified against SymPy:** all four inputs integrate (each confirmed by
  differentiating back), and `∫1/cosh(x)` now equals `∫sech(x)` structurally.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[7][integrate][hyperbolic][oracle][regression]` (INT-RECIP-2).

### INT-RECIP-1 — `∫1/cos(x)` (reciprocal trig as a `Pow`) wasn't integrated
- **Input:** `integrate(1/cos(x))`, `integrate(1/sin(x))`, and affine variants
  (`1/cos(2x+1)`, `1/sin(3x)`).
- **Was:** the unevaluated `Integral(cos(x)**(-1), x)` — even though
  `integrate(sec(x))` and `integrate(cos(x)**(-2))` both worked. `1/cos(x)`
  parses as `Pow(cos(x), -1)`: the `exp == -1` log rule only fires for an
  *affine* base (not `cos(x)`), and the reciprocal-trig branch only matched
  `exp == -2`, so the first power fell through.
- **Expected (SymPy):** the same antiderivative as `∫sec(x)` /
  `∫csc(x)` (a half-angle log form, ≡ `log|sec x + tan x|`).
- **Fix (`src/integrals/integrate.cpp`):** a reciprocal-first-power branch in the
  `Pow` case — `Pow(cos(u), -1) → sec(u)` and `Pow(sin(u), -1) → csc(u)` for an
  affine `u` route to the exact antiderivatives the Sec/Csc function table
  already used.
- **Verified against SymPy:** all four inputs integrate (each confirmed by
  differentiating back to the integrand), and `∫1/cos(x)` now equals `∫sec(x)`
  structurally.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[7][integrate][trig][oracle][regression]` (INT-RECIP-1).
- **Scope:** `1/cos`, `1/sin` of an affine argument. The hyperbolic reciprocals
  written as a `Pow` are done in INT-RECIP-2 above.

### INT-IMPROPER-1 — improper rational functions over a linear denominator weren't integrated
- **Input:** `integrate(x/(x+1))`, `integrate(x**2/(x+1))`,
  `integrate((x**2+1)/(x-1))`, `integrate((x+1)/x)`.
- **Was:** the unevaluated `Integral(...)` marker. `try_rational` does the
  polynomial division (`x/(x+1) → 1 + (−1)/(x+1)`), but when `apart` left the
  proper remainder as a single `c/(x+a)` term, the `apart_form == proper` branch
  only handled a **degree-2** denominator (`if (den_p.degree() != 2) return
  nullopt`) and dropped everything else — so a linear denominator fell through to
  the marker even though its integral is an elementary log.
- **Expected (SymPy):** `x - log(x + 1)`, `x**2/2 - x + log(x + 1)`,
  `x**2/2 + x + 2*log(x - 1)`, `x + log(x)`.
- **Fix (`src/integrals/integrate.cpp`):** in that branch (with `q ≠ 0`) the
  proper remainder is now closed by the **general** integrator
  `integrate(proper, var)` instead of the degree-2-only helpers. That reaches the
  affine-log rule for a linear denominator and the arctan / linear-over-quadratic
  helpers for a quadratic one; the remainder is a *proper* fraction so its own
  `try_rational` bails immediately (`q' == 0`), giving no recursion.
- **Verified against SymPy:** all four inputs now integrate (each confirmed by
  differentiating the antiderivative back to the integrand), and the
  previously-working quadratic-denominator improper cases (`x**3/(x**2-1)`,
  `x**2/(x**2+1)`) are unchanged.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[7][integrate][rational][oracle][regression]` (INT-IMPROPER-1): the linear
  family verified by differentiation, plus a quadratic-denominator
  no-regression guard.
- **Scope:** rational integrands. (`1/cos(x)` written as `Pow(cos(x), -1)` is
  fixed in INT-RECIP-1 above.)

### SPECVAL-1 — `gamma` poles and `polygamma(n, 1)` special values weren't evaluated
- **Input:** `gamma(0)`, `gamma(-1)`, `gamma(-3)`; `polygamma(0, 1)` /
  `digamma(1)`, `polygamma(1, 1)`, `polygamma(2, 1)`, `polygamma(3, 1)`.
- **Was:** `gamma(0)` etc. stayed an unevaluated `gamma(0)` (the factory comment
  even read *"= pole; keep symbolic"*), and every `polygamma(n, 1)` stayed
  symbolic — including `digamma(1)`, which the DIGAMMA-1 entry had flagged as a
  follow-up (`ψ(1) = −γ`).
- **Expected (SymPy):** `gamma(non-positive integer) = zoo`;
  `polygamma(0,1) = -EulerGamma`, `polygamma(1,1) = pi**2/6`,
  `polygamma(2,1) = -2*zeta(3)`, `polygamma(3,1) = pi**4/15`.
- **Fix (`src/functions/combinatorial.cpp`):**
  - `gamma`: a non-positive integer argument now returns `S::ComplexInfinity()`
    (a simple pole) instead of an unevaluated node.
  - `polygamma`: at `x = 1` with non-negative integer order `n`,
    `ψ⁽⁰⁾(1) = −γ` and `ψ⁽ⁿ⁾(1) = (−1)^(n+1)·n!·ζ(n+1)` (the `(−1)^(n+1)` folds
    to ±1 through the existing parity rule in the `pow` factory; `ζ(even)` then
    closes to a `π` power via the existing `zeta` evaluation).
- **Verified against SymPy:** all listed inputs match; `gamma` of positive
  integers / half-integers and `polygamma` of a non-unit argument
  (`polygamma(1, 2)`, `polygamma(1, x)`) are unaffected (the rule does not
  over-fire).
- **Regression test:** `tests/functions/combinatorial_test.cpp`
  — `[3i][gamma][regression]` and `[3i][polygamma][oracle][regression]`
  (SPECVAL-1).
- **Scope:** `gamma` poles and the `x = 1` polygamma family. The
  `polygamma(n, x)` recurrence for other integer `x` (`ψ⁽¹⁾(2) = −1 + π²/6`) and
  `harmonic`/`li`/Bessel special values (those functions have no `FunctionId`
  implementation yet) are follow-ups.

### SOLVE-RAD-1 — `solve` couldn't handle radical equations (`sqrt(x) = 2`)
- **Input:** `solve(sqrt(x) - 2, x)`, `solve(x**(1/3) - 2, x)`,
  `solve(x**(2/3) - 4, x)`, `solve(sqrt(x+1) - 2, x)`, `solve(sqrt(x) - y, x)`.
- **Was:** `[]` for all of them. The polynomial path can't build a `Poly` over a
  fractional power, and the `invert_solveset` chain only peeled `Function` heads
  (log/exp/sin/…), bailing on a `Pow` head — and `solve` only routed to
  `solveset` when the expression contained a `Function` of the variable, never a
  radical.
- **Expected (SymPy):** `[4]`, `[8]`, `[8]`, `[3]`, `[y**2]`.
- **Fix (`src/solvers/solve.cpp`):**
  - `invert_solveset` gained a `Pow` branch: `g^p = c` with `p` a **non-integer**
    rational inverts to `g = c^(1/p)`, recursing on `g` when it isn't the bare
    variable. Integer powers are deliberately left to the polynomial solver so
    `x**2 = 4` still yields **both** `±2` (not just the principal root).
  - Principal-branch convention (matches SymPy): a provably-negative real RHS
    gives `∅` (`sqrt(x) = −2`, `x**(1/3) = −2`).
  - `solve` now also routes to `solveset` when the equation carries a radical of
    the variable (new `has_radical_of_var`), not only a `Function`.
- **Verified against SymPy:** all five inputs match, the negative-RHS cases give
  `[]`, the symbolic RHS gives `[y**2]`, and integer powers (`x²−4 → [−2, 2]`,
  `x³−8`, `x²−1`) are unchanged.
- **Regression test:** `tests/solvers/solve_test.cpp`
  — `[10][solve][regression]` (SOLVE-RAD-1): each radical form, the
  no-real-solution branch, and the integer-power no-regression guard.
- **Scope:** single radical head reachable through the invert chain. Equations
  mixing a radical with polynomial terms (`sqrt(x) + x − 6`) still need the
  general radical-elimination machinery and are a follow-up.

### POLE-SIGN-1 — `limit` at a finite pole returned unsigned `zoo` instead of `±oo`
- **Input:** `limit(1/x**2, x, 0)`, `limit(1/x**4, x, 0)`,
  `limit(-1/x**2, x, 0)`, `limit(1/(x-1)**2, x, 1)`.
- **Was:** `zoo` for all of them. After ZERODIV-1, direct substitution at the
  pole correctly produced `zoo` (the unsigned 1/0), but `limit` returned that as
  the answer without analysing the approach direction.
- **Expected (SymPy):** `oo`, `oo`, `-oo`, `oo` — an even-order pole diverges
  with the same sign from both sides.
- **Fix (`src/calculus/limit.cpp`):** new `signed_pole` — when direct
  substitution at a finite numeric target yields `zoo`, sample the sign of the
  expression at `target ± 1e-6` (exact substitution + `evalf`, reusing the
  infinity atoms for already-infinite samples). Matching signs ⇒ `+oo` / `-oo`;
  opposite signs ⇒ the limit is genuinely two-sided and stays `zoo`; an
  inconclusive sample (non-real, no definite sign) also stays `zoo`.
- **Verified against SymPy:** all four even-pole inputs match (`±oo`), plus
  scaled/shifted variants (`2/(x-3)**2 → oo`, `-5/x**4 → -oo`). An **odd** pole
  (`1/x`, `1/x**3`) keeps `zoo`: it is `+∞` from the right and `−∞` from the
  left, so the two-sided limit is genuinely the unsigned `zoo`. SymPy reports
  `oo` there only because its `limit` defaults to `dir='+'` (one-sided);
  SymPP's `limit` is two-sided, so `zoo` is the correct answer — a principled,
  documented divergence.
- **Regression test:** `tests/calculus/series_limit_test.cpp`
  — `[6][limit][infinity][regression]` (POLE-SIGN-1): even poles → `±oo`,
  shifted pole, odd pole stays `zoo`.
- **Scope:** finite numeric targets. Symbolic targets and essential
  singularities are out of scope; the one-sided `limit(f, x, c, dir)` API itself
  remains a separate feature gap.

### SOLVE-VAR-1 — `solve` returned a "solution" still containing the solve variable
- **Input:** `solve(x*exp(x) - 1, x)`, `solve(x*exp(x) - 2, x)`,
  `solve(exp(x) + x, x)`, `solve(x*log(x) - 1, x)`.
- **Was:** `[exp(x)**(-1)]`, `[2*exp(x)**(-1)]`, `[-exp(x)]`, `[log(x)**(-1)]` —
  every one a *rearrangement* that still contains `x`, i.e. not a solution at
  all. `solve_poly` builds a polynomial in `x` and treats a var-dependent
  "coefficient" (the `exp(x)` in `x·exp(x) − 1`) as a constant, then solves the
  apparent linear equation `x = 1/exp(x)`.
- **Expected (SymPy):** `LambertW(1)`, `LambertW(2)`, `-LambertW(1)`,
  `exp(LambertW(1))`. SymPP has no Lambert-W solver, so the correct answer is
  *"none found"* (empty) — never a `x`-containing value.
- **Fix (`src/solvers/solve.cpp`):** a correctness guard — a genuine solution
  `x = c` must be free of `x`. `solve()` now drops any candidate with
  `has(root, var)` (both from the `solve_poly` path and from the `solveset`
  finite-set path), and `solveset_impl`'s polynomial fallback rejects the same
  rearrangements before building its `FiniteSet`.
- **Verified against SymPy:** the four inputs now return `[]` (no false
  solution); every genuine solve is preserved — `x²−1 → [−1, 1]`,
  `x²−5x+6 → [2, 3]`, `log(x)−1 → [E]`, `exp(x)−2 → [log 2]`, `x³−x → [−1,0,1]`,
  `x²−y → [±√y]` (free of `x`).
- **Regression test:** `tests/solvers/solve_test.cpp`
  — `[10][solve][regression]` (SOLVE-VAR-1): the four Lambert-W inputs yield no
  `x`-containing root, plus sanity guards that genuine polynomial/parametric
  solves still return their roots.
- **Scope:** this removes the *wrong* answers. Actually solving these (Lambert-W)
  and radical equations like `sqrt(x) − 2 = 0 → 4` (still `[]`, since `sqrt` is a
  `Pow` the invert chain doesn't peel) are separate missing-feature follow-ups.

### ZERODIV-1 — `0^(negative)` (i.e. `1/0`) escaped as the malformed `0**(-1)`
- **Input:** `1/0`, `pow(0, -1)`, `0**(-2)`, `0**(-1/2)`, `2/0`, and
  `limit(1/x**2, x, 0)`.
- **Was:** the literal unevaluated `0**(-1)` — a malformed object (neither a
  number nor an infinity). The `pow` factory fell through to `number_pow(0, -1)`
  which returns `nullopt` for division by zero (`src/core/number_arith.cpp`
  already carried the comment *"0^negative — ComplexInfinity (Phase 1e+)"* but
  never produced it), so the unevaluated `Pow` leaked out. It then poisoned
  downstream consumers: `limit(1/x**2, x, 0)` returned `0**(-1)` instead of an
  infinity.
- **Expected (SymPy):** `zoo` (ComplexInfinity) for every `0**(negative)` — SymPy
  gives `0**-1 == 0**-2 == 0**Rational(-1,2) == zoo`.
- **Fix (`src/core/pow.cpp`):** an explicit early rule — `base == 0` and a
  provably-negative exponent → `S::ComplexInfinity()`. Placed after
  `pow_with_infinity` (so `0**(-oo)` is still handled there) and after the
  `x**0 → 1` rule (so `0**0 = 1` wins). `0**(positive)` (→ 0) and symbolic /
  unknown-sign exponents are untouched.
- **Verified against SymPy:** `1/0`, `0**-1`, `0**-2`, `0**(-1/2)`, `2/0`,
  `1/(x-x)` all give `zoo`; `0**2 → 0`, `0**0 → 1`, `x**-1`, `5**-1 → 1/5`
  unchanged.
- **Regression test:** `tests/core/infinity_test.cpp`
  — `[1][infinity][regression]` (ZERODIV-1): the `0**negative` family, `1/0`,
  the unaffected non-singular cases, and `limit(1/x**2, x, 0)` no longer leaking
  `0**(-1)`.
- **Scope:** this fixes the malformed-output bug. Refining the pole `zoo` to the
  signed `±oo` (so `limit(1/x**2, x, 0) = oo`) is done in POLE-SIGN-1 above.

### TRIG-PWR — `trigsimp` didn't apply the power-reduction / half-angle identities
- **Input:** `(1 − cos 2x)/2`, `(1 + cos 2x)/2`, `1 − cos 2x`, `1 + cos 2x`,
  `3·(1 − cos 2x)/2`.
- **Was:** unchanged (`1/2 − cos(2x)/2`, …). `trigsimp_add` collapsed sums of
  `a·sin²x + b·cos²x`, but a `cos(2x)` term (cosine to the first power) was not
  recognised, so the reverse power-reduction direction never fired.
- **Expected (SymPy):** `sin²x`, `cos²x`, `2·sin²x`, `2·cos²x`, `3·sin²x`.
- **Fix (`src/simplify/simplify.cpp`):** a `q·cos(2·g)` term is now folded into
  the per-argument sin²/cos² buckets via `cos(2g) = cos²g − sin²g`
  (`as_cos_double_term` / `cos_double_arg`, restricted to a literal integer-2
  factor in the argument). A third **cos-based Pythagorean** candidate
  (`a + (b − a)·cos²x`) was added alongside the existing sin-based and
  double-angle candidates; `trigsimp_add` returns whichever of the three has the
  fewest leaves. Because the selection always takes the global minimum, there is
  no oscillation: a bare `cos(2x)` stays `cos(2x)`, and the existing
  `1 − 2·sin²x → cos 2x` collapse is preserved (the `cos 2x` form has fewer
  leaves there).
- **Verified against SymPy:** all five inputs plus the no-oscillation guards
  match `sympy.trigsimp`. `(1 − cos 4x)/2` stays unreduced in **both** SymPP and
  SymPy (the literal-`cos(2·g)` restriction mirrors SymPy's own behaviour).
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[5][trigsimp][oracle][regression]` (TRIG-PWR): the (1 ∓ cos 2x)/2 family,
  scaled/unhalved variants, and the `cos(2x)` / `1 − 2·sin²x` no-oscillation
  guards.
- **Scope:** the `cos(2·g)` ↔ sin²/cos² family. `sin⁴x − cos⁴x → −cos 2x` (a
  4th-power difference) is a separate, narrower follow-up not covered here.

### SIMP-4 — `simplify` could return a form *larger than its input*
- **Input:** `simplify((x+1)**3)`, `simplify((x+1)**2)`,
  `simplify((exp(x)-1)/(exp(x/2)+1))`.
- **Was:** `3*x + 3*x**2 + x**3 + 1`, `2*x + x**2 + 1`, and a 14-node nested
  fraction (`((exp(x/2)+1)*exp(x) - (exp(x/2)+1))*(exp(x/2)+1)**(-2)`) — each
  *bigger* than the input. The pipeline expands eagerly (`expand()` at step 2)
  and never compared the final result against the original, so already-compact
  forms got inflated. (SIMP-1 had added a *local* strictly-simpler guard for the
  univariate rational-cancel step only; the global pipeline had none.)
- **Expected (SymPy):** `(x + 1)**3`, `(x + 1)**2`, `(exp(x) - 1)/(exp(x/2) + 1)`
  — SymPy's `simplify` guarantees it never returns something more complicated
  than the input (it scores candidates by a complexity measure).
- **Fix (`src/simplify/simplify.cpp`):** a global anti-bloat guard at the end of
  `simplify()` — when `node_count(current) > node_count(canon)` (the canonical
  input), return `canon`. Genuine reductions are unaffected because they shrink
  the node count: `(x+1)*(x-1) → x**2-1`, `sin²+cos² → 1`, `(x²-1)/(x-1) → x+1`
  all still fire. Rationalization that legitimately grows the count
  (`1/√2 → √2/2`, `√8 → 2√2`) is preserved (radsimp's form is not larger by
  `node_count`).
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[5][simplify][regression]` / `[5][simplify][oracle][regression]`
  ((x+1)**2/(x+1)**3 stay factored; genuine reductions still fire; exp fraction
  no longer bloats).
- **Scope:** the guard prevents *growth*; it does not add new reductions. Cases
  where SymPy reduces *below* the input but SymPP cannot yet — `exp(x/2)-1` from
  the exp fraction (needs generator-aware `cancel` in `exp(x/2)`), or `(a+c)/b`
  from `a/b + c/b` (a `together` defect that emits `b**-2` at equal node count)
  — remain separate follow-ups; SymPP now returns the input form rather than a
  bloated one in those cases.

### DERIV-1 — derivatives of undefined / untabulated functions were silently `0`
- **Input:** `diff(f(x))` for an undefined `f`, `diff(x*g(x))`, `diff(f(x)**2)`,
  `diff(besselj(0,x))`, `diff(zeta(x))`, `diff(li(x))`, `diff(beta(x,y))`,
  `diff(fresnels(x))`.
- **Was:** `0` (and, in products/sums, the function term was *dropped*:
  `diff(x*g(x)) = g(x)`). `Function::diff_arg`'s default returned `S::Zero()`,
  so **every** function without an explicit derivative rule — including all
  undefined functions `f(x)` — differentiated to a silently-wrong `0`. The
  header already documented a "Derivative marker" contract that the
  implementation never honoured. DIGAMMA-1 had patched `gamma`/`loggamma`
  case-by-case but left the unsafe default and the architectural hole (no
  `Derivative` node existed).
- **Expected (SymPy):** `Derivative(f(x), x)`, `x*Derivative(g(x), x) + g(x)`,
  `2*f(x)*Derivative(f(x), x)`, `Derivative(besselj(0, x), x)` (or the closed
  form), `Derivative(zeta(x), x)`, etc. — never `0` for a var-dependent argument.
- **Fix:**
  - New unevaluated **`Derivative`** node (`include/sympp/core/derivative.hpp`,
    `src/core/derivative.cpp`, `TypeId::Derivative` which already existed in the
    enum). Holds `(expr, var, order)`; prints `Derivative(f(x), x)` and
    `Derivative(f(x), (x, n))` for higher order, matching SymPy. The `derivative()`
    factory folds same-variable nesting into a single bumped order.
  - `Function::diff_arg` default now returns `derivative(shared_from_this(),
    args_[i])` — the unevaluated partial w.r.t. the i-th argument slot — instead
    of `0`. Subclasses with a closed form (sin, exp, gamma, Si, Ei, erf, …) still
    override and are unaffected.
  - `diff()` computes the inner derivative *before* the partial (so an argument
    independent of `var` short-circuits to `0` with no spurious
    `Derivative(f, const)`), and handles `TypeId::Derivative` for higher orders.
  - With this, the chain/product/power/sum rules compose the node correctly:
    `diff(f(x))`, `diff(x*g(x))`, `diff(f(x)**2)`, `diff(f(x)+x**2)` all match
    SymPy **exactly**; the previously-wrong specials (besselj, li, beta,
    fresnels, …) now return a correct unevaluated `Derivative` rather than `0`.
- **Regression test:** `tests/calculus/diff_test.cpp`
  — `[6a][diff][derivative][regression]` (undefined-function Derivative,
  product/power/sum carry, independent-variable → 0, second-order order bump).
- **Scope:** the node makes the result *correct* (never wrong). Adding the
  remaining closed-form reductions SymPy applies (besselj recurrence,
  `li'(x)=1/log(x)`, `fresnels'(x)=sin(πx²/2)`, beta via polygamma) is a
  mechanical follow-up — `li`/`fresnels`/`fresnelc` first need a `FunctionId`
  (they currently parse as undefined functions, for which `Derivative` already
  matches SymPy). `zeta` already matches SymPy (both keep it unevaluated).

### DIGAMMA-1 — `gamma`/`loggamma` derivatives were silently `0`
- **Input:** `diff(gamma(x))`, `diff(loggamma(x))`, `diff(gamma(x²))`.
- **Was:** `0` — `GammaFn`/`LogGamma` had no `diff_arg` override, so they fell
  through to `Function::diff_arg`'s default of `0`. Differentiating either gave a
  wrong answer with no error.
- **Expected (SymPy):** `gamma(x)·polygamma(0, x)`, `polygamma(0, x)`,
  `2·x·gamma(x²)·polygamma(0, x²)`.
- **Fix:** new `polygamma(n, x)` special function (`src/functions/combinatorial.cpp`,
  `FunctionId::PolyGamma`), kept symbolic for symbolic arguments as SymPy does,
  with `∂/∂x polygamma(n,x) = polygamma(n+1, x)`. `GammaFn::diff_arg` now returns
  `Γ(x)·polygamma(0,x)` and `LogGamma::diff_arg` returns `polygamma(0,x)`.
  `digamma(x)` is provided as sugar for `polygamma(0, x)` (SymPy's canonical form
  for `ψ`). Parser accepts `polygamma` and `digamma`.
- **Regression test:** `tests/functions/combinatorial_test.cpp`
  — `[gamma][diff][oracle][regression]`.
- **Scope:** the derivative chain (Γ, logΓ, ψ⁽ⁿ⁾). `polygamma` is left symbolic —
  numeric special values (`ψ(1) = −γ`, etc.) and `factorial`'s derivative remain
  follow-ups.

### TRIG-PYTH — `trigsimp` didn't apply the additive trig Pythagorean identities
- **Input:** `1 + tan²x`, `sec²x − tan²x`, `csc²x − cot²x`, `1 + cot²x`,
  `tan²x − sec²x`, `3 + 3tan²x`.
- **Was:** unchanged — `trigsimp` had the `sin²+cos²` collapse and (after
  TRIG-HYP-4) the hyperbolic analogue, but no `tan/cot/sec/csc` Pythagorean.
- **Expected (SymPy):** `cos⁻²x`, `1`, `1`, `sin⁻²x`, `−1`, `3·cos⁻²x`.
- **Fix (`src/simplify/simplify.cpp`):** new `trig_pyth_add` (run inside
  `trigsimp_node`) — the analogue of `tanh_coth_pyth_add` with the opposite sign
  (`sec² − tan² = 1`): rewrites each squared `tan/cot/sec/csc` term into the
  `cos⁻²`/`sin⁻²` basis via `tan² = cos⁻² − 1`, `cot² = sin⁻² − 1`,
  `sec² = cos⁻²`, `csc² = sin⁻²`, kept only when it shrinks the number of
  additive terms (so a bare `tan²x` or `2 + tan²x` is left untouched).
- **Follow-on (`src/integrals/integrate.cpp`):** because `simplify(d/dx tan x)`
  now folds `1 + tan²x → cos⁻²x`, heurisch's `u = tan x` substitution lost its
  rational-in-`g` form and fell through to the (latent-buggy) Weierstrass path,
  which hung on `∫1/(1 + tan x)`. heurisch now tries both `simplify(g')` and the
  raw `diff(g)`; the raw `1 + tan²x` keeps the substitution closed, and the
  integral once again resolves directly (no Weierstrass, no hang).
- **Regression tests:** `tests/simplify/simplify_test.cpp`
  (`[trigsimp][oracle][regression]`) and the updated Weierstrass guard test in
  `tests/integrals/integrate_test.cpp`.
- **Scope:** the additive squared-identity family; surviving-constant sums are
  left as SymPy leaves them.

### TRIG-RATIO — `trigsimp` didn't cancel trigonometric ratio products
- **Input:** `tan x·cos x`, `cot x·sin x`, `sec x·cos x`, `csc x·sin x`,
  `cot x·tan x`, `3·tan x·cos x`, `tan²x·cos²x`.
- **Was:** unchanged — the hyperbolic ratio-cancel pass (TRIG-HYP-3) had no
  trigonometric counterpart, so a `tan/cot/sec/csc` factor was never cancelled
  against the `sin`/`cos` it shared a product with.
- **Expected (SymPy):** `sin x`, `cos x`, `1`, `1`, `1`, `3·sin x`, `sin²x`.
- **Fix (`src/simplify/simplify.cpp`):** new `trig_ratio_mul` (run inside
  `trigsimp_node`, before `trigsimp_mul`) rewrites each `tan/cot/sec/csc` factor
  (to any power) as the equivalent `sin`/`cos` power(s) and lets `Mul` recombine
  same-base powers. Kept only when it lowers the leaf count, so a lone `tan x`
  (or `2·tan x`) is left untouched.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[trigsimp][oracle][regression]`.
- **Scope:** multiplicative ratio cancellation; the additive trig Pythagorean
  identities (`1 + tan² → cos⁻²`, etc.) are a separate follow-up.

### TRIG-HYP-4 — `trigsimp` didn't apply the additive tanh/coth Pythagorean identities
- **Input:** `1 − tanh²x`, `coth²x − 1`, `sech²x + tanh²x`,
  `csch²x − coth²x`, `3 − 3tanh²x`.
- **Was:** unchanged — `trigsimp` had the hyperbolic Pythagorean for `sinh`/`cosh`
  (TRIG-HYP-1) but no analogue for the `tanh`/`coth`/`sech`/`csch` squares.
- **Expected (SymPy):** `cosh⁻²x`, `sinh⁻²x`, `1`, `−1`, `3·cosh⁻²x`.
- **Fix (`src/simplify/simplify.cpp`):** new `tanh_coth_pyth_add` (run inside
  `trigsimp_node`) rewrites each squared `tanh/coth/sech/csch` term into the
  `cosh⁻²`/`sinh⁻²` basis via `tanh² = 1 − cosh⁻²`, `coth² = 1 + sinh⁻²`,
  `sech² = cosh⁻²`, `csch² = sinh⁻²`, accumulating the loose constants. The
  rewrite is kept only when it lowers the number of additive terms, so a bare
  `tanh²x` — or `2 − tanh²x`, where the constant survives — is left untouched.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[trigsimp][oracle][regression]`.
- **Scope:** the additive squared-identity family; mixed/surviving-constant sums
  are left as SymPy leaves them.

### TRIG-HYP-3 — `trigsimp` didn't cancel hyperbolic ratio products
- **Input:** `tanh x·cosh x`, `coth x·sinh x`, `sech x·cosh x`,
  `csch x·sinh x`, `coth x·tanh x`, `3·tanh x·cosh x`, `tanh²x·cosh²x`.
- **Was:** unchanged — `trigsimp` had no rule to cancel a `tanh/coth/sech/csch`
  factor against the `sinh`/`cosh` in the same product.
- **Expected (SymPy):** `sinh x`, `cosh x`, `1`, `1`, `1`, `3·sinh x`, `sinh²x`.
- **Fix (`src/simplify/simplify.cpp`):** new `hyp_ratio_mul` (run inside
  `trigsimp_node`) rewrites each `tanh/coth/sech/csch` factor (to any power) as
  the equivalent `sinh`/`cosh` power(s) and lets `Mul` recombine same-base
  powers. The rewrite is kept only when it lowers the leaf count, so a lone
  `tanh x` (or `2·tanh x`) — which would only grow — is left untouched.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[trigsimp][oracle][regression]`.
- **Scope:** multiplicative ratio cancellation; the additive tanh/coth
  Pythagorean identities (`1 − tanh² → cosh⁻²`, etc.) remain a separate gap.

### TRIG-HYP-2 — `trigsimp` didn't rewrite `cosh±sinh` as `e^±x`
- **Input:** `cosh x + sinh x`, `cosh x − sinh x`, `3cosh x + 3sinh x`,
  `2cosh x − 2sinh x`.
- **Was:** unchanged — `trigsimp` had no rule to collapse the linear
  combination `cosh ± sinh` even though it equals a single exponential.
- **Expected (SymPy):** `eˣ`, `e⁻ˣ`, `3eˣ`, `2e⁻ˣ`.
- **Fix (`src/simplify/simplify.cpp`):** new `hyp_to_exp_add` (run inside
  `trigsimp_node`, after `hypsimp_add`) collects, per argument, the linear
  coefficients of `cosh(x)` and `sinh(x)`; when they are equal it emits
  `c·eˣ`, when opposite `c·e⁻ˣ`, otherwise leaves the term untouched.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[trigsimp][oracle][regression]`.
- **Scope:** equal/opposite cosh & sinh coefficients per argument; mixed
  coefficients (e.g. `cosh + 2sinh`) are left unchanged.

### TRIG-HYP-1 — `trigsimp` didn't apply the hyperbolic Pythagorean identity
- **Input:** `cosh²x − sinh²x`, `1 + sinh²x`, `cosh²x − 1`, `3cosh²x − 3sinh²x`.
- **Was:** unchanged — `trigsimp` collapsed `sin² + cos² → 1` but had no
  hyperbolic analogue, so `cosh² − sinh²` stayed a two-term sum.
- **Expected (SymPy):** `1`, `cosh²x`, `sinh²x`, `3`.
- **Fix (`src/simplify/simplify.cpp`):** new `hypsimp_add` (run inside
  `trigsimp_node`) collects `a·sinh²(x) + b·cosh²(x)` per argument and, via
  `cosh² − sinh² = 1`, produces both the sinh form `b + (a+b)·sinh²` and the cosh
  form `−a + (a+b)·cosh²`, keeping whichever (with the rest of the sum) has the
  fewest leaves.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[trigsimp][oracle][regression]`.
- **Scope:** the hyperbolic Pythagorean and its scaled forms; the trig
  Pythagorean (`sin²+cos²`) path is unchanged.

### SIMP-3 — `simplify` didn't pull `log` of a positive base out of `exp`
- **Input:** `exp(x + log p)`, `exp(2·log p + x)`, `exp(log p + log q + x)` for
  positive `p, q`.
- **Was:** unchanged — only the whole-argument `exp(c·log p)` folded (ASSUME-6),
  not a `log` term living inside a larger sum.
- **Expected (SymPy):** `p·eˣ`, `p²·eˣ`, `p·q·eˣ`.
- **Fix (`src/simplify/simplify.cpp`):** new `exp_log_sum` pass — for
  `exp(Add(…))`, any addend that is `c·log(p)` with `p` positive is pulled out as
  the factor `p^c`, leaving `exp(rest)`.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[simplify][assumptions][regression]`.
- **Scope:** positive log bases; a sum with no positive-log addend is unchanged.

### SIMP-2 — `simplify` didn't combine exponential products
- **Input:** `simplify(eˣ·eʸ)`, `eˣ·e⁻ˣ`, `(eˣ)²·eʸ`, `e²·e³`.
- **Was:** unchanged (`exp(x)*exp(y)`) — the canonical `Mul` keeps `exp` factors
  separate (SymPP models `exp` as a `Function`, not `Pow(E, ·)`, so the same-base
  power-merge never fires).
- **Expected (SymPy):** `eˣ⁺ʸ`, `1`, `e^(2x+y)`, `e⁵`.
- **Fix (`src/simplify/simplify.cpp`):** new `combine_exp` pass (after `powsimp`):
  in a product, sum the arguments of all `exp(a)` / `(exp(a))^k` factors into a
  single `exp(Σ)`; `e⁰` folds to `1`.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[simplify][oracle][regression]` (verified against the oracle).
- **Scope:** `simplify`-level (matches SymPy's `simplify`/`powsimp`); the core
  `Mul` still keeps `exp` products separate by default.

### ASSUME-14 — `Mod(n, 1)` not simplified for integer `n`
- **Input:** `Mod(n, 1)` for an integer symbol `n`.
- **Was:** unevaluated — `mod` folded numeric arguments and `Mod(0,q)`/`Mod(x,x)`,
  but not the integer-modulo-1 identity for a symbolic integer.
- **Expected (SymPy):** `Mod(n, 1) = 0`; a non-integer argument keeps `Mod(x, 1)`
  (= the fractional part).
- **Fix (`src/functions/integers.cpp`):** `Mod(p, 1) → 0` when `is_integer(p)`.
- **Regression test:** `tests/functions/integers_test.cpp`
  — `[mod][assumptions][regression]`.

### ASSUME-13 — `floor(n + 1/2)` not simplified for integer `n`
- **Input:** `floor(n + 1/2)`, `ceiling(n + 1/2)`, `floor(2n + x)` for integer `n`.
- **Was:** unevaluated — floor/ceiling folded an integer *symbol* (`floor(n)=n`)
  and numeric/constant arguments, but not an integer-plus-remainder sum.
- **Expected (SymPy):** `floor(n + 1/2) = n`, `ceiling(n + 1/2) = n + 1`,
  `floor(2n + x) = 2n + floor(x)`.
- **Fix (`src/functions/integers.cpp`):** floor/ceiling are integer-shift
  invariant — a new `pull_integer_shift` splits an `Add` into its
  provably-integer terms and the remainder, returning `(Σ int) + floor(rest)`.
- **Regression test:** `tests/functions/integers_test.cpp`
  — `[floor][ceiling][assumptions][regression]`.
- **Scope:** sums with a provably-integer part; a purely non-integer argument
  stays under floor/ceiling.

### ASSUME-12 — parity not inferred through Mul / Add / Pow at the `ask` level
- **Was:** after ASSUME-11 added the even/odd key, `is_even(2·n)` (the `ask`
  query) was still Unknown for an integer `n` — only the structural
  `is_provably_even` helper knew it. The two disagreed.
- **Fix:** implement parity in the node `ask`s:
  - `Mul::ask` — an integer product is even iff some factor is even, odd iff every
    factor is odd (requires all factors known integer).
  - `Add::ask` — fold the terms' parities (XOR; Unknown if any term's parity is).
  - `Pow::ask` — `base^n` for a positive integer `n` keeps the base's parity
    (`oddⁿ` odd, `evenⁿ` even).
- **Expected (SymPy):** `is_even(2n)=True`, `is_odd(2n+1)=True`,
  `is_even(e+e)=True`, `is_odd(o²)=True`, `is_odd(n·m)=None`.
- **Regression test:** `tests/core/assumptions_test.cpp`
  — `[assumptions][regression]` (product/sum/power parity, unknown cases stay).
- **Scope:** `ask`-level parity now matches `is_provably_even/odd` for compound
  integer expressions.

### ASSUME-11 — no `even` / `odd` assumption (symbol-declared parity)
- **Was:** the assumption vocabulary had no parity predicate, so `Symbol("n",
  even=True)`-style declarations were impossible and the parity consumers
  (ASSUME-7/8/9/10) only fired on *structurally* even/odd exponents (`2n`,
  `2n+1`), never on a symbol simply declared even/odd.
- **Fix:** added `Even` / `Odd` to `AssumptionKey`, mask fields + builders
  (`set_even`/`set_odd`), hash, and deductive closure:
  `even ⇒ integer (⇒ rational ⇒ real)`, `odd ⇒ integer + nonzero`,
  `zero ⇒ even`, even/odd mutually exclusive, `integer ∧ ¬even ⇒ odd`
  (and `¬integer ⇒ ¬even ∧ ¬odd`). `Integer`/`Rational` literals answer parity by
  value; `is_even`/`is_odd` query wrappers added; `is_provably_even/odd` now
  consult the declared/derived `ask(Even/Odd)` first, so the existing consumers
  (`(−1)^n`, `cos(nπ)`, …) fire for declared-parity symbols too.
- **Regression test:** `tests/core/assumptions_test.cpp`
  — `[assumptions]` closure cases + `[assumptions][regression]` (declared
  even/odd predicates, integer-literal parity, `(−1)^even=1`, `(−1)^odd=−1`).
- **Scope:** parity as a first-class assumption. Parity *inference* through
  `Mul`/`Add` at the `ask` level is still deferred to the structural
  `is_provably_even/odd` helper (which already covers `2n`, `2n+1`).

### ASSUME-10 — `cot/sec/csc` at integer / half-integer multiples of π weren't evaluated
- **Input:** `cot(nπ)`, `csc(nπ)`, `sec(nπ)`, and the odd-half-integer forms, for
  integer `n`.
- **Was:** unevaluated — the reciprocal trio reduced only numeric rational
  multiples; symbolic integer / half-integer multiples fell through (the
  ASSUME-7/9 work covered only sin/cos/tan).
- **Expected (SymPy):** `cot(nπ)=zoo`, `csc(nπ)=zoo`, `sec(nπ)=(−1)^n`;
  `sec((2n+1)π/2)=zoo`, `cot((2n+1)π/2)=0`, `csc((2n+1)π/2)=(−1)^n`.
- **Fix (`src/functions/trigonometric.cpp`):** the cot/sec/csc factories now use
  the `pi_factor` + `is_integer` / `is_provably_odd(2k)` checks: integer `k`
  poles for cot/csc (`sin=0`) and gives `(−1)^k` for sec (`1/cos`); an odd
  half-integer poles for sec (`cos=0`), gives `0` for cot, `(−1)^(k−1/2)` for csc.
- **Regression test:** `tests/functions/trigonometric_test.cpp`
  — `[trig][reciprocal][assumptions][regression]`.
- **Scope:** symbolic integer / odd-half-integer multiples of π; numeric
  multiples keep their exact path.

### ASSUME-9 — `cos((2n+1)*pi/2)` / `sin((2n+1)*pi/2)` weren't evaluated
- **Input:** `cos((2n+1)·π/2)`, `sin((2n+1)·π/2)` for integer `n`.
- **Was:** unevaluated — only integer multiples of π (ASSUME-7) and numeric
  rational multiples were handled; an odd half-integer multiple fell through.
- **Expected (SymPy):** `cos = 0`, `sin = (−1)^n`.
- **Fix:** lifted the structural parity helpers from ASSUME-8 into
  `core/queries` as `is_provably_even` / `is_provably_odd` (and refactored
  `pow.cpp` to use them). In `sin`/`cos`, when the π-coefficient `k` has `2k` a
  provable odd integer (an odd half-integer), `cos(kπ)=0` and
  `sin(kπ)=(−1)^(k−1/2)` (the exponent `expand`s to `n`).
- **Regression test:** `tests/functions/trigonometric_test.cpp`
  — `[trig][assumptions][regression]` (cos=0, sin=(−1)^n; literal `π/2`, `3π/2`
  still precise).
- **Scope:** odd half-integer multiples of π with a structurally-odd numerator.

### ASSUME-8 — `(-1)**(2*n)` wasn't folded for an integer `n`
- **Input:** `(−1)^(2n)`, `(−1)^(2n+1)`, `(−1)^(4n+3)` for integer `n`.
- **Was:** unevaluated — only a *literal* integer exponent folded (via
  `number_pow`); a symbolic exponent with a determinable parity did not.
- **Expected (SymPy):** `1` (even exponent), `−1` (odd exponent); a bare
  integer of unknown parity or a non-integer coefficient stays.
- **Fix (`src/core/pow.cpp`):** structural `provably_even` / `provably_odd`
  helpers (Integer by value; a `Mul` of integers is even iff some factor is even;
  an `Add` folds term parities) drive `(−1)^k = 1`/`−1` in the `pow` factory.
  Conservative — `2·n` is even only when `n` is a known integer (else `2·n` need
  not even be an integer, e.g. `n = 1/2`).
- **Regression test:** `tests/core/assumptions_test.cpp`
  — `[assumptions][pow][regression]` (`2n`, `2n+1`, `2n+2`, `4n+3`; unknown-parity
  and non-integer coefficient stay).
- **Scope:** base `−1`, exponent with structurally-determinable parity. (A full
  `even`/`odd` assumption-key with symbol-declared parity remains a larger
  follow-up; this covers the common `2n`/`2n+1` forms.)

### ASSUME-7 — `sin(n*pi)` / `cos(n*pi)` / `tan(n*pi)` weren't evaluated for integer `n`
- **Input:** `sin(n·π)`, `cos(n·π)`, `tan(n·π)`, `sin(2n·π)` for integer `n`.
- **Was:** unevaluated — the trig factories reduced only a *numeric* rational
  multiple of π (`pi_coefficient`), so a symbolic integer coefficient fell
  through.
- **Expected (SymPy):** `sin(n·π)=0`, `tan(n·π)=0`, `cos(n·π)=(−1)^n`; a
  non-integer / generic coefficient stays unevaluated.
- **Fix (`src/functions/trigonometric.cpp`):** new `pi_factor` helper returns the
  (possibly symbolic) coefficient `k` of `k·π`; when `is_integer(k)` the
  factories return `0` (sin/tan) or `(−1)^k` (cos). Covers `2n·π` etc. since `2n`
  is integer.
- **Regression test:** `tests/functions/trigonometric_test.cpp`
  — `[trig][assumptions][regression]` (integer `n`, `2n`, `n+1`; generic
  coefficient stays).
- **Scope:** integer-valued coefficient of π; numeric multiples keep their exact
  special-value path.

### ASSUME-6 — `exp(c*log(p))` didn't fold to `p^c` for positive `p`
- **Input:** `exp(2·log p)`, `exp(log(p)/2)`, `exp(−log p)`, `exp(x·log p)` for
  positive `p`.
- **Was:** unevaluated — the `exp` factory folded only the bare `exp(log p) = p`,
  not a scaled `c·log(p)`.
- **Expected (SymPy):** `p²`, `√p`, `1/p`, `p^x` for positive `p`; a non-positive
  (generic) base stays unevaluated (branch-cut conservative, matching the
  existing `exp(log x)` gate).
- **Fix (`src/functions/exponential.cpp`):** in `exp`, an argument that is a
  product of a single `log(p)` (with `p` positive) and a constant coefficient `c`
  folds to `pow(p, c)`.
- **Regression test:** `tests/functions/exponential_test.cpp`
  — `[exp][log][assumptions][regression]` (integer/fractional/negative/symbolic
  `c`, plus a generic base that must stay).
- **Scope:** positive base; same positivity gate as `exp(log p) = p`.

### ASSUME-5 — `Abs(p*x)` didn't pull out a positive symbolic factor
- **Input:** `Abs(p·x)`, `Abs(p·q·x)`, `Abs(n·x)` for positive `p, q` / negative `n`.
- **Was:** `Abs(p·x)` stayed — the `abs` factory pulled out only a leading
  *numeric* coefficient, not a positive (or negative) *symbol*.
- **Expected (SymPy):** `p·Abs(x)`, `p·q·Abs(x)`, `−n·Abs(x)`; a fully generic
  product `Abs(x·y)` stays (the modulus of each factor is unknown).
- **Fix (`src/functions/miscellaneous.cpp`):** using `|∏fᵢ| = ∏|fᵢ|`, the `abs`
  factory now pulls out *every* factor whose modulus is known — numeric (`|c|`),
  positive (`= f`), or negative (`= −f`) — leaving the rest under a single `Abs`.
- **Regression test:** `tests/functions/miscellaneous_test.cpp`
  — `[abs][assumptions][regression]` (positive/negative factor pulls; generic
  product stays; positive factor pulled with the rest kept under one `Abs`).
- **Scope:** factors with a provable sign; genuinely-unknown factors stay inside.

### ASSUME-4 — `expand(log(x*y))` didn't split for positive symbols
- **Input:** `expand(log(p·q))`, `expand(log(p³))` for positive `p, q`.
- **Was:** unchanged (`log(p*q)`) — `expand` only distributed `Mul`/`Pow`/`Add`
  and never touched `log`. SymPy's `expand` (with `force=False`) splits a log
  whose argument is provably positive.
- **Expected (SymPy):** `log(p)+log(q)`, `3·log(p)`; unchanged when any factor is
  not provably positive (e.g. `log(p·z)` with generic `z`).
- **Fix (`src/core/expand.cpp`):** new `expand_log_arg`, gated on positivity —
  `log(∏ aᵢ) → Σ log(aᵢ)` when every factor is positive, `log(b^k) → k·log(b)`
  when `b` is positive. Without provable positivity it is left intact (e.g.
  `log((−1)(−1)) ≠ log(−1)+log(−1)`).
- **Regression test:** `tests/core/expand_test.cpp`
  — `[expand][assumptions][regression]` (positive product/power split; generic
  factor blocks the split).
- **Scope:** `expand` of `log` of a positive product/power. (`logcombine` — the
  reverse direction — remains a separate follow-up.)

### ASSUME-3 — `simplify(Abs(x)**2)` stayed `Abs(x)**2` for a real symbol
- **Input:** `simplify(|x|²)`, `simplify(|x|⁴)` for a real `x`.
- **Was:** unchanged — the power-of-power rule only handled `(bᵖ)^q`, not a power
  of an `Abs(...)`.
- **Expected (SymPy):** `x²`, `x⁴` for real `x`; unchanged for an odd exponent
  (`|x|³`) or a generic (possibly complex) `x`.
- **Fix (`src/simplify/simplify.cpp`):** extend `pow_of_pow_node` — `|y|^(2m)` with
  `y` real and `2m` a positive even integer folds to `y^(2m)` (the only case where
  `|y|² = y²`; for complex `y`, `|y|² = y·ȳ ≠ y²`).
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[simplify][assumptions][regression]` (real `|x|²`/`|x|⁴`, plus odd-exponent
  and generic-base cases that must stay), asserted structurally.
- **Scope:** even powers of `Abs` of a real argument.

### ASSUME-2 — `is_nonnegative(x**2)` was Unknown for a real symbol
- **Input:** `is_nonnegative(x²)`, `is_positive(x²+1)` for a real `x`.
- **Was:** Unknown — `Pow::ask` derived sign facts only from a *positive* base, so
  an even power of a merely-real base inferred nothing, and the Unknown
  propagated up through `Add` (so `x²+1` wasn't provably positive either).
- **Expected:** `x²≥0` (and `x⁴≥0`) for real `x`; `x²+1>0`; `x²>0` for real
  *nonzero* `x`; odd powers and generic (possibly complex) bases stay Unknown.
- **Fix (`src/core/pow.cpp`):** in `Pow::ask`, a base that is `Real` raised to a
  positive **even integer** exponent answers `Nonnegative = true` (and
  `Positive = true` when the base is also `Nonzero`; `Nonpositive = false` for a
  nonzero base). The `Add` sign rules already cascade, so `x²+1>0` falls out.
- **Why it matters:** foundational inference — every downstream consumer of the
  sign queries (simplify's assumption-gated rules, abs, limits, integration
  domain choices) now sees `x²`, `x²+c`, `x⁴`, … as nonnegative/positive for real
  symbols.
- **Regression test:** `tests/core/assumptions_test.cpp`
  — `[assumptions][pow][regression]` (even power nonneg, nonzero⇒positive,
  `x²+1>0` via Add, odd-power and generic-base stay Unknown).

### ASSUME-1 — `simplify(sqrt(x**2))` ignored symbol assumptions
- **Input:** `simplify(√(x²))` for `x` positive / real / generic.
- **Was:** `(x²)^(1/2)` in all three cases — the canonical `Pow` leaves a
  power-of-power alone (branch cuts), and `simplify` never consulted the symbol's
  assumptions to recover the safe cases.
- **Expected (SymPy):** `x` for `x ≥ 0`, `Abs(x)` for `x` real, unchanged for a
  generic (possibly complex) `x`.
- **Fix (`src/simplify/simplify.cpp`):** new `pow_of_pow_node` in the simplify
  pipeline (after `powsimp`). For `(bᵖ)^q` it consults the assumption queries:
  `b ≥ 0 ⇒ b^(p·q)` (valid for all real `p,q`), and the `√(b²) ⇒ Abs(b)` case
  for real `b`. The existing `exp(log x)→x` / `log(exp x)→x` folds were already
  assumption-gated; this extends the same idea to roots of squares.
- **Regression test:** `tests/simplify/simplify_test.cpp`
  — `[simplify][assumptions][regression]` (positive / real / generic `√(x²)`,
  plus a nonnegative-base power-of-power), asserted structurally (the oracle uses
  generic symbols and cannot verify assumption-dependent results).
- **Scope:** first consumer-side use of the assumption engine in `simplify`; the
  broader gated rules (`log(a·b)` split, `|x|→x`) remain follow-ups.

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
- **Scope:** convergent non-geometric series stay as `Sum` markers — closing
  them needs `zeta` / series recognition. The even p-series subset is now closed
  by ZETA-EVEN; `Σ 1/k! = E` and odd-`p` zeta stay deferred. The fix guarantees
  correctness (no dropped sum), not closure.

### ZETA-EVEN — convergent even p-series `Σ 1/k^(2n)` not closed
- **Input:** `Σ_{k=1}^∞ 1/k²`, `1/k⁴`, `1/k⁶`, …, `1/k¹⁴`.
- **Was:** an unevaluated `Sum` marker (SUM-3 preserved it but couldn't close).
- **Expected (SymPy):** `ζ(2n) = rₙ·π^(2n)` — `π²/6`, `π⁴/90`, `π⁶/945`,
  `π⁸/9450`, `π¹⁰/93555`, `691·π¹²/638512875`, `2·π¹⁴/18243225`.
- **Fix (`src/calculus/summation.cpp`):** a branch matching `lo=1`, `hi=∞`,
  summand `var^m` with integer `m ≤ -2` returns `zeta(-m)` (ZETA-FN). Even
  exponents close to a `π`-power; odd `p>1` close to a symbolic `ζ(s)` (matching
  SymPy's `Sum(1/k**3).doit() = zeta(3)`); the divergent harmonic `p=1` (m=-1)
  is excluded and falls through to the `Sum` marker.
- **Regression test:** `tests/calculus/series_limit_test.cpp`
  — `[summation][zeta][regression]`.
- **Scope:** integer `s ≥ 2`. Non-power summands stay deferred.

### ZETA-FN — Riemann `zeta` was not a function type
- **Input:** `zeta(2)`, `zeta(0)`, `zeta(1)`, `zeta(-1)`, `zeta(-2)`, `zeta(3)`,
  `zeta(s)`.
- **Was:** `FunctionId::Zeta` existed in the enum but had no class/factory/parser
  — the parser made a generic node.
- **Now:** a `Zeta` function type (`functions/special.{hpp,cpp}`). Exact values:
  `zeta(0)=-1/2`, `zeta(1)=zoo` (pole), even positives `zeta(2n)=rₙ·π^(2n)`
  (`zeta(2)=π²/6` … `zeta(14)`), negative integers rational
  (`zeta(-1)=-1/12`, `zeta(-3)=1/120`, …, trivial zeros `zeta(-2k)=0`). Odd
  positive (`zeta(3)`) and symbolic args stay. Parser accepts `zeta`; `str()`
  round-trips. The even-`p` summation closure now routes through this.
- **Regression test:** `tests/functions/special_test.cpp` — `[zeta]`.
- **Scope:** integer arguments fold (even ≤14, odd-negatives ≤9); the
  derivative and non-integer/complex evaluation stay deferred.

### LAMBERT-W — `LambertW` was not a function type
- **Input:** `LambertW(0)`, `LambertW(E)`, `LambertW(-1/E)`, `LambertW(oo)`,
  `LambertW(x)`, `diff(LambertW(x))`.
- **Was:** `FunctionId::LambertW` existed in the enum but had no
  class/factory/parser — the parser made a generic node.
- **Now:** a `LambertWFn` principal-branch type (`functions/special.{hpp,cpp}`)
  — the inverse of `x·eˣ`. Exact values `W(0)=0`, `W(e)=1`, `W(-1/e)=-1`
  (the branch point, matched as the canonical `-E^(-1)`), `W(oo)=oo`; other
  arguments stay symbolic. Derivative `W'(x)=W(x)/(x·(1+W(x)))`. Parser accepts
  `LambertW`/`lambertw`; `str()` round-trips.
- **Regression test:** `tests/functions/special_test.cpp` — `[lambertw]`.
- **Scope:** principal branch, the four exact values + derivative. Numeric
  (Float) evaluation, other branches `W(x,k)`, and `W(x·eˣ)=x` inverse folding
  (branch-cut sensitive) stay deferred.

### EXPINT — Si/Ci/Ei integral functions, and ∫sin(x)/x, ∫cos(x)/x, ∫eˣ/x
- **Input:** `∫sin(x)/x`, `∫cos(x)/x`, `∫eˣ/x`, `∫sin(3x)/x`, and the functions
  `Si(x)`, `Ci(x)`, `Ei(x)`.
- **Was:** these integrands returned the unevaluated `Integral` marker (`∫eˣ/x`
  was the INT-15 by-parts hang, since fixed to bail), and `Si`/`Ci`/`Ei` weren't
  function types.
- **Now:** three special-integral function types (`Si`/`Ci`/`Ei` in
  `functions/special.{hpp,cpp}`, new `FunctionId` values). `Si(0)=0`,
  `Si(±oo)=±π/2`, `Si` odd, `Ci(oo)=0`, `Ei(0)=-oo`, `Ei(oo)=oo`; derivatives
  `Si'=sin(x)/x`, `Ci'=cos(x)/x`, `Ei'=eˣ/x`. A `try_expint_integral` helper in
  `integrate.cpp` maps `∫sin(c·x)/x → Si(c·x)`, `∫cos(c·x)/x → Ci(c·x)`,
  `∫exp(c·x)/x → Ei(c·x)` (monomial argument `c·x`, constant prefactors pulled
  out). Parser accepts `Si`/`Ci`/`Ei`; `str()` round-trips.
- **Regression test:** `tests/functions/special_test.cpp` — `[expint]`;
  `tests/integrals/integrate_test.cpp` — `[integrate][expint][regression]`
  (incl. the updated INT-15 case, now closing to `Ei(x)`).
- **Scope:** monomial argument `c·x` (no constant term); `∫sin(x²)/x`-style and
  the two-argument `Ei(x,k)`/`Eₙ` generalisations stay deferred. The hyperbolic
  analogues `Shi`/`Chi` ship in EXPINT-HYP.

### EXPINT-HYP — Shi/Chi hyperbolic integral functions, and ∫sinh(x)/x, ∫cosh(x)/x
- **Input:** `∫sinh(x)/x`, `∫cosh(x)/x`, `∫sinh(2x)/x`, `Shi(x)`, `Chi(x)`.
- **Was:** unevaluated markers; `Shi`/`Chi` weren't function types.
- **Now:** the hyperbolic mirror of EXPINT — `Shi`/`Chi` types (new `FunctionId`
  values) with `Shi(0)=0`, `Shi` odd, `Shi(±oo)=±oo`, `Chi(oo)=oo`, derivatives
  `Shi'=sinh(x)/x`, `Chi'=cosh(x)/x`. `try_expint_integral` extended:
  `∫sinh(c·x)/x → Shi(c·x)`, `∫cosh(c·x)/x → Chi(c·x)`. Parser + `str()`
  round-trip.
- **Regression test:** `tests/functions/special_test.cpp` — `[expint]`;
  `tests/integrals/integrate_test.cpp` — `[integrate][expint][regression]`.
- **Scope:** monomial argument `c·x`, as for EXPINT.

### POLYLOG — `polylog` (polylogarithm) was not a function type
- **Input:** `polylog(s,0)`, `polylog(s,1)`, `polylog(2,1)`, `polylog(2,-1)`,
  `polylog(2,z)`, `diff(polylog(s,z), z)`.
- **Was:** the parser made a generic node — no evaluation.
- **Now:** a two-argument `Polylog` type (`functions/special.{hpp,cpp}`, new
  `FunctionId`). `Li_s(0)=0`, `Li_s(1)=ζ(s)` (routed through the `zeta`
  function — so `Li_2(1)=π²/6`, `Li_3(1)=zeta(3)`), `Li_2(-1)=-π²/12`; other
  arguments stay symbolic (`Li_1(z)` is *not* folded to `-log(1-z)`, matching
  SymPy). The z-derivative is `Li_{s-1}(z)/z`. Parser accepts `polylog`; `str()`
  round-trips.
- **Regression test:** `tests/functions/special_test.cpp` — `[polylog]`.
- **Scope:** the clean special values + z-derivative. `∫log(1-x)/x` is *not*
  wired (SymPy's own answer is branch-cut-sensitive); the order-derivative
  (`d/ds`), `Li_2(1/2)`, and series expansion stay deferred.

### ERFI — `erfi` (imaginary error function), and ∫exp(+x²)
- **Input:** `erfi(x)`, `erfi(0)`, `diff(erfi(x))`, `∫exp(x²)`, `∫exp(2x²)`.
- **Was:** `FunctionId::Erfi` existed in the enum but had no class — `erfi`
  parsed generically; and `∫exp(c·x²)` for **positive** c returned the marker
  (`try_gaussian` only handled negative c → erf).
- **Now:** an `Erfi` function type (`functions/special.{hpp,cpp}`): `erfi(0)=0`,
  `erfi(±oo)=±oo`, odd, derivative `2·eˣ²/√π`. `try_gaussian` extended to
  positive c: `∫exp(a·x²) dx = √π·erfi(√a·x)/(2√a)` (so `∫exp(x²)=√π·erfi(x)/2`).
  Parser accepts `erfi`; `str()` round-trips.
- **Regression test:** `tests/functions/special_test.cpp` — `[erfi]`;
  `tests/integrals/integrate_test.cpp` — `[integrate][erfi][regression]` (the
  `manualintegrate` "intractable" case moved to `exp(x³)`, since `exp(x²)` now
  closes).
- **Scope:** pure `c·x²` exponent (no linear/constant term). No MPFR `erfi`, so
  `Float` arguments stay symbolic.

### BETA — `beta` (Euler Beta) was not a function type
- **Input:** `beta(2,3)`, `beta(5,2)`, `beta(1/2,1/2)`, `beta(a,b)`.
- **Was:** `FunctionId::Beta` existed in the enum but had no class — `beta`
  parsed generically.
- **Now:** a two-argument `Beta` type (`functions/combinatorial.{hpp,cpp}`)
  defined through `gamma`: `B(a,b)=Γ(a)·Γ(b)/Γ(a+b)`. It folds to the gamma
  ratio only when all three gammas evaluate to a closed form (positive integers,
  half-integers) — `beta(2,3)=1/12`, `beta(5,2)=1/30`, `beta(1/2,1/2)=π` — and
  stays `Beta(a,b)` symbolic otherwise. Parser accepts `beta`; `str()`
  round-trips.
- **Regression test:** `tests/functions/combinatorial_test.cpp` — `[beta]`.
- **Scope:** args where `gamma` closes. (SymPy's own auto-eval is inconsistent —
  `beta(2,3)` folds but `beta(5,2)` stays — but every folded value is
  numerically equal to SymPy's, so the oracle agrees.) The derivative
  (digamma-based) stays deferred.

### EXPINT-BYPARTS — ∫ of a special-integral function (erf, Si, Ei, …)
- **Input:** `∫erf(x)`, `∫erfi(x)`, `∫erfc(x)`, `∫Si(x)`, `∫Ci(x)`, `∫Ei(x)`,
  `∫Shi(x)`, `∫Chi(x)`.
- **Was:** the unevaluated `Integral` marker — by-parts had no standalone-function
  rule beyond `log`.
- **Now:** a whitelisted by-parts branch in `try_integration_by_parts`: for a
  standalone `f(affine)` with `f ∈ {erf, erfc, erfi, Si, Ci, Ei, Shi, Chi}`,
  `∫f dx = x·f − ∫x·f'`. Each `x·f'` is elementary (`x·erf' = 2x·e^(−x²)/√π`
  integrates; `x·Si' = sin(x)`; `x·Ei' = eˣ`; …), so it closes:
  `∫erf = x·erf + e^(−x²)/√π`, `∫Si = x·Si + cos(x)`, `∫Ei = x·Ei − eˣ`, etc.
  The whitelist is essential — a function with the default 0-derivative (gamma,
  zeta, …) would otherwise yield a bogus `x·f`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][expint][regression]`.
- **Scope:** the eight special-integral functions with affine argument. The
  inverse-trig analogues (`∫asin`, `∫atan`, …) would work by the same identity
  but are not yet whitelisted.

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

### INT-29 — `integrate(asin/acos/atan/acot/asinh/acosh/atanh)` returned the marker
- **Input:** `∫asin(x)`, `∫acos(x)`, `∫atan(x)`, `∫acot(x)`, `∫asinh(x)`,
  `∫acosh(x)`, `∫atanh(x)`, and `∫x/√(a·x²+c)`.
- **Was:** the marker — by-parts was only enabled for `log` and the
  special-integral functions (erf, Si, Ci, Ei, Shi, Chi); a standalone inverse
  trig/hyperbolic function fell through. The `asin`/`acos`/`asinh`/`acosh`
  cases additionally need `∫x/√(quadratic)`, which had no rule.
- **Expected (SymPy):** `x·asin(x) + √(1−x²)`, `x·atan(x) − log(x²+1)/2`,
  `x·asinh(x) − √(x²+1)`, `x·atanh(x) + log(x²−1)/2`, etc.
- **Fix (`src/integrals/integrate.cpp`):**
  - Extend the `by_parts_fn` whitelist in `try_integration_by_parts` to
    `Asin/Acos/Atan/Acot/Asinh/Acosh/Atanh`. By parts gives
    `∫f = x·f − ∫x·f'`, where `x·f'` is a rational (atan/acot/atanh) or
    `x/√(quadratic)` (the rest).
  - New `try_x_over_sqrt_quadratic`: `∫x/√(a·x²+c) = √(a·x²+c)/a`, matching a
    lone `var` times a `(quadratic)^(−1/2)` factor with no linear term.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][invtrig][regression]` (seven inverse-function integrals plus
  `∫x/√(quadratic)`, each verified by differentiation against the oracle).
- **Scope:** the seven inverse functions whose `x·f'` the table/heurisch
  closes. `acsc`/`asec`/`acoth`/`asech`/`acsch` reduce to integrands still out
  of scope and remain unevaluated.

### INT-30 — `integrate(tanh**n)` / `integrate(coth**n)` returned the marker or an ugly form
- **Input:** `∫tanh²`, `∫tanh³`, `∫tanh⁴`, `∫coth²`, `∫coth³`, `∫coth⁴`,
  `∫tanh(2x)³`.
- **Was:** `∫coth³` (and higher odd powers) fell through to the unevaluated
  marker; `∫tanhⁿ` was caught by `try_heurisch` (`u = tanh` substitution) and
  came out as an ugly `log(tanh ± 1)` partial-fraction expansion rather than the
  clean reduction. Only the `coth²` square (INT-26, via `try_trig_reduction`)
  was handled directly; `tanh²` had no square case there either.
- **Expected (SymPy):** `x − tanh`, `−tanh²/2 + log(cosh)`,
  `x − tanh³/3 − tanh`, `x − coth`, `−coth²/2 + log(sinh)`,
  `x − coth³/3 − coth`.
- **Fix (`src/integrals/integrate.cpp`):** new `try_tanh_coth_power`, the
  hyperbolic analogue of `try_tan_power`, dispatched after it (before heurisch).
  Both functions share one reduction (tanh from `tanh² = 1 − sech²`, coth from
  `coth² = 1 + csch²`):
  - `∫tanhⁿ = ∫tanh^(n-2) − tanh^(n-1)/((n-1)·g')`
  - `∫cothⁿ = ∫coth^(n-2) − coth^(n-1)/((n-1)·g')`
  Recurses through `integrate` to the `n=1` table case (`∫tanh = log(cosh)/g'`,
  `∫coth = log(sinh)/g'`) and the `n=0` case `∫1 = x`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][hyperbolic][regression]` (tanh/coth powers 2–4 plus an affine
  argument; each asserts no `Integral` marker leaks and verifies by
  differentiation against the oracle).
- **Scope:** integer `tanhⁿ`/`cothⁿ` with an affine argument. SymPP's `simplify`
  does not always reduce the `diff − integrand` residual to a structural 0
  (tanh/coth ↔ sinh/cosh rewrites are incomplete), but the oracle's numeric
  fallback confirms equivalence.

### INT-31 — `∫1/√(quadratic)` / `∫√(quadratic)` / `∫(linear)/√(quadratic)` with a linear term returned the marker
- **Input:** `∫1/√(x²+x+1)`, `∫1/√(2x−x²)`, `∫√(x²+2x+5)`,
  `∫(2x+3)/√(x²+x+1)`, `∫(x−1)/√(x²+4x+8)`.
- **Was:** the marker — `try_sqrt_quadratic` and `try_x_over_sqrt_quadratic`
  only matched a *pure* quadratic (no linear term, INT-20). The rational
  analogues (INT-16, `try_arctan_quadratic` / `try_linear_over_quadratic`)
  already complete the square, but the square-root branches did not.
- **Expected (SymPy):** e.g. `∫1/√(x²+x+1) = asinh(√3·(2x+1)/3)`,
  `∫1/√(2x−x²) = asin(x−1)`, `∫(2x+3)/√(x²+x+1) = 2√(x²+x+1) +
  2·asinh(√3·(2x+1)/3)`.
- **Fix (`src/integrals/integrate.cpp`):**
  - `try_sqrt_quadratic`: when `b ≠ 0`, substitute `u = x + b/(2a)` (so
    `Q = a·u² + (c − b²/(4a))`, `du = dx`) and reuse the pure-quadratic branch
    on the shifted radicand, then back-substitute `x ← x + b/(2a)`. Works for
    both the `+1/2` and `−1/2` exponents.
  - `try_x_over_sqrt_quadratic`: generalised to a linear numerator `N = p·x + q`
    over a general quadratic. Using `d/dx √Q = (2a·x+b)/(2√Q)`,
    `∫N/√Q = (p/a)·√Q + (q − p·b/(2a))·∫1/√Q`, the reciprocal term handled by
    the completing-the-square branch above.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][invtrig][regression]` (five cases incl. `a < 0` and two linear
  numerators, each verified by differentiation against the oracle).
- **Scope:** rational coefficients. The `diff − integrand` residual is not
  always a structural 0 (SymPP does not pull the completed-square constant out
  from under the radical, e.g. `√(4/3·Q) = (2/√3)√Q`), but the oracle's numeric
  fallback confirms equivalence.

### INT-32 — improper rational over an irreducible quadratic, and `∫P(x)·atan/atanh/acot` returned the marker
- **Input:** `∫x²/(x²+1)`, `∫x³/(x²+1)`, `∫x·atan(x)`, `∫x²·atan(x)`,
  `∫x·atanh(x)`, `∫x·acot(x)`, `∫(x+1)·atan(x)`.
- **Was:** the marker. Two linked causes:
  1. `try_rational` did polynomial division, but if the *proper* remainder's
     `apart()` could not split further (a single irreducible quadratic, e.g.
     `−1/(x²+1)`), it dropped the **whole** result — so `∫x²/(x²+1)` failed even
     though the quotient `1` and remainder `−1/(x²+1)` are each integrable.
     (`∫x²/(1−x²)` worked only because `1−x²` factors over ℝ.)
  2. Integration by parts had no `polynomial × inverse-function` branch, so
     `∫x·atan(x)` fell through — and even with the branch, its remaining
     `∫(x²/2)/(1+x²)` is exactly the improper-rational case above.
- **Expected (SymPy):** `∫x²/(x²+1) = x − atan(x)`,
  `∫x·atan(x) = x²·atan(x)/2 − x/2 + atan(x)/2`, etc.
- **Fix (`src/integrals/integrate.cpp`):**
  - `try_rational`: when `apart` leaves the proper part intact *and* the quotient
    is non-zero (improper input), integrate the remainder directly via the
    quadratic helpers — `try_arctan_quadratic` for a constant numerator,
    `try_linear_over_quadratic` for a linear one — rather than bailing. The bare
    proper case (`q == 0`) still defers downstream, avoiding re-entry.
  - New `polynomial × f` by-parts branch (`u = f`, `dv = rest dx`), reusing the
    `is_by_parts_fn` whitelist (factored out of the standalone branch). For
    atan/acot/atanh the remaining integral is rational; for asin/acos/asinh/acosh
    it is a polynomial over `√(quadratic)`, closed when low-degree, else the
    marker guard bails.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][rational][regression]` (improper rationals over `x²+1`) and
  `[integrate][invtrig][regression]` (poly × atan/atanh/acot), verified by
  differentiation against the oracle.
- **Scope:** `∫P(x)·asin/acos/asinh/acosh` for `deg P ≥ 1` and `∫P(x)·erf/Si/…`
  still bail (the remaining `∫P/√(quad)` or `∫P·e^(−x²)` needs trig-sub /
  Gaussian-moment reductions not yet implemented).

### INT-33 — `∫1/(a + b·cos x)`, `∫1/(a + b·sin x)` and other rational-in-trig integrands returned the marker
- **Input:** `∫1/(2+cos x)`, `∫1/(3+5cos x)`, `∫1/(1+sin x)`,
  `∫1/(2+cos x+sin x)`, etc.
- **Was:** the marker — there was no general strategy for a rational function of
  `sin x` / `cos x`; only the specific table forms and the power-reduction
  helpers were tried.
- **Expected (SymPy):** e.g. `∫1/(2+cos x) = (2√3/3)·atan(√3·tan(x/2)/3)`.
- **Fix (`src/integrals/integrate.cpp`):** new `try_weierstrass`, the half-angle
  substitution `t = tan(x/2)` (`sin x = 2t/(1+t²)`, `cos x = (1−t²)/(1+t²)`,
  `dx = 2/(1+t²) dt`). It first rewrites `tan/cot/sec/csc(x)` into `sin/cos(x)`,
  substitutes the half-angle forms, and — if no `var` survives (confirming the
  integrand was rational in the trig functions of the *bare* argument) —
  integrates the resulting rational function of `t` (closed by `try_rational`,
  including the INT-32 improper/irreducible-quadratic fix) and back-substitutes
  `t = tan(x/2)`. Dispatched **last**, after by-parts: its `tan(x/2)` output is
  uglier than the dedicated trig integrators, which still win for `∫sin`, `∫tan`,
  `∫sin²`, etc.
- **Follow-up fix 1 (hang):** the substituted integrand must be *rational* in `t`
  before integrating it — otherwise a non-rational trig integrand such as
  `√(tan x)` substitutes to `√(2t/(1−t²))`, a non-elementary algebraic integral
  that sent `integrate` into an unbounded search (a true hang, worse than the
  marker). Added an `is_rational_in(integrand, t)` guard; non-rational cases now
  bail to the marker. Regression: `∫√(tan x)`, `∫√(sin x)` must terminate.
- **Follow-up fix 2 (hang):** a *trig function raised to a power* (`∫1/(1+tan²x)`,
  `∫sec²x/(1+tan²x)`) substitutes to a high-degree nested rational in `t` whose
  normalisation (`cancel`) or integration (`try_rational`'s Poly GCD, cf. the
  CANCEL-1 family) runs away — `is_rational_in` passes it through because it *is*
  structurally rational. Added a `has_trig_power_of(expr, var)` guard that
  excludes any integrand containing `sin/cos/tan/cot/sec/csc(…var…)` as the base
  of a `Pow`; trig appearing only to the first power inside a polynomial
  denominator (the classic family, and `∫1/(1+tan x)`) is unaffected.
  Regression: `∫1/(1+tan²x)`, `∫sec²x/(1+tan²x)` must terminate.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][weierstrass][regression]` (six denominators spanning the atan,
  log, and rational sub-cases, verified by differentiation against the oracle;
  plus `∫√(tan x)` / `∫√(sin x)` asserting termination to the marker).
- **Scope:** rational functions of `sin x`/`cos x`/`tan x`/`cot x`/`sec x`/`csc x`
  with the **bare** argument `x` (affine arguments like `sin 2x`, and any
  polynomial factor, bail). Many results are correct but print in a `tan(x/2)`
  form whose derivative SymPy's `simplify` cannot reduce to the integrand, so
  the regression set is the oracle-confirmable subset.

### INT-34 — `∫sin·sinh`, `∫cos·cosh`, `∫e^x·sinh`, … (trig/exp × hyperbolic) returned the marker
- **Input:** `∫sin x·sinh x`, `∫cos x·cosh x`, `∫sin x·cosh x`, `∫cos x·sinh x`,
  `∫e^x·sinh x`, `∫e^(2x)·cosh x`, `∫sin 2x·sinh 3x`.
- **Was:** the marker — by-parts on these recurses (sinh/cosh don't terminate the
  way a polynomial factor does) and no rule rewrote the hyperbolics.
- **Expected (SymPy):** e.g. `∫sin x·sinh x = (sin x·cosh x − cos x·sinh x)/2`,
  `∫e^x·sinh x = e^(2x)/4 − x/2`.
- **Fix (`src/integrals/integrate.cpp`):** new `try_hyperbolic_to_exp`, gated on a
  product containing **both** a `sinh/cosh(affine)` factor and a
  `sin/cos/exp(affine)` factor. It rewrites `sinh g = (e^g − e^−g)/2`,
  `cosh g = (e^g + e^−g)/2`, expands, and integrates term by term: each term is a
  `c·e^(·)·sin/cos(·)` (the existing exp·trig cyclic closed form) or, after a
  local exp-merge step (`e^a·e^b → e^(a+b)`, which the canonical Mul does not do),
  a pure exponential. Pure `sinh·cosh` products (no trig/exp partner) are left to
  `try_hyperbolic_power`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][hyperbolic][regression]`. The antiderivatives print in
  exponential form while the integrand is in `sinh/cosh` form, so the test
  verifies **deterministically** by evaluating `diff(F) − e` to ~0 at fixed
  rational points (SymPy's `simplify` can't bridge the forms and its numeric
  `.equals` sampling is flaky here).
- **Scope:** affine arguments, products mixing the two families. A standalone
  hyperbolic or a pure trig product is handled by the existing dedicated rules.

### INT-35 — `∫P(x)·cosⁿ(x)`, `∫P(x)·sinⁿ(x)` (polynomial × trig/hyperbolic power) returned the marker
- **Input:** `∫x·cos²x`, `∫x·sin²x`, `∫x·sin³x`, `∫x²·cos²x`, `∫x·cosh²x`,
  `∫x·cos²(2x)`.
- **Was:** the marker — the polynomial × function by-parts branch only matched a
  *bare* `sin/cos/exp/sinh/cosh(affine)` factor, not a power of one. `∫cos²x` etc.
  already integrate (INT-18/trig-reduction), so the by-parts `v = ∫dv` step had a
  closed form available but was never tried.
- **Expected (SymPy):** e.g. `∫x·cos²x = x²/4 + x·sin(2x)/4 + cos(2x)/8`.
- **Fix (`src/integrals/integrate.cpp`):**
  - Extend the by-parts target test to accept a positive-integer power of
    `sin/cos/sinh/cosh(affine)` (an `is_byparts_target` lambda), since
    `integrate` already supplies the antiderivative of those powers. `u` is still
    the polynomial rest, so the by-parts recursion terminates as `deg u` drops.
  - `expand` the by-parts remainder `v·u'` before integrating it: for `deg u ≥ 2`
    the product `(x/2 + sin 2x/4)·2x` is a `Mul`-of-`Add` that no single strategy
    matches; expanding distributes it into a sum the linearity path integrates
    term by term (fixes `∫x²·cos²x`).
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][parts][regression]` (six cases incl. an odd power, a quadratic
  polynomial, a hyperbolic power, and an affine argument), verified
  deterministically by evaluating `diff(F) − e` to ~0 at fixed rational points.
- **Scope:** polynomial × integer power of `sin/cos/sinh/cosh(affine)`.

### INT-36 — `∫g'/(1+g²)` (heurisch substitution into an irreducible quadratic) returned the marker
- **Input:** `∫cos x/(1+sin²x)`, `∫sin x/(1+cos²x)`, `∫eˣ/(1+e^(2·)x²)`
  (`∫eˣ/(1+(eˣ)²)`), `∫1/(x(1+log²x))`.
- **Was:** the marker — `try_heurisch` correctly finds the substitution
  `u = g(x)` (g = sin, cos, exp, log) and reduces the integrand to `c/(1+u²)`,
  but its inner integration was table + `try_rational` only, and neither closes a
  bare/scaled irreducible quadratic (`try_rational` defers it; cf. INT-32).
- **Expected (SymPy):** `atan(sin x)`, `−atan(cos x)`, `atan(eˣ)`, `atan(log x)`.
- **Fix (`src/integrals/integrate.cpp`):** after the table and `try_rational`
  attempts, `try_heurisch` now pulls any leading numeric factor and falls back to
  `try_arctan_quadratic` / `try_linear_over_quadratic` on the substituted
  integrand, so `∫g'/(1+g²) = atan(g)` closes.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][heurisch][regression]` (sin/cos/exp/log substitutions), verified
  by differentiation against the oracle.
- **Scope / known limitation:** the `g = exp(x)` *denominator* cases such as
  `∫1/(eˣ+e⁻ˣ)` and `∫x/(x⁴+1)` still return the marker — there the substitution
  itself fails because SymPP does not fold `e^(2x)`/`e^(−x)` to `(eˣ)²`/`(eˣ)⁻¹`
  (the `exp(a)·exp(b)` non-combination gap) nor recognise `x⁴` as `(x²)²`, so the
  substituted integrand still depends on `x`. (The `∫sec²x/(1+tan²x)` hang noted
  here earlier is fixed by INT-33 follow-up fix 2 above — it now bails cleanly.)

### INT-37 — `∫1/(a·x²+b·x+c)ⁿ` (power of an irreducible quadratic) returned the marker
- **Input:** `∫1/(x²+1)²`, `∫1/(x²+1)³`, `∫1/(x²+4)²`, `∫1/(2x²+3)²`,
  `∫1/(x²+2x+5)²`.
- **Was:** the marker — `try_arctan_quadratic` handled only `n = 1`, and `apart`
  does not split a repeated irreducible-quadratic denominator, so `try_rational`
  bailed for `n ≥ 2`.
- **Expected (SymPy):** e.g. `∫1/(x²+1)² = atan(x)/2 + x/(2(x²+1))`.
- **Fix (`src/integrals/integrate.cpp`):** new `try_quadratic_power`, the standard
  reduction `Iₙ = x/(2(n−1)c·Qⁿ⁻¹) + (2n−3)/(2(n−1)c)·Iₙ₋₁` with `Q = a·x²+c`,
  recursing through `integrate` down to `I₁ = ∫1/(a·x²+c)` (atan / log). The
  leading coefficient `a` cancels in the derivation (`x² = (Q−c)/a`), so it does
  not appear in the formula — an earlier draft that kept an `a` factor gave a
  wrong answer for `a ≠ 1`, caught by the regression test. A linear term is
  removed first by completing the square (`u = x + b/(2a)`).
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][rational][regression]` (squares and a cube, a non-unit leading
  coefficient, and a completed square), verified deterministically by evaluating
  `diff(F) − e` to ~0 at fixed rational points.
- **Scope:** constant numerator over an integer power of an irreducible
  quadratic, rational coefficients. A non-constant numerator over a
  repeated-quadratic denominator still needs `apart` repeated-factor support.

### INT-38 — rational functions with repeated factors returned the marker
- **Input:** `∫1/((x−1)²(x+1))`, `∫1/(x²(x+1))`, `∫x³/(x²+1)²`.
- **Was:** the marker. Two linked causes:
  1. `partial_fractions_squarefree` (the undetermined-coefficients engine behind
     `apart`) bailed on any repeated factor (`m ≠ 1`).
  2. `try_rational` only recognised a denominator written as a `Pow` with
     exponent exactly `−1`, so `(x²+1)^(-2)` was not seen as a denominator at all.
- **Expected (SymPy):** `∫1/((x−1)²(x+1)) = −1/(2(x−1)) − log(x−1)/4 + log(x+1)/4`,
  `∫x³/(x²+1)² = log(x²+1)/2 + 1/(2(x²+1))`.
- **Fix:**
  - `src/polys/poly.cpp`: generalise `partial_fractions_squarefree` to repeated
    factors — a factor `fᵢ` of multiplicity `mᵢ` contributes terms `Pᵢⱼ/fᵢʲ` for
    `j = 1..mᵢ`, still an `N×N` undetermined-coefficient system
    (`N = deg den`). It now returns `nullopt` when the result is a single term
    (nothing actually split, e.g. `1/(x²+1)²` is already a partial fraction) so
    the integration pipeline does not loop on an unchanged fraction.
  - `src/integrals/integrate.cpp`: `try_rational` accepts any `base^(−n)` factor
    (`n ≥ 1`) as a denominator contribution, not just exponent `−1`.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][rational][regression]` (repeated linear factors, an `x²` factor,
  an improper repeated-quadratic, plus distinct-factor / irreducible-cubic
  regressions), verified by differentiation against the oracle.
- **Scope:** a linear numerator over a *single repeated irreducible quadratic*
  (`(x+1)/(x²+1)²`) is handled by INT-39 below.

### INT-39 — `∫(p·x+q)/(a·x²+b·x+c)ⁿ` (linear numerator over a quadratic power) returned the marker
- **Input:** `∫(x+1)/(x²+1)²`, `∫(2x+3)/(x²+1)²`, `∫(3x+2)/(x²+1)³`,
  `∫(x−1)/(x²+2x+5)²`.
- **Was:** the marker — `try_quadratic_power` (INT-37) matched only a bare
  `(quadratic)^(−n)` (constant numerator), and `apart` leaves such a fraction
  intact (it is already a partial fraction), so nothing split the linear
  numerator.
- **Expected (SymPy):** the antiderivative is a rational term plus an `atan`/`log`
  term, e.g. `∫(x+1)/(x²+1)² = −1/(2(x²+1)) + x/(2(x²+1)) + atan(x)/2`.
- **Fix (`src/integrals/integrate.cpp`):** generalise `try_quadratic_power` to a
  linear numerator. It now also matches a `Mul` of a degree-≤1 numerator with a
  `(quadratic)^(−n)` factor and splits using `d/dx Q = 2a·x+b`:
  `∫(p·x+q)/Qⁿ = (p/2a)·Q^(1−n)/(1−n) + (q − p·b/(2a))·∫1/Qⁿ`, the constant
  remainder reduced by the existing bare-power recursion.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][rational][regression]` (linear numerators over `(x²+1)²`,
  `(x²+1)³`, and a completed square), verified deterministically by evaluating
  `diff(F) − e` to ~0 at fixed rational points.
- **Scope:** rational coefficients; numerator degree ≤ 1 over an integer power of
  a quadratic. With INT-37/38 this closes the proper-rational-with-quadratic
  story.

### INT-40 — `∫P(x)/√(quadratic)` (polynomial numerator) and `∫xⁿ·asin/asinh` returned the marker
- **Input:** `∫x²/√(1−x²)`, `∫x³/√(x²+1)`, `∫x²/√(x²+2x+5)`, `∫x²·asin(x)`,
  `∫x²·asinh(x)`.
- **Was:** the marker — `try_x_over_sqrt_quadratic` handled only a *linear*
  numerator (INT-31), so a degree-≥2 numerator over a root fell through. This
  also blocked the INT-32 `∫P(x)·asin/acos/asinh/acosh` by-parts for `deg P ≥ 1`,
  whose remaining integral is exactly `∫(polynomial)/√(quadratic)`.
- **Expected (SymPy):** e.g. `∫x²/√(1−x²) = −x√(1−x²)/2 + asin(x)/2`,
  `∫x²·asin(x) = x³·asin(x)/3 + x²√(1−x²)/9 + 2√(1−x²)/9`.
- **Fix (`src/integrals/integrate.cpp`):** new `try_poly_over_sqrt_quadratic`,
  the reduction `∫xᵏ/√Q = [xᵏ⁻¹√Q − (k−1)c·∫xᵏ⁻²/√Q]/(k·a)` (pure quadratic
  `Q = a·x²+c`), recursing through `integrate` to the `k=1` (√Q/a) and `k=0`
  (asin/asinh/log) base cases. A linear term is removed first by completing the
  square; a multi-term numerator is distributed so linearity handles it.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][invtrig][regression]` (monomials over `√(1−x²)`, `√(x²+1)`, a
  completed square, plus `∫x²·asin`, `∫x²·asinh`), verified deterministically by
  evaluating `diff(F) − e` to ~0 at fixed rational points.
- **Scope:** rational coefficients; polynomial numerator over `√(quadratic)`.

### INT-41 — `∫log(polynomial)` returned the marker
- **Input:** `∫log(x²+1)`, `∫log(x²−1)`, `∫log(x²+x+1)`.
- **Was:** the marker — the standalone-log by-parts had a closed form only for
  `log(affine)`; a non-affine argument fell through.
- **Expected (SymPy):** `∫log(x²+1) = x·log(x²+1) − 2x + 2·atan(x)`,
  `∫log(x²−1) = x·log(x²−1) − 2x + log(x+1) − log(x−1)`.
- **Fix (`src/integrals/integrate.cpp`):** add a general `log(g)` standalone
  branch — by parts with `u = log(g)`, `dv = dx`, `v = x`:
  `∫log(g) = x·log(g) − ∫x·g'/g`. The remaining integrand `x·g'/g` is rational
  when `g` is a polynomial, so `try_rational` (with the INT-32 improper /
  irreducible-quadratic handling) closes it; the marker guard bails otherwise.
- **Regression test:** `tests/integrals/integrate_test.cpp`
  — `[integrate][parts][regression]` (`log` of an irreducible quadratic, a
  reducible quadratic, and an irreducible `x²+x+1`, plus `log(affine)`
  regressions), verified by differentiation against the oracle.
- **Scope:** the remaining `∫x·g'/g` must close (polynomial `g`); `log` of a
  transcendental argument bails to the marker.

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

### MUL-RAD — radical base collection left an un-collapsed numeric factor
- **Input:** `√2·√8`, `√3·√12`, `√8·√8`.
- **Was:** `2·2`, `2·3`, `2·4` — a `Mul` of bare numbers instead of `4`, `6`,
  `8`. Mul base collection summed the `½` exponents and called `pow(2, 1) = 2`,
  but that numeric result (and the numeric part of e.g. `2^(3/2) → 2·√2`) was
  pushed alongside the already-finalised numeric coefficient without merging.
  (Same root cause as the old `cot(π/4) → 1/2·2` cosmetic wart.)
- **Fix (`src/core/mul.cpp`):** a Step-4b sweep after base collection folds any
  numeric factor — and the numeric part of any `Mul` factor — back into the
  running product via `number_mul`. Cross-base radicals are deliberately **not**
  merged (`√2·√3` stays a two-factor `Mul`); that is a separate feature.
- **Regression test:** `tests/core/arithmetic_test.cpp` — `[mul][regression]`
  (collapse cases + a `2^(3/2)` extract-but-keep-radical + a `√2·√3`
  no-merge guard).
- **Scope:** numeric collapse only. `√2·√3 → √6`-style cross-base radical
  merging stays deferred.

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
- **Scope:** positive integer **or rational** base, `p ≥ 2`. A rational base
  `a/b` reduces via `(a/b)^(p/q) = (a·b^(q-1))^(p/q)/b^p` (the `b^p` division
  rationalises the denominator: `(2/3)^(2/3) = 2^(2/3)·3^(1/3)/3`,
  `(1/2)^(3/2) = √2/4`). Negative bases and negative exponents (`16^(-2/3)`)
  stay deferred. Distinct prime radicals are left unmerged (`√2·√3`, not `√6`).

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

### FIB-CAT — `fibonacci` / `catalan` were not function types
- **Input:** `fibonacci(10)`, `catalan(5)`, `fibonacci(x)`.
- **Was:** the parser turned `fibonacci`/`catalan` into generic
  undefined-function nodes — no evaluation at all.
- **Now:** two distinct function types (`Fibonacci`/`Catalan` in
  `functions/combinatorial.{hpp,cpp}`, new `FunctionId` values in the
  combinatorial 700-block). Evaluate on non-negative integers —
  `fibonacci` via GMP `mpz_fib_ui` (`fibonacci(10)=55`), `catalan` via
  `binomial(2n,n)/(n+1)` (`catalan(5)=42`) — and stay symbolic for symbolic
  args. Parser accepts both; `str()` round-trips. Safety bounds: `fibonacci`
  n ≤ 1e6, `catalan` n ≤ 1e5.
- **Regression test:** `tests/functions/combinatorial_test.cpp`
  — `[fibonacci]`, `[catalan]` (values, parse round-trip, subs, symbolic
  guards).
- **Scope:** non-negative integer arguments. Negafibonacci (`fibonacci(-n)`),
  negative `catalan`, and derivatives stay deferred (the latter matches
  `factorial`'s existing 0-derivative convention). `gcd`/`lcm` as functions
  remain a separate item.

### RFF-SUBF — `RisingFactorial`, `FallingFactorial`, `subfactorial` unimplemented
- **Input:** `RisingFactorial(3,2)`, `FallingFactorial(5,2)`, `RisingFactorial(x,2)`,
  `subfactorial(4)`.
- **Was:** their `FunctionId` values existed in the combinatorial 700-block but
  had no class/factory/parser, so the parser made generic nodes.
- **Now:** three function types in `functions/combinatorial.{hpp,cpp}`:
  - `rising_factorial(x,n)` (Pochhammer) = `x·(x+1)·…·(x+n-1)` and
    `falling_factorial(x,n)` = `x·(x-1)·…·(x-n+1)` — for a non-negative integer
    `n` they expand to the product even for symbolic `x` (`rf(x,2)=x·(x+1)`,
    matching SymPy); `n=0 → 1`; symbolic `n` stays.
  - `subfactorial(n)` = derangement count via the recurrence
    `!0=1, !1=0, !k=(k-1)(!(k-1)+!(k-2))` (`!4=9`, `!5=44`).
  - Parser accepts `RisingFactorial`/`FallingFactorial`/`subfactorial`; `str()`
    round-trips. Safety bounds (n ≤ 1e5).
- **Regression test:** `tests/functions/combinatorial_test.cpp`
  — `[rising]`, `[falling]`, `[subfactorial]`.
- **Scope:** non-negative integer order/argument. `binomial`-style negative or
  rational extensions stay deferred.

### GCD-LCM — `gcd` / `lcm` were not function types
- **Input:** `gcd(12,18)`, `lcm(4,6)`, `gcd(-12,8)`, `gcd(0,5)`, `gcd(x,y)`.
- **Was:** the parser made generic undefined-function nodes — no evaluation.
- **Now:** two distinct two-argument function types (`Gcd`/`Lcm` in
  `functions/combinatorial.{hpp,cpp}`, new `FunctionId` values). Evaluate on
  integer pairs via GMP `mpz_gcd`/`mpz_lcm` — non-negative results, with the
  edge cases (`gcd(0,0)=0`, `gcd(±n,0)=|n|`, `lcm(0,n)=0`) matching SymPy — and
  stay symbolic otherwise. Registered in the parser's two-arg table; `str()`
  round-trips. Distinct from the polynomial `gcd(Poly,Poly)` (different
  overload).
- **Regression test:** `tests/functions/combinatorial_test.cpp`
  — `[gcd]`, `[lcm]` (values incl. sign/zero edge cases, parse round-trip, subs,
  symbolic guards).
- **Scope:** integer arguments. Rational `gcd` (`gcd(1/2,1/3)=1/6`) and
  polynomial/symbolic gcd stay deferred (the latter is CANCEL-1 territory).

### AINV-RECIP — `acot`, `asec`, `acsc` were not function types
- **Input:** `acot(1)`, `asec(2)`, `acsc(2)`, `acot(0)`, `asec(0)`, `acot(x)`,
  `diff(acot(x))`.
- **Was:** the parser made generic undefined-function nodes — no evaluation,
  no derivatives.
- **Now:** three distinct inverse-reciprocal function types (`Acot`/`Asec`/
  `Acsc` in `functions/trigonometric.{hpp,cpp}`, new `FunctionId` values).
  Each folds exact values through the reciprocal-argument identity —
  `acot(x)=atan(1/x)` (`acot(0)=π/2`, odd), `asec(x)=acos(1/x)` (`asec(0)=zoo`),
  `acsc(x)=asin(1/x)` (`acsc(0)=zoo`, odd) — keeping its own node when the inner
  inverse stays unevaluated (so `acot(2)`, `asec(√3)` print symbolically, as in
  SymPy). Derivatives: `acot'=-1/(1+x²)`, `asec'=1/(x²√(1-1/x²))`,
  `acsc'=-1/(x²√(1-1/x²))`. Numeric args evalf through the inner inverse. Parser
  + LaTeX (`\operatorname{acot}`/…) updated; `str()` round-trips.
- **Regression test:** `tests/functions/inverse_trig_test.cpp`
  — `[acot]`, `[asec]`, `[acsc]` (canonical values, poles, parity, derivatives,
  parse round-trip).
- **Scope:** the inverse-reciprocal *hyperbolic* analogues
  (`acoth`/`asech`/`acsch`) ship in HYP-AINV-RECIP; `asec`/`acsc` real-domain
  assumptions (|x|≥1) stay agnostic.

### HYP-AINV-RECIP — `acoth`, `asech`, `acsch` were not function types
- **Input:** `acoth(oo)`, `asech(1)`, `acsch(0)`, `acoth(-x)`, `diff(acoth(x))`.
- **Was:** the parser made generic undefined-function nodes — no evaluation.
- **Now:** three distinct inverse-reciprocal hyperbolic types (`Acoth`/`Asech`/
  `Acsch` in `functions/hyperbolic.{hpp,cpp}`, new `FunctionId` values) — the
  hyperbolic mirror of AINV-RECIP. Fold via the reciprocal-arg identity:
  `acoth(x)=atanh(1/x)` (odd), `asech(x)=acosh(1/x)` (`asech(0)=oo`),
  `acsch(x)=asinh(1/x)` (`acsch(0)=zoo`, odd), keeping the node when the inner
  inverse stays unevaluated. Clean folded values: `acoth(±oo)=0`, `asech(1)=0`,
  `acsch(±oo)=0`. Derivatives: `acoth'=1/(1-x²)`, `asech'=-1/(x√(1-x²))`,
  `acsch'=-1/(x²√(1+1/x²))`. Parser + LaTeX updated; `str()` round-trips.
- **Regression test:** `tests/functions/hyperbolic_test.cpp` — `[reciprocal]`
  (values/poles, odd parity, derivatives, parse round-trip).
- **Scope:** SymPP's minimal inverse-hyperbolics mean complex/log special values
  (`acoth(0)=iπ/2`, `asech(2)=iπ/3`, `acsch(1)=log(1+√2)`) stay unevaluated —
  correct (just less reduced than SymPy). This completes the full trig +
  hyperbolic function family (forward, reciprocal, and both inverse sets).

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

### FRAC-1 — `frac` (fractional part) had an enum value but no implementation
- **Input:** `frac(7/2)`, `frac(-7/2)`, `frac(5)`, `frac(pi)`, `frac(x)`.
- **Was:** `FunctionId::Frac` existed but had no class/factory/parser entry, so
  the parser produced a generic undefined-function node.
- **Now:** a `Frac` function type (`functions/integers.{hpp,cpp}`) for the
  fractional part `frac(x)=x−floor(x)`, always in `[0,1)`: `frac(7/2)=1/2`,
  `frac(-7/2)=1/2` (not `−1/2`), `frac(int)=0`, `frac(pi)=pi−3`. Reuses `floor`'s
  numeric/constant folding — when `floor` evaluates, returns `x−floor(x)`, else
  keeps `Frac`. Parser accepts `frac`; `str()` round-trips. `frac(real)` is real
  and nonnegative.
- **Regression test:** `tests/functions/integers_test.cpp` — `[frac]`.
- **Scope:** numeric/constant args fold; the derivative is left unevaluated
  (matching SymPy, which returns `Derivative(frac(x), x)`).

### MOD-1 — `Mod` was not a function type
- **Input:** `Mod(7,3)`, `Mod(-7,3)`, `Mod(7,-3)`, `Mod(1/2,1/3)`, `Mod(x,x)`,
  `Mod(x,0)`.
- **Was:** the parser made a generic undefined-function node — no evaluation.
- **Now:** a two-argument `Mod` function type (`FunctionId::Mod` in the integer-
  rounding 500-block) implementing **floored** modulo `Mod(p,q)=p−q·⌊p/q⌋`, so
  the result takes the divisor's sign: `Mod(-7,3)=2`, `Mod(7,-3)=-2`,
  `Mod(-7,-3)=-1` (matching SymPy, not C's truncated `%`). Integer and rational
  pairs evaluate (via `mpz_fdiv_q` on `p/q`); the identities `Mod(0,q)=0`,
  `Mod(x,x)=0` fold symbolically; `Mod(x,0)` is kept unevaluated (SymPy raises —
  we avoid throwing). Parser accepts `Mod` and `mod`; `str()` emits `Mod(p, q)`
  and round-trips.
- **Regression test:** `tests/functions/integers_test.cpp` — `[mod]` (sign
  cases, rationals, identities, zero-divisor guard, parse round-trip).
- **Scope:** integer/rational arguments. `Mod` of floats and richer symbolic
  reductions (e.g. `Mod(2x, x)`) stay deferred.

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
