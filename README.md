# SymPP

[![CI](https://github.com/leonardoaraujosantos/SymPP/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/leonardoaraujosantos/SymPP/actions/workflows/ci.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-00599C?logo=cplusplus&logoColor=white)](https://en.cppreference.com/w/cpp/20)
[![CMake](https://img.shields.io/badge/CMake-3.25%2B-064F8C?logo=cmake&logoColor=white)](https://cmake.org/)
[![Tests](https://img.shields.io/badge/tests-1554%20passing-brightgreen)](#)
[![Oracle](https://img.shields.io/badge/oracle-SymPy%201.13%2B-3B5526?logo=python&logoColor=white)](https://www.sympy.org/)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)](#)
[![Last commit](https://img.shields.io/github/last-commit/leonardoaraujosantos/SymPP)](https://github.com/leonardoaraujosantos/SymPP/commits/main)

Modern C++20 symbolic mathematics library. Clean-room port of SymPy
algorithms with SymPy itself wired in as the validation oracle.

## Status

```
1554 tests / 5488 assertions  all passing
672 cases (2724 assertions) oracle-validated against SymPy
14 of 15 phases shipped
```

On textbook-shaped inputs (calculus, algebra, transforms, solvers, sets)
SymPP and SymPy are now effectively interchangeable; the remaining gaps are
research-grade algorithms (Gruntz limits, full Risch integration, Wang
multivariate factorization, transcendental `solve`) and whole SymPy modules
not yet ported. See [docs/04-roadmap.md](docs/04-roadmap.md#how-far-are-we-from-sympy)
for the multi-metric parity breakdown.

| Phase | Title | Status |
|---|---|---|
| 0  | Foundation & oracle harness            | ✅ |
| 1  | Core expression tree                   | ✅ |
| 2  | Assumptions                            | 🟡 full sign/number ontology with closure + Add/Mul/Pow propagation, consumed by simplify/integrate/refine; **at SymPy parity on the common predicates** (see [Assumptions](#assumptions) below). Only prime/irrational/algebraic/commutative predicates + a SAT-based `ask` solver remain deferred |
| 3  | Elementary & special functions         | ✅ |
| 4  | Polynomials                            | ✅ |
| 5  | Simplification                         | ✅ |
| 6  | Calculus                               | ✅ |
| 7  | Integration                            | ✅ |
| 8  | Transforms                             | ✅ |
| 9  | Linear algebra                         | ✅ |
| 10 | Equation solvers                       | ✅ |
| 11 | ODE / PDE                              | ✅ |
| 12 | Units                                  | ✅ |
| 13 | Code generation                        | 🟡 printers + function emission |
| 15 | Parser & MATLAB facade                 | ✅ |
| 16 | Hardening & v1.0                       | 🟡 find_package + install + CI shipped; v1.0 tag pending |

> Phase 14 (Plotting bridge) was dropped from scope. Plot via the
> code-gen pipeline or `evalf` / `vpaintegral` numeric output piped
> to the consumer's renderer of choice.

See [docs/04-roadmap.md](docs/04-roadmap.md) for the per-phase
shipped/deferred breakdown and the path to SymPy parity.

## Quick example

```cpp
#include <sympp/sympp.hpp>
#include <iostream>

int main() {
    using namespace sympp;

    auto x = symbol("x");
    auto y = symbol("y");

    // Calculus
    auto f = pow(x, integer(3)) + integer(2) * x * y;
    std::cout << "df/dx = " << diff(f, x)->str() << "\n";
    std::cout << "∫f dx = " << integrate(f, x)->str() << "\n";

    // Solve
    auto roots = solve(pow(x, integer(2)) - integer(5) * x + integer(6), x);
    for (auto& r : roots) std::cout << "root: " << r->str() << "\n";

    // Linear algebra
    Matrix A = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    std::cout << "det(A) = " << A.det()->str() << "\n";

    // Code gen
    std::cout << printing::ccode(diff(sin(x) * exp(x), x)) << "\n";
}
```

Output:
```
df/dx = 2*y + 3*x**2
∫f dx = 1/4*x**4 + x**2*y
root: 3
root: 2
det(A) = -2
exp(x)*sin(x) + exp(x)*cos(x)
```

More worked examples: [docs/08-tutorial.md](docs/08-tutorial.md).

## More examples

Every call below returns an `Expr`; the comment shows its `->str()`. All
outputs are produced by the current build and cross-checked against SymPy.

**Calculus, algebra, solving**

```cpp
using namespace sympp;
auto x = symbol("x");

diff(pow(x, integer(3)), x);                          // → 3*x**2
integrate(log(x), x);                                 // → -x + x*log(x)
integrate(pow(integer(1) + x * x, integer(-1)), x);   // → atan(x)
limit(pow(integer(1) + integer(1) / x, x),
      x, S::Infinity());                              // → E
factor(pow(x, integer(4)) - integer(1), x);           // → (x + 1)*(x - 1)*(x**2 + 1)
solve(pow(x, integer(2)) - integer(2), x);            // → [2**(1/2), -2**(1/2)]
```

**Sums & simplification** — including closed forms recently brought to parity:

```cpp
auto n = symbol("n");
auto oo = S::Infinity();

summation(pow(n, integer(-2)), n, integer(1), oo);    // → 1/6*pi**2  (Basel)
summation(pow(n * n + integer(1), integer(-1)),
          n, integer(1), oo);                         // → 1/2*(pi*coth(pi) - 1)

simplify(pow(sin(x), integer(2))
         + pow(cos(x), integer(2)));                  // → 1
simplify(sin(integer(2) * x)
         * pow(sin(x), integer(-1)));                 // → 2*cos(x)
```

**Transforms & ODEs**

```cpp
auto s = symbol("s"), t = symbol("t");

laplace_transform(sin(t), t, s);                      // → (s**2 + 1)**(-1)
inverse_laplace_transform(
    pow(pow(s + integer(1), integer(2)) + integer(1),
        integer(-1)), s, t);                          // → exp(-t)*sin(t)
inverse_laplace_transform(
    pow(pow(s, integer(2)) + integer(1), integer(-2)),
    s, t);                                            // → 1/2*(-t*cos(t) + sin(t))

FunctionSymbol y("y");
auto ode = diff(y(x), x, 2) + y(x);                   // y'' + y = 0
dsolve(ode, y(x), x);                                 // → __C0*cos(x) + __C1*sin(x)
```

**Assumptions & linear algebra**

```cpp
auto p = symbol("p", AssumptionMask{}.set_positive(true));
is_positive(pow(p, integer(2)) + integer(1));         // → std::optional<bool>{true}

Matrix M = {{integer(2), integer(0)}, {integer(0), integer(3)}};
M.eigenvals();                                        // → [3, 2]
```

## Assumptions

Symbols carry an **assumptions ontology** — the same idea as SymPy's
`Symbol('x', positive=True)`. You attach predicates at construction, the system
closes them under their logical consequences, propagates them through
`Add`/`Mul`/`Pow`, and feeds the result back into `simplify`, `integrate`, and
`refine`. Queries are **three-valued** (`std::optional<bool>`: `true` / `false` /
`unknown` (`std::nullopt`)) — an honest "I can't prove it" is never confused with
"false".

```cpp
using namespace sympp;

// Declare predicates with the AssumptionMask builder.
auto p = symbol("p", AssumptionMask{}.set_positive(true));   // p > 0
auto n = symbol("n", AssumptionMask{}.set_integer(true));    // n ∈ ℤ
auto x = symbol("x");                                        // generic (complex)

// Closure: positive ⇒ real, finite, nonzero; integer ⇒ rational ⇒ real.
is_real(p);                          // → true
is_nonzero(p);                       // → true

// Propagation through Add / Mul / Pow (incl. sign-by-parity).
is_positive(pow(p, integer(2)) + integer(1));   // → true   (p² ≥ 0, +1 > 0)
is_negative(mul(S::NegativeOne(), p));          // → true   (−p < 0)
is_integer(integer(2) * n + integer(1));        // → true   (2n+1 ∈ ℤ)
is_positive(mul(p, n));                          // → unknown (n's sign is open)

// Abs is always real and nonnegative — even for a generic complex x.
is_nonnegative(abs(x));              // → true   (SymPy returns None here)

// simplify / integrate consume the assumptions.
simplify(sqrt(pow(p, integer(2))));  // → p          (nonneg base)
simplify(sqrt(pow(x, integer(2))));  // → sqrt(x**2)  (generic — Abs, not x)
```

**Scoped `assuming` context.** Beyond construction-time facts, a thread-local
`assuming(...)` RAII scope asserts a relational fact about an *existing* symbol for
a region of code — SymPy's `with assuming(Q.positive(x)): …`. The central `ask`
consults the active scopes, so the fact propagates through `Add`/`Mul`/`Pow` and
unlocks `refine`; it retracts at scope exit (and adds no overhead when no scope is
active).

```cpp
auto x = symbol("x");                                   // no built-in assumptions
{
    assuming pos(x, AssumptionMask{}.set_positive(true));
    is_positive(mul(x, x));        // → true  (propagates through Mul)
    refine(sqrt(pow(x, integer(2))));  // → x   (√(x²) = |x| = x for x > 0)
    refine(abs(x));                    // → x
}                                       // fact retracted here
```

**Declarable predicates** (via `AssumptionMask::set_*`): `real`, `positive`,
`negative`, `zero`, `nonzero`, `nonnegative`, `nonpositive`, `finite`, `integer`,
`rational`, `even`, `odd`, `complex`, `imaginary`. They are closed under the
obvious implications (e.g. `zero ⇒ integer ∧ nonnegative ∧ nonpositive`,
`real ⇒ finite`, `nonnegative ⇒ real ∧ ¬negative`) before propagation.

**Parity with SymPy.** On a battery of sign/number queries the subsystem matches
SymPy's `ask` on every common predicate; the only divergences are cases where
SymPP is *strictly more decisive* — e.g. `Abs(x)` is reported `real`/`nonnegative`
for a generic `x`, where SymPy returns `None`. Still **deferred** (and tracked in
the roadmap): the `prime` / `irrational` / `algebraic` / `commutative` predicates
and a general SAT-based `ask` solver for arbitrary boolean combinations.

## What's in the box

- **Symbolic algebra** — Add/Mul/Pow with full canonical form,
  hash-cons cache, structural equality.
- **Number tower** — `Integer` / `Rational` (GMP) / `Float` (MPFR
  arbitrary-precision) / `ImaginaryUnit` / `Pi` / `E` / `oo` / `-oo` /
  `zoo` / `nan` (full infinity arithmetic).
- **Assumptions** — declarable predicates (`real`, sign, `nonzero`,
  `nonnegative`/`nonpositive`, `finite`, `integer`/`rational`, parity,
  `complex`/`imaginary`) with logical closure and Add/Mul/Pow propagation;
  three-valued `is_positive` / `is_real` / … queries consumed by
  `simplify`/`integrate`/`refine`. At SymPy parity on the common predicates —
  see [Assumptions](#assumptions).
- **30+ named functions** — sin/cos/tan/exp/log/sqrt/abs/floor/
  factorial/gamma/erf/heaviside/dirac_delta plus `lowergamma`/`uppergamma`
  (incomplete gamma, with positive-integer and half-integer erf/erfc closed
  forms), `Hyper` and `MeijerG` (proper Function classes with auto-eval) and the
  rest of the elementary + special + combinatorial canon.
- **Calculus** — `diff` (closed-form parameter derivatives for the gamma family —
  `∂ₐΒ(a,b)=Β(a,b)(ψ(a)−ψ(a+b))`, `H′(x)=ψ⁽¹⁾(x+1)` — with an unevaluated
  `Derivative` node for the genuinely non-elementary directions (`∂ₛγ(s,x)`,
  `∂ₙψ⁽ⁿ⁾`) and for undefined/untabulated functions — never a silent `0`),
  `integrate`
  (table + trig + reciprocal-trig/hyperbolic + parts + arctan/asin/asinh +
  rational incl. improper/linear-denominator + heurisch + Weierstrass
  half-angle + incomplete-gamma forms `∫xˢ⁻¹e⁻ˣ = γ(s,x)`, `∫₀^∞ xˢ⁻¹e⁻ᶜˣ = Γ(s)/cˢ`;
  definite integrals evaluate boundaries as limits, so `∫₀^∞ xⁿe^(-x) = n!`),
  `series`, **`limit` with a Gruntz-style leading-term engine** (L'Hôpital +
  composable asymptotic stages: power-as-exp, dominant-term `∞−∞`, power
  continuity, log-exp reduction for nested transcendentals, numerically-verified
  small-angle/Maclaurin-head for `0·∞` analytic forms, harmonic-number and
  log-Stirling (`log(n!)/(n·log n)→1`) asymptotics, a gamma-ratio asymptotic with
  Legendre/Gauss multiplication, exponential rates, unbalanced gamma powers and
  combined-denominator flattening, the full Stirling prefactor for cancelling
  products (`n!/(nⁿe⁻ⁿ)→∞`), constant-base-exponential rationals
  (`(2ˣ+3ˣ)/3ˣ→1`), hyperbolic→exp rewriting (`(sinh x+cosh x)/eˣ→1`),
  inverse-hyperbolic two-term asymptotics (`asinh x/log x→1`), an MRV rewrite for
  exponential differences — including a product carrying one, `M·(eᵖ−eᵍ)` with
  `p−q→0` (Gruntz's flagship `eˣ·(e^{1/x−e⁻ˣ}−e^{1/x})→−1`) — special-function
  asymptotic series for `erfc`/`erf` (Gaussian tail), `Ei`, Riemann `ζ`, digamma /
  higher polygamma, log-Stirling `loggamma`, the Tricomi–Erdélyi gamma-ratio
  subleading term, and the arctangent (`x·(atan x−π/2)→−1`, with subleading),
  dominant-summand extraction from a `log` of an exponential- or logarithmic-rate
  sum (`log(x+eˣ)/log(x+e²ˣ)→1/2`, `log(log x+log log x)−log log x→0`), the
  log-exp reduction for tower products (`xˣ/(x+1)^{x+1}→0`), and a
  leading-term-by-series stage for `1^∞` corrections at `±∞` and finite points
  (`x·((1+1/x)ˣ−e)→−e/2`, `(xˣ−1)/(x·log x)→1`), generalized to a series core
  times a non-polynomial multiplier (`eˣ·((1+1/x)ˣ−e)→−∞`) — resolves `Γ(2n)/nⁿ`,
  `eˣ·(sin(1/x+e⁻ˣ)−sin(1/x))→1`; size-bounded so it never hangs),
  `summation` (p-series → ζ, exponential/geometric series, telescoping,
  geometric/arithmetic-geometric, irreducible-quadratic
  `Σ1/(n²+a²) → (π√a²·coth(π√a²)−1)/(2a²)`), Padé, Euler-Lagrange, asymptotes.
- **Polynomials** — div/gcd/sqf, factor over ℤ (univariate + homogeneous
  bivariate: `x²−y² → (x−y)(x+y)`, cubes, perfect-square trinomials),
  Cardano cubic, Ferrari quartic, `RootOf`, partial fractions, Gröbner basis.
- **Simplification** — `simplify` orchestrator (anti-bloat-guarded so it never
  returns a form larger than its input) chaining trigsimp (Pythagorean +
  double-angle / power-reduction `(1∓cos2x)/2 → sin²/cos²` / `2sin·cos → sin2x`
  collapses + double-angle ratio reduction `sin(2x)/sin(x) → 2cos(x)`),
  `expand_trig` (circular **and hyperbolic** angle-addition + multiple-angle
  `sin(nx)`/`sinh(nx)`), `fu`, powsimp, exp folds (`(eˣ)ⁿ → e^(nx)`,
  `eˣ+e⁻ˣ → 2cosh x`), hyperbolic double-angle (`2sinh·cosh → sinh2x`,
  `cosh²+sinh² → cosh2x`, `sinh2x/sinh x → 2cosh x`), imaginary-argument folds
  (`cos(I·x) → cosh x`, `|exp(I·x)| → 1`), `Abs`/`sign` identities
  (`sign(x)·|x| → x`, `|x|·|y| → |x·y|`), `log(√x) → ½log x`, radsimp, sqrtdenest,
  combsimp, gammasimp, cse, nsimplify, `together` over the LCM of denominators
  (incl. nested compound fractions `1/(1+1/x) → x/(x+1)`),
  `hyperexpand` (₀F₀→exp, ₁F₀, ₁F₁(1;2;z), ₂F₁(1,1;2;z), parameter cancellation).
- **Linear algebra** — det, inverse, eigendecomposition, LU/QR/
  Cholesky, rref / rank / nullspace, jacobian/hessian/wronskian,
  `Matrix::jordan_form()` (chains ≤ 2), `Matrix::exp(t)` matrix
  exponential via Jordan, `MatrixSymbol` expression tree.
- **Equation solvers** — `solve` (distinct roots; expands factored
  polynomials; radical equations `√x = 2`), `solveset` with `_invert` chain
  (peels log/exp/sin/cos/tan/sinh/cosh/tanh/abs and integer powers from the
  LHS; emits ImageSet over ℤ for periodic trig, including scaled arguments
  `cos(2x)=1`), `nsolve` / `vpasolve` (Newton in MPFR), inequality solver
  (exact endpoints, real `±∞`), `rsolve`, `nonlinsolve` via resultants and
  via Gröbner. **Set algebra** computes interval `∪` / `∩` / complement
  (`[1,3]∩[2,4] → [2,3]`, `ℝ\[1,3] → (−∞,1)∪(3,∞)`).
- **ODE / PDE** — `dsolve` for separable, linear, Bernoulli,
  exact, Riccati, homogeneous, autonomous, constant-coefficient,
  Cauchy-Euler, hypergeometric. Variation of parameters for
  nonhomogeneous 2nd-order linear (constant-coef and Cauchy-Euler).
  Linear systems via Jordan-form matrix exponential (handles
  defective A). PDE for first-order linear, heat, wave. IVP
  application + `checkodesol`.
- **Transforms** — Laplace (table + s-shift theorem: `sinh`/`cosh`,
  damped oscillations `e^(at)sin/cos`, `tⁿe^(at)`, and the multiplication-by-`tⁿ`
  rule `L{t·cos t} = (s²−1)/(s²+1)²`) and inverse Laplace
  (incl. completed-square quadratics `1/((s+1)²+1) → e^(-t)sin t`, linear
  numerators, and repeated irreducible quadratics
  `1/(s²+1)² → (sin t − t·cos t)/2`), Fourier, Mellin, Z, sine/cosine —
  forward and inverse, table-driven with linearity.
- **Units** — SI / CGS / US customary, prefixes, 12 physical
  constants with exact post-2019-redef values, affine
  temperature conversion.
- **Code generation** — C / C++ / Fortran / LaTeX / Octave
  printers + function emission.
- **MATLAB facade** — `sympp::matlab::*` namespace with `syms`,
  `sym(string)` / `str2sym` (parser-routed), `assume` /
  `assumeAlso` / `refresh`, `dsolve` / `linsolve` / `nsolve` /
  `vpasolve` / `pdsolve` overloads under MATLAB-friendly names.
- **SymPy oracle** — every numeric and structural assertion
  cross-checked against SymPy (1.13+) via a long-lived Python
  subprocess; 670 oracle-validated test cases. New `hyperexpand` op
  cross-validates SymPP's hypergeometric rewrites against SymPy's
  reference implementation.

## Documentation

- [docs/README.md](docs/README.md) — documentation index
- [docs/01-vision-and-scope.md](docs/01-vision-and-scope.md)
- [docs/02-architecture.md](docs/02-architecture.md)
- [docs/03-feature-mapping.md](docs/03-feature-mapping.md) — MATLAB → SymPy reference table
- [docs/04-roadmap.md](docs/04-roadmap.md) — phase-by-phase status + forward plan
- [docs/05-validation-strategy.md](docs/05-validation-strategy.md)
- [docs/06-build-and-tooling.md](docs/06-build-and-tooling.md)
- [docs/07-coding-standards.md](docs/07-coding-standards.md)
- [docs/08-tutorial.md](docs/08-tutorial.md) — **worked examples**
- [docs/09-known-issues.md](docs/09-known-issues.md) — SymPy-parity bug log
- [openspec/](openspec/) — spec-driven change specs (see
  `openspec/changes/` and `openspec/specs/`; `openspec list`)

## Build

### Prerequisites

- CMake ≥ 3.25
- C++20 compiler (gcc 11+, clang 14+, AppleClang 14+)
- GMP + GMPXX (`brew install gmp` / `apt install libgmp-dev libgmpxx4ldbl`)
- MPFR (`brew install mpfr` / `apt install libmpfr-dev`)
- Python 3.10+ with SymPy 1.13+ for the test oracle (`pip install sympy`)

CMake fetches Catch2 v3.5.4 and nlohmann/json v3.11.3 automatically.

### Build & test

```bash
git clone https://github.com/leonardoaraujosantos/SymPP.git
cd SymPP
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/tests/sympp_tests
```

The test runner spawns a Python subprocess on demand for the oracle —
SymPy must be importable from your `python3`. Tests with the `[oracle]`
Catch2 tag exercise that path; everything else runs without Python.

### Use as a dependency

```bash
cmake --install build --prefix /your/prefix
```

```cmake
find_package(SymPP REQUIRED)
target_link_libraries(your_target PRIVATE SymPP::sympp)
```

`SymPPConfig.cmake` and the per-target export are installed by the
build's install rules. A worked consumer build (with `demo.cpp` +
`CMakeLists.txt`) lives in the CI workflow's `find_package consumer
build` step.

## License

BSD 3-Clause. See [LICENSE](LICENSE). Matches SymPy upstream.

## Repository

https://github.com/leonardoaraujosantos/SymPP
