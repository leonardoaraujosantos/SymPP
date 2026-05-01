# 05 — Validation Strategy

The single most important architectural decision in this project: **SymPy is the source of truth.** Every algorithm SymPP ships has at least one test that runs the same input through SymPy and asserts equivalence. Without this, regressions are invisible and we'll never catch up to a CAS that's had a thousand contributors over twenty years.

## The oracle harness

A **long-lived Python subprocess** running SymPy. The C++ test process talks to it via stdin/stdout (line-delimited JSON). Started once per test binary, reused across thousands of tests.

### Why subprocess (not embed CPython)

1. **Build-time isolation.** SymPP itself never links against Python. The oracle is a *test-only* dependency.
2. **No GIL contention** in test runners.
3. **Crash isolation.** If SymPy hits an internal error, our test process keeps going.
4. **Trivial setup.** Any developer with `python3 -m pip install sympy` can run the suite.

### Protocol

Line-delimited JSON. Each request is a single line; each response is a single line.

```jsonc
// C++ → Python
{"id": 42, "op": "simplify", "expr": "sin(x)**2 + cos(x)**2"}
{"id": 43, "op": "diff",     "expr": "x**3 + sin(x)", "var": "x"}
{"id": 44, "op": "evalf",    "expr": "pi",            "prec": 50}
{"id": 45, "op": "integrate","expr": "1/(x**2+1)",    "var": "x"}

// Python → C++ (success)
{"id": 42, "ok": true, "result": "1"}
{"id": 43, "ok": true, "result": "3*x**2 + cos(x)"}
{"id": 44, "ok": true, "result": "3.1415926535897932384626433832795028841971693993751"}

// Python → C++ (failure)
{"id": 45, "ok": false, "error": "TimeoutError", "detail": "exceeded 10s"}
```

### Equivalence checks

Two C++ Expr → SymPy expressions are "equivalent" by one of:

1. **Structural** — `srepr(a) == srepr(b)` after both go through `sympy.simplify`. Strict; catches reordering only when commutativity fails.
2. **Symbolic-zero** — `sympy.simplify(a - b) == 0`. The default check.
3. **Numeric** — `abs((a - b).evalf(50)) < 1e-40` with random substitutions for free symbols. Used when symbolic check times out or returns inconclusive.
4. **Set equality** — for solver outputs: `sympy.Set(a) == sympy.Set(b)` (handles different ordering of solutions).

C++ side picks the appropriate mode per test.

### C++ API

```cpp
namespace sympp::testing {

class Oracle {
public:
    static Oracle& instance();                           // singleton, starts subprocess on first use
    
    // Send request, get SymPy's answer as string.
    std::string ask(std::string_view op,
                    std::string_view expr,
                    const json& kwargs = {});
    
    // High-level equivalence check: simplify(a - b) == 0
    bool equivalent_to_sympy(const Expr& cpp_result,
                             std::string_view sympy_input,
                             std::string_view sympy_op,
                             const json& kwargs = {});
};

#define ASSERT_MATCHES_SYMPY(cpp_expr, sympy_input, sympy_op, ...) \
    do { \
        auto _r = (cpp_expr); \
        REQUIRE(::sympp::testing::Oracle::instance() \
            .equivalent_to_sympy(_r, sympy_input, sympy_op, ##__VA_ARGS__)); \
    } while (0)

}
```

### Catch2 integration

```cpp
TEST_CASE("diff: polynomial chain rule", "[calculus][oracle]") {
    auto x = symbol("x");
    auto e = pow(x, 3) + sin(x);
    auto d = diff(e, x);
    ASSERT_MATCHES_SYMPY(d, "x**3 + sin(x)", "diff", {{"var", "x"}});
}
```

## Test taxonomy

Five distinct test layers; each phase produces tests in all five.

### 1. Unit tests
Smallest unit: one function, one input, one expected output (literal). No oracle.
**Purpose:** catch regressions in known-good cases. Fast.

### 2. Oracle tests
C++ result vs SymPy result. The bulk of the suite.
**Purpose:** catch divergence from SymPy on any input we test.

### 3. Property tests
QuickCheck-style. Random expression generators. Properties checked:
- Round-trip: `parse(str(e)) == e`.
- Identity: `simplify(e) == e` modulo equivalence; `expand(factor(p)) == p` for polynomials; `diff(int(f, x), x) == f`.
- Algebraic: `(a + b) - b == a`; `a * 1 == a`; commutativity, associativity hold.
- vs SymPy: random inputs through both, results must match.
**Purpose:** find cases we'd never write by hand.

### 4. Golden tests
Snapshot-based. Run a curated example, write its serialized output to disk, future runs diff.
**Purpose:** catch unintended output changes (printer changes, simplification heuristic shifts).

### 5. MATLAB regression suite
The 20 worked examples chosen from the *Symbolic Math Toolbox User's Guide*. Each has:
- Original MATLAB input (preserved).
- Expected MATLAB output (captured once, by hand or by running MATLAB).
- SymPP equivalent expression.
- Pass criterion: SymPP output == MATLAB output up to equivalence.
**Purpose:** the actual product acceptance test.

## Selecting the MATLAB regression suite

Pick examples that span features and are non-trivial. Concrete starter list (page numbers from R2026a User's Guide):

| # | Topic | Page | Features exercised |
|---|---|---|---|
| 1 | Damped harmonic oscillator | 2-78 | dsolve, simplify, plot |
| 2 | Wind turbine power | 2-87 | symbolic ints, simplify |
| 3 | Image undistortion | 2-93 | matrix algebra, solve |
| 4 | Electric dipole | 2-99 | diff, vector calculus |
| 5 | Algebraic equations | 3-3 | solve |
| 6 | System of algebraic | 3-7 | linsolve / solve |
| 7 | ODE with IC | 3-49 | dsolve |
| 8 | ODE system | 3-53 | dsolve_system |
| 9 | PDE — tsunami | 3-59 | pdsolve |
| 10 | DAE manipulation | 3-87 | symbolic DAE prep |
| 11 | Inverse kinematics 2-link arm | 3-94 | jacobian, solve |
| 12 | Black-Scholes | 3-118 | symbolic + evalf |
| 13 | partfrac, expand, factor | 3-123 | expression rearrangement |
| 14 | Harmonic analysis | 3-160 | Laplace, simplify |
| 15 | Differentiation | 3-181 | diff (high order, partial) |
| 16 | Integration with complex params | 3-189 | int |
| 17 | Taylor series | 3-192 | taylor |
| 18 | Fourier transform — beam | 3-194 | fourier |
| 19 | Z-transform difference eq | 3-205 | ztrans |
| 20 | Padé of time-delay | 3-261 | pade |

Each example is captured as:
```
tests/matlab_regression/01_damped_harmonic/
├── input.m              — original MATLAB
├── expected.txt         — captured MATLAB output
├── sympp.cpp            — SymPP equivalent
└── README.md            — what's being tested
```

## Performance baselines

Two parallel test passes:
- **Correctness pass** — Catch2 + oracle. Slow OK; runs in CI on every push.
- **Performance pass** — Google Benchmark. SymPP vs SymEngine on shared operations. Runs nightly. Regression budget: 10% on main, alert on Slack.

## Coverage targets per phase

| Phase | Min oracle tests | Notes |
|---|---|---|
| 0 | 5 | Smoke. |
| 1 | 200+ | Arithmetic + subs + evalf coverage. |
| 2 | 50+ | Per assumption × per node type. |
| 3 | 1000+ | Each function × {eval, evalf, diff, rewrite, series}. |
| 4 | 500+ | Polynomial ops × domain matrix. |
| 5 | 300+ | Standard simplification patterns. |
| 6 | 300+ | Limits (Gruntz hard cases), Taylor, sums. |
| 7 | 500+ | Wester integration test set. |
| 8 | 200+ | Standard transform tables. |
| 9 | 200+ | Matrix operations, decompositions. |
| 10 | 200+ | Equation classes. |
| 11 | 100+ | ODE classifier coverage. |
| 12 | 50+ | Unit conversions. |
| 13 | 100+ | Round-trip codegen → execute → match. |
| 14–15 | 50+ | End-to-end. |
| 16 | regression sweep | All phases combined. |

Total expected: **~3500+ oracle tests** at v1.0.

## CI matrix

```yaml
matrix:
  os: [ubuntu-latest, macos-latest, windows-latest]
  compiler: [gcc, clang, msvc]   # filtered per OS
  build_type: [Debug, Release]
  sanitizer: [none, asan, ubsan]
  python: ["3.11"]
  sympy: ["1.13", "latest"]      # detect SymPy version drift
```

Every push: full matrix. Nightly: include `sympy:latest` and report new divergences.

## What we don't oracle-validate

- **Performance** — different metric.
- **Memory layout / cache behaviour** — opaque to SymPy.
- **API ergonomics** — design judgment.
- **Cases SymPy gets wrong** — we maintain a `known_sympy_bugs.json` allowlist with rationale; tests in this list assert against our chosen correct answer, not SymPy.

## What if SymPy disagrees with itself across versions?

Pin a baseline SymPy version (e.g. 1.13). The CI runs against `latest` only as a *warning*, not a failure. Drift goes into `sympy_drift.md` for review.

## Bootstrapping order

The first thing built in Phase 0 is the oracle harness — *before* any algorithm code. The first test ever written is "create Symbol x, send to oracle, receive `Symbol('x')` back". Everything else builds on this.
