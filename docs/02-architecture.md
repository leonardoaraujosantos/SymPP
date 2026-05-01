# 02 — Architecture

## Core principles

1. **Immutability.** Expression nodes are immutable after construction. All "modifications" return a new node. Aligns with how mathematical expressions actually work and enables safe sharing.
2. **Hash-consing.** Identical subtrees share storage via a global expression cache. Critical for performance (SymPy uses this implicitly via `__hash__` + interning; we do it explicitly).
3. **Value semantics on top of `shared_ptr`.** Users pass `Expr` by value; internally it's a `shared_ptr<const Basic>`. Cheap copies, no manual lifetime management.
4. **Visitor-driven transforms.** Each algorithm is a visitor; nodes don't carry algorithm logic.
5. **No exceptions in tight loops.** Errors at API boundaries via exceptions; internal hot paths use `std::expected` (C++23) or `tl::expected` shim (C++20).
6. **Compile-time dispatch where possible.** Concepts + CRTP for hot paths; virtual dispatch for the public node hierarchy.

## Type hierarchy

```
Basic                                    (abstract: hash, equality, args, type tag)
├── Expr                                 (abstract: arithmetic-bearing nodes)
│   ├── AtomicExpr                       (abstract: leaves)
│   │   ├── Symbol                       ("x", with assumptions)
│   │   ├── Dummy                        (unique-per-construction symbol)
│   │   ├── Wild                         (pattern-matching wildcard)
│   │   ├── Number                       (abstract)
│   │   │   ├── Integer                  (GMP mpz)
│   │   │   ├── Rational                 (GMP mpq)
│   │   │   ├── Float                    (MPFR mpfr_t)
│   │   │   └── ComplexBase              (a + b*I where a,b are Number)
│   │   └── NumberSymbol                 (singletons: Pi, E, EulerGamma, GoldenRatio, …)
│   ├── Add                              (commutative; canonical args)
│   ├── Mul                              (commutative; canonical args)
│   ├── Pow                              (binary: base, exp)
│   ├── Function                         (abstract; named function applied to args)
│   │   ├── ElementaryFunction           (sin, cos, exp, log, sqrt, abs, sign, …)
│   │   ├── SpecialFunction              (gamma, beta, erf, Bessel, Mathieu, …)
│   │   └── UndefinedFunction            (user f(x))
│   ├── Derivative                       (unevaluated d/dx form)
│   ├── Integral                         (unevaluated ∫ form)
│   ├── Sum / Product                    (unevaluated Σ / Π)
│   ├── Limit                            (unevaluated lim form)
│   ├── Order                            (Big-O / Landau)
│   ├── Piecewise                        (cases)
│   └── MatrixBase / MatrixExpr          (matrix-valued expressions)
├── Boolean                              (abstract; for assumptions and Piecewise conditions)
│   ├── BooleanAtom (True, False)
│   ├── And / Or / Not / Xor / Implies / Equivalent
│   └── Relational                       (Eq, Ne, Lt, Le, Gt, Ge)
├── Set                                  (abstract; for solveset, intervals, unions)
│   ├── Interval, FiniteSet, Union, Intersection, Complement, …
└── Tuple / List / Dict                  (containers as Basic; needed for solver outputs)
```

### Public alias

```cpp
namespace sympp {
    using Expr = std::shared_ptr<const Basic>;
}
```

Users always handle `Expr`. Concrete node types are construction details.

## Hash-consing

Every `Basic` subclass produces a structural hash combining `type_tag` + arg hashes + (for atomic types) value hash. A global thread-safe cache (concurrent hash map, e.g. `tbb::concurrent_hash_map` or a sharded `std::unordered_map`) deduplicates subtrees.

```cpp
template <typename T, typename... Args>
Expr make(Args&&... args) {
    auto candidate = std::make_shared<T>(std::forward<Args>(args)...);
    return ExprCache::instance().intern(std::move(candidate));
}
```

`ExprCache::intern` looks up the candidate's hash; returns the existing pointer if a structurally-equal node exists, otherwise inserts the candidate. Equality check is structural.

**Lifetime:** cache holds `weak_ptr`. When the last external `shared_ptr` drops, the cache slot is reclaimed lazily.

## Canonical form

Commutative operations canonicalize arguments at construction:
- `Add` and `Mul` flatten nested same-type args, combine numeric coefficients, and sort symbolic args by an internal stable order (type-tag, then printed-name, then hash).
- `Pow` applies basic auto-evaluation: `x^0 → 1`, `x^1 → x`, `1^x → 1`, `0^x → 0` (under assumptions), `(x^a)^b → x^(a*b)` (under assumptions).
- `Mul` extracts numeric coefficient: `2*x*3*y → 6*x*y`.

This is the canonicalization layer; semantic simplification (`expand`, `simplify`, `trigsimp`) is separate and explicit.

## Number tower

| Type | Backend | Notes |
|---|---|---|
| `Integer` | GMP `mpz_t` | Arbitrary precision. Cached singletons for {-1, 0, 1, 2, …, 32}. |
| `Rational` | GMP `mpq_t` | Auto-reduced; denominator > 0. `Rational(a, 1)` returns `Integer`. |
| `Float` | MPFR `mpfr_t` | Variable precision. Backs `vpa` and `evalf(prec)`. |
| `Complex` | Pair of `Number` | If imaginary part is zero, returns the real part; collapses to `Number`. |

### Coercion ladder
`Integer ⊂ Rational ⊂ Float ⊂ Complex`. Mixed-arithmetic rules in a single dispatch table.

### Singletons
`Zero`, `One`, `NegativeOne`, `Half`, `Infinity`, `NegativeInfinity`, `ComplexInfinity`, `NaN`, `Pi`, `E`, `EulerGamma`, `Catalan`, `GoldenRatio`, `ImaginaryUnit`. Each is a `Number` or `NumberSymbol` with a stable address.

## Visitors

```cpp
class Visitor {
public:
    virtual ~Visitor() = default;
    virtual void visit(const Symbol&) = 0;
    virtual void visit(const Integer&) = 0;
    // ... one overload per concrete node type
};
```

Algorithms (`expand`, `diff`, `subs`, every printer, `expand_log`, `trigsimp` rules) are visitors. `Basic::accept(Visitor&)` is the dispatch hook.

For high-frequency operations (hash, equality), we use a tag-based switch over `Basic::type_id()` to avoid virtual call overhead.

## Module layout

> **Note**: the tree below is the *target* layout — the design intent
> for SymPy parity. The current `include/sympp/` layout differs in
> places: `ode/dsolve.hpp` (not `solvers/ode.hpp`), no separate
> `polys/factor.hpp` or `polys/domain.hpp` yet (factor lives in
> `polys/poly.hpp`), `simplify/hyperexpand.hpp` is shipped,
> `functions/hypergeometric.hpp` is shipped (`Hyper`, `MeijerG`),
> sparse matrices and `codegen/` are not yet present, and the
> `matlab/` directory has been split into per-area sub-headers
> (`parsing.hpp`, `assumptions.hpp`, `solvers.hpp`, `ode.hpp`,
> `pde.hpp`) under the umbrella `matlab.hpp`.

```
sympp/
├── include/sympp/
│   ├── sympp.hpp                 — umbrella header
│   ├── core/
│   │   ├── basic.hpp             — Basic, Expr alias
│   │   ├── symbol.hpp            — Symbol, Dummy, Wild
│   │   ├── number.hpp            — Number tower (Integer, Rational, Float, Complex)
│   │   ├── add.hpp / mul.hpp / pow.hpp
│   │   ├── function.hpp          — Function base, UndefinedFunction
│   │   ├── operators.hpp         — operator+, -, *, /, ==, etc.
│   │   ├── subs.hpp / diff.hpp / expand.hpp
│   │   └── assumptions.hpp
│   ├── functions/
│   │   ├── elementary.hpp        — sin, cos, exp, log, sqrt, abs, …
│   │   └── special.hpp           — gamma, erf, Bessel, …
│   ├── polys/
│   │   ├── poly.hpp / domain.hpp
│   │   ├── factor.hpp / gcd.hpp / partfrac.hpp
│   │   └── roots.hpp
│   ├── simplify/
│   │   ├── simplify.hpp / trigsimp.hpp / radsimp.hpp / powsimp.hpp
│   │   └── cse.hpp
│   ├── calculus/
│   │   ├── series.hpp / limit.hpp / summation.hpp
│   │   └── pade.hpp / euler.hpp
│   ├── integrals/
│   │   ├── integrate.hpp         — orchestrator
│   │   ├── manual.hpp            — rule-based ports of manualintegrate.py
│   │   ├── rational.hpp          — Lazard-Rioboo-Trager
│   │   ├── heurisch.hpp          — Risch heuristic
│   │   └── transforms.hpp        — fourier, laplace, ztrans
│   ├── solvers/
│   │   ├── solve.hpp / nsolve.hpp / linsolve.hpp
│   │   ├── ode.hpp / pde.hpp / recurrence.hpp
│   ├── matrices/
│   │   ├── matrix.hpp            — dense
│   │   ├── sparse.hpp
│   │   ├── decompositions.hpp    — LU, QR, SVD, Cholesky, Jordan
│   │   └── expressions.hpp       — MatrixSymbol, MatMul, MatAdd
│   ├── sets/                     — Interval, FiniteSet, Union, …
│   ├── physics/units.hpp
│   ├── printing/
│   │   ├── str.hpp / latex.hpp / pretty.hpp
│   │   └── ccode.hpp / fcode.hpp / octave.hpp
│   ├── codegen/
│   │   ├── ast.hpp / cnodes.hpp / fnodes.hpp
│   │   └── lambdify.hpp          — JIT via LLVM ORC
│   ├── parsing/parse.hpp
│   └── matlab/                   — MATLAB-named facade: syms, solve, simplify, …
├── src/                          — implementations
├── tests/                        — Catch2 tests, oracle-validated
├── oracle/                       — SymPy oracle Python harness
├── examples/                     — sample consumer projects
├── benchmarks/                   — perf vs SymEngine
├── cmake/                        — Find/Config modules
└── docs/                         — this directory
```

## Error model

- **API boundary errors** (parse failure, undefined operation, domain mismatch): exceptions derived from `sympp::SymPPError`.
- **Algorithm "not implemented" cases** (e.g. integral SymPy itself returns unevaluated): return the unevaluated form (`Integral`, `Limit`, `Sum`) — not an error.
- **Internal invariant violations**: assertions; in release builds, undefined behavior. Don't validate what can't happen.
- **Hot paths**: `std::expected<Expr, ErrorCode>` to avoid throw cost.

## Concurrency

- Expression cache is thread-safe (sharded mutex or concurrent hash map).
- Expressions themselves are immutable, so reads are lock-free.
- Algorithms are pure functions on `Expr`; safe to parallelize across expressions.
- No global mutable state outside the cache and a debug-print toggle.

## Memory model

- Per-expression: `shared_ptr` overhead (16 bytes) + node fields. Atomic ref-count costs ~10ns per copy.
- Cache: bounded by program working set; expression nodes drop when unreferenced.
- Numbers: GMP/MPFR allocate on heap; `Integer` for small values uses a singleton table to avoid allocation entirely.

## Compile-time concepts

```cpp
template <typename T>
concept ExpressionLike = std::derived_from<T, Basic> && requires(T t) {
    { t.hash() } -> std::convertible_to<std::size_t>;
    { t.args() } -> std::ranges::range;
};
```

Used to constrain templated visitors and helpers.
