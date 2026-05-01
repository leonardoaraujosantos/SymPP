# 07 — Coding Standards

Modern C++ conventions for SymPP. Pragmatic over dogmatic. The point is consistency and readability across a large codebase that will be maintained for years.

## Language & compiler

- **C++20 minimum.** Concepts, `<ranges>`, `std::span`, `std::format` (or fmtlib fallback), `consteval`, `[[likely]]`/`[[unlikely]]`.
- C++23 features used **only with feature-test guards** (e.g. `std::expected` falls back to `tl::expected` shim).
- Supported compilers: gcc ≥ 12, clang ≥ 16, MSVC ≥ 19.36 (VS 2022 17.6).

## File layout

- One primary type per header. Free helpers may share a header if tightly related.
- Header guards: `#pragma once`. No `#ifndef` guards.
- Include order, blank-line separated:
  1. Corresponding header (in `.cpp`).
  2. C system headers.
  3. C++ standard library.
  4. Third-party (GMP, MPFR, fmtlib).
  5. SymPP public headers.
  6. SymPP internal headers (paths starting with `src/`).
- Forward-declare in headers; full include in `.cpp`.
- Every header self-contained (must compile alone).

## Naming

| Element | Convention | Example |
|---|---|---|
| Namespace | `lower_snake` | `sympp::core`, `sympp::matlab` |
| Type | `UpperCamel` | `Basic`, `Symbol`, `RationalDomain` |
| Function | `lower_snake` | `make_symbol`, `as_numer_denom` |
| Variable | `lower_snake` | `arg_count` |
| Member variable | trailing `_` | `args_`, `cache_` |
| Constant | `kUpperCamel` | `kPi`, `kZero` |
| Template parameter | `UpperCamel` | `template <typename Expr>` |
| Concept | `UpperCamel` | `ExpressionLike` |
| Macro (rare) | `SYMPP_UPPER_SNAKE` | `SYMPP_EXPORT`, `SYMPP_LIKELY` |
| File | `lower_snake.hpp` / `.cpp` | `basic.hpp`, `manual_integrate.cpp` |

Avoid Hungarian, `m_`-prefixed members, hungarian on function args. Don't shorten beyond recognition (`expression`, not `expr` for a member; `expr` is fine for a local).

## Memory & ownership

- Public surface: **`Expr` (= `std::shared_ptr<const Basic>`) by value.** Cheap, safe.
- Internals: `const Basic&` for read-only access without bumping the ref-count.
- `std::unique_ptr<T>` for sole-ownership internals (rarely needed; most things go in the cache).
- **No raw owning pointers.** Raw non-owning pointers OK for hot paths after profiling shows shared_ptr cost matters.
- **No manual `new`/`delete`.** Construction via factory `make<T>(...)` which routes through the hash-cons cache.

## Const-correctness

- All `Basic` methods that don't mutate are `const`.
- All node fields are `const` after construction (immutability invariant).
- `mutable` only on cache structures and rarely on memoized computations (e.g. cached hash).

## Error handling

| Situation | Mechanism |
|---|---|
| User-supplied invalid input (parse, undefined op) | Exception derived from `sympp::SymPPError`. |
| Algorithm cannot find closed form (e.g. `integrate` returns unevaluated) | Return unevaluated form (`Integral(...)`), **not** an error. |
| Internal invariant violation | `SYMPP_ASSERT(cond)` — `assert` in debug, undefined behaviour in release. Don't validate things that can't happen. |
| Hot-path expected-failure path | `std::expected<Expr, ErrorCode>` (tl::expected shim on C++20). |

Exception hierarchy:
```cpp
class SymPPError : public std::runtime_error { /* ... */ };
class ParseError       : public SymPPError {};
class DomainError      : public SymPPError {};
class NotImplemented   : public SymPPError {};   // when feature truly unsupported
class AssumptionError  : public SymPPError {};
```

## Templates & concepts

- Use **concepts** at every templated boundary. Compile errors should point at the concept, not at template instantiation hell.
- Prefer concepts over SFINAE.
- Hide template implementation in `.tpp`/`-inl.hpp` files included from headers if the implementation is non-trivial.

## Visitor pattern

- Public `Basic::accept(Visitor&) const` is virtual.
- Concrete visitors derive from `Visitor` and override per-node `visit(...)` methods.
- For double-dispatch over arithmetic (e.g. `Number * Number`), use a single dispatch table keyed on `(type_id_a, type_id_b)`.

## RAII and resource management

- All resources (mpz/mpq/mpfr handles, file handles, JIT modules) wrapped in RAII types.
- No `init`/`finalize` patterns at namespace scope. Static initialization order fiasco is avoided by Meyer-singleton (`static T& instance()`).

## Threading

- Expression nodes are immutable → freely shared across threads.
- The expression cache is the only mutable shared state; it uses sharded mutex or `tbb::concurrent_hash_map`.
- No global mutable state besides the cache and a debug-print toggle.
- Algorithms are pure functions on `Expr`. Public algorithms must be thread-safe.

## Performance

- **Profile before optimizing.** Premature optimization in a CAS is especially dangerous because algorithmic improvements dominate constant factors.
- **Hot paths: no exceptions, no virtuals.** Use type-tag switch.
- **Cold paths: prefer clarity.** Most code is cold.
- Mark hot single-translation-unit helpers `inline` or `[[gnu::always_inline]]` where measured to matter.
- `[[likely]]` / `[[unlikely]]` annotations only with profile-driven evidence.
- Avoid premature `constexpr` proliferation. Use it where it actually helps users (singletons, simple algebra).

## Comments and documentation

- **Default to no comments.** Well-named identifiers do the job.
- Comments only when the *why* is non-obvious: a workaround, a subtle invariant, a citation to a paper or SymPy file you ported.
- Public headers: brief one-paragraph doc per type, one line per public method. Doxygen-friendly but no `@param`/`@return` ceremony unless adding information beyond signature.
- For algorithms ported from SymPy: top-of-function comment naming the SymPy file. Example:
  ```cpp
  // Ported from sympy/integrals/heurisch.py::heurisch
  // Reference: Bronstein, "Symbolic Integration I" §5
  Expr heurisch(const Expr& integrand, const Symbol& var) {
      ...
  }
  ```

## Logging

- No `std::cout` in library code.
- Optional debug logging via a compile-time-disabled `SYMPP_DEBUG_LOG(...)` macro.
- Tests may print on failure via Catch2's facilities.

## Tests

- One Catch2 file per source file: `src/foo.cpp` ↔ `tests/foo_test.cpp`.
- Each test tagged: `[<phase>][<feature>][oracle?]`. Examples: `[1][add]`, `[7][integrate][oracle]`.
- Use `SECTION` to share setup across related cases.
- Oracle tests use `ASSERT_MATCHES_SYMPY` macro from [05-validation-strategy](05-validation-strategy.md).

## Forbidden / discouraged

- `using namespace std;` — never in headers, avoid in .cpp.
- C-style casts — use `static_cast` / `dynamic_cast` / `reinterpret_cast`.
- `goto` — except in cleanup-loop refactors that can't reasonably be expressed otherwise.
- Macros for code generation (function-like) — prefer templates and concepts.
- Static lifetime mutable globals — except the expression cache.
- Output parameters via non-const reference — return tuples/structs instead. Const-ref is fine for read-only.
- Implicit conversions on user-defined types — mark single-arg constructors `explicit`.

## Deprecation

When a public API changes:
1. Mark old API `[[deprecated("use X instead")]]`.
2. Keep for at least one minor version.
3. Document in CHANGELOG.
4. Remove on next major.

## Commit & PR conventions

(See user's global git rules for commit messages.)

- One logical change per commit.
- Tests + implementation in the same commit when introducing a feature.
- PRs reference the phase from [04-roadmap](04-roadmap.md): `[Phase 7] Add manual integration rule for trigonometric powers`.
- Every PR introducing a new public API includes oracle tests and at least one example in `examples/`.
