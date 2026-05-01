# SymPP — Documentation

**SymPP** is a modern C++20 symbolic mathematics library. Originally
scoped to feature parity with the MATLAB Symbolic Math Toolbox, the
project's current target is broader: a real C++ alternative to SymPy.

The implementation is a clean-room port informed by SymPy's algorithms,
with SymPy itself wired in as the validation oracle for every numeric
and structural test.

## Index

| # | Document | Purpose |
|---|---|---|
| 01 | [Vision and Scope](01-vision-and-scope.md)         | What we're building, what we're not, success criteria |
| 02 | [Architecture](02-architecture.md)                  | Class hierarchy, memory model, hash-consing, dispatch |
| 03 | [Feature Mapping](03-feature-mapping.md)            | MATLAB function → SymPy source file reference |
| 04 | [Roadmap](04-roadmap.md)                            | Per-phase status (shipped / deferred) + forward plan |
| 05 | [Validation Strategy](05-validation-strategy.md)    | SymPy oracle harness, test taxonomy |
| 06 | [Build and Tooling](06-build-and-tooling.md)        | CMake layout, dependencies, CI |
| 07 | [Coding Standards](07-coding-standards.md)          | C++20 conventions, style, error model |
| 08 | [Tutorial](08-tutorial.md)                          | Worked examples covering each shipped phase |

## Status

14 of 15 phases shipped (Phase 16 partial). 962 tests / 1872
assertions, 307 oracle-validated against SymPy 1.14, all passing.
See [Roadmap](04-roadmap.md) for the
shipped/deferred breakdown per phase and the path to full SymPy
parity (Phases 17-24 cover modules originally listed as anti-scope:
vector calculus, tensors, statistics, number theory, combinatorics,
physics submodules, cryptography, geometry).

## Reference checkout

The SymPy source is checked out at `sympy/` in the repo root for
direct algorithm reading. Treat it as **reference-only**: never link
against it, never copy code verbatim. All ports are clean-room
reimplementations in modern C++.
