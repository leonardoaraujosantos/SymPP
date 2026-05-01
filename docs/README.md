# SymPP — Documentation

**SymPP** is a modern C++ symbolic mathematics library aiming for feature parity with the MATLAB Symbolic Math Toolbox (excluding Simulink, Live Editor, MATLAB Coder/Compiler integration).

The implementation is a clean-room port informed by SymPy's algorithms, with SymPy itself wired in as the validation oracle for every numeric and structural test.

## Index

| # | Document | Purpose |
|---|---|---|
| 01 | [Vision and Scope](01-vision-and-scope.md) | What we're building, what we're not, success criteria |
| 02 | [Architecture](02-architecture.md) | Class hierarchy, memory model, hash-consing, dispatch |
| 03 | [Feature Mapping](03-feature-mapping.md) | MATLAB function → SymPy source file reference |
| 04 | [Roadmap](04-roadmap.md) | 16-phase delivery plan with dependencies |
| 05 | [Validation Strategy](05-validation-strategy.md) | SymPy oracle harness, test taxonomy |
| 06 | [Build and Tooling](06-build-and-tooling.md) | CMake layout, dependencies, CI |
| 07 | [Coding Standards](07-coding-standards.md) | C++20 conventions, style, error model |

## Status

Pre-development. See [Roadmap](04-roadmap.md) — Phase 0 (foundation + oracle harness) is the next concrete deliverable.

## Reference checkout

The SymPy source is checked out at `sympy/` in the repo root for direct algorithm reading. Treat it as **reference-only**: never link against it, never copy code verbatim. All ports are clean-room reimplementations in modern C++.
