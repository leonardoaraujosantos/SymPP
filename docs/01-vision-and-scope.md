# 01 — Vision and Scope

## Mission

Build a modern, idiomatic C++ symbolic mathematics library — **SymPP** — that any C++ project can drop in via CMake `find_package`. The original target was MATLAB Symbolic Math Toolbox parity; the current target is broader — a real C++ alternative to SymPy itself, with SymPy's algorithm catalog ported and validated against SymPy as the oracle.

## Why

The C++ ecosystem has gaps:
- **GiNaC** — solid algebra, weak on simplification, no transforms, no ODE solver, dated C++.
- **SymEngine** — fast core, but a thin slice of features compared to a full CAS; designed primarily as a backend.
- **SymPy** — feature-rich, but Python-only, slow, GIL-bound, ships a runtime.
- **Mathematica/Maple/MATLAB** — closed source, license-locked.

We want SymPy's breadth, SymEngine's speed, and a license/runtime profile that lets it ship inside arbitrary C++ products (engineering, robotics, embedded simulation, optimization, custom DSLs).

## Goals

1. **MATLAB Symbolic Math Toolbox feature parity** for the math features documented in *Symbolic Math Toolbox User's Guide R2026a* — Chapters 1–5 except Simulink/Simscape integration. See [03-feature-mapping](03-feature-mapping.md). **Status: 13/16 phases shipped, on track.**
2. **SymPy parity** for the modules a research/engineering user expects: vector calculus, tensors, statistics, number theory, combinatorics, physics submodules, geometry. Originally listed as anti-scope; promoted to in-scope (Category C in the [Roadmap](04-roadmap.md)).
3. **Clean-room port from SymPy.** Read the SymPy algorithm; reimplement in modern C++. Never copy code verbatim — license-compatible (BSD-3) but we want idiomatic C++, not Python-shaped C++.
4. **SymPy as validation oracle.** Every algorithm has a test that runs the same input through SymPy and asserts equivalence. See [05-validation-strategy](05-validation-strategy.md).
5. **Modern C++.** C++20 minimum: concepts, `std::span`, `<ranges>`, `std::variant`, `std::format`. Header-only public API where viable; compiled internals.
6. **Deployable.** No Python runtime in production. CMake `find_package(SymPP)`, vcpkg/Conan recipes, single shared library.
7. **License: BSD-3.** Match SymPy. Compatible with commercial use.

## Non-goals

- **Not a SymPy or SymEngine wrapper.** SymPy is a reference and an oracle, not a runtime dependency of the shipped library. SymEngine is studied for its design but not vendored.
- **Not Python bindings first.** A pybind11 layer may come later for adoption, but the C++ API is the product.
- **Not Simulink/Simscape.** Proprietary block formats are out of scope. If we generate code, it will be C, C++, Fortran, LaTeX, or Octave-syntax (the latter is incidentally MATLAB-readable for `matlabFunction`).
- **Not the Live Editor.** That's an IDE feature, not a math feature.
- **Not a numeric ODE/PDE solver.** We solve *symbolically* (`dsolve`, `pdsolve` analytic methods). For numeric DAE/ODE the consumer wires their own integrator (Sundials/CasADi/Boost.Numeric.Odeint). We may emit code that drives such a solver.

## Target users

- C++ engineering codebases needing inline symbolic math (control systems, robotics IK/dynamics, signal processing, circuit simulation).
- Teams currently shelling out to MATLAB or Python for one-shot symbolic computations who want it native.
- Custom DSL/PL projects needing a CAS backend.
- Authors of larger numeric libraries who need symbolic Jacobians/Hessians/code-gen.

## Success criteria for v1.0

- A consumer can `find_package(SymPP REQUIRED)` and link.
- The 20-example regression suite (chosen from MATLAB toolbox worked examples — see [05-validation-strategy](05-validation-strategy.md)) passes with output equivalent to SymPy.
- Public headers documented, every function has at least one oracle-validated test.
- CI green on Linux + macOS + Windows, with ASan/UBSan/clang-tidy clean.
- Performance: within 5× of SymEngine on the SymEngine benchmark suite for shared operations; faster than SymPy on every shared operation.

## Explicit scope by chapter (MATLAB R2026a User's Guide)

| Chapter | In scope | Notes |
|---|---|---|
| 1 — Getting Started | ✅ Full | Symbol/Number/Matrix/Function/Assumptions API |
| 2 — Symbolic Computations | ✅ Mostly | Excluding Live Editor UI features (subscripts/copy-paste). Units full. |
| 3 — Mathematics | ✅ Full | Solvers, calculus, transforms, simplification, linear algebra, polynomials, prime factorization |
| 4 — Graphics | ❌ Out of scope | Consumer plots SymPP output via their renderer of choice (gnuplot/matplotlib/etc.) |
| 5 — Code Generation | ✅ Mostly | C/Fortran/LaTeX/Octave-syntax. No Simulink, no Simscape, no MATLAB Coder app. |
| 6 — Function Reference | ✅ Full | All math functions; MATLAB-named facade in `sympp::matlab::` |

## Anti-scope (genuinely out of scope)

- **Simulink / Simscape** — proprietary block formats, no math content.
- **MATLAB Live Editor** — IDE feature, not a math feature.
- **Plotting / rendering backends** — consumers pipe SymPP output (`ccode`, `latex`, `evalf`, `vpaintegral`) to their plotting layer of choice. A SymPP-specific plotting bridge would couple the library to a renderer without adding mathematical capability.
- **Numeric ODE/PDE solvers** — we solve symbolically. Numeric integration is the consumer's problem (Sundials / CasADi / Boost.Numeric.Odeint). We may emit code that drives them.
- **Domain-specific languages** beyond plain math notation. Categorical-style algebra packages (HOL, Coq-style proof assistants) are a different product class.

## Scope reclamation

The earlier draft of this document listed tensor algebra, crypto,
combinatorics, statistics, geometry, vector calculus, and mechanics
modules as anti-scope. Those have been **promoted to in-scope** —
they are explicit goals tracked as **Phases 17-24 in the
[Roadmap](04-roadmap.md)**. They land after the v0.5 release
(Phases 14-16) and the deep-deferred items within already-shipped
phases close out.
