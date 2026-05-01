# SymPP

[![CI](https://github.com/leonardoaraujosantos/SymPP/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/leonardoaraujosantos/SymPP/actions/workflows/ci.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-00599C?logo=cplusplus&logoColor=white)](https://en.cppreference.com/w/cpp/20)
[![CMake](https://img.shields.io/badge/CMake-3.25%2B-064F8C?logo=cmake&logoColor=white)](https://cmake.org/)
[![Tests](https://img.shields.io/badge/tests-962%20passing-brightgreen)](#)
[![Oracle](https://img.shields.io/badge/oracle-SymPy%201.14-3B5526?logo=python&logoColor=white)](https://www.sympy.org/)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)](#)
[![Last commit](https://img.shields.io/github/last-commit/leonardoaraujosantos/SymPP)](https://github.com/leonardoaraujosantos/SymPP/commits/main)

Modern C++20 symbolic mathematics library. Clean-room port of SymPy
algorithms with SymPy itself wired in as the validation oracle.

## Status

```
962 tests / 1872 assertions  all passing
14 of 15 phases shipped
```

| Phase | Title | Status |
|---|---|---|
| 0  | Foundation & oracle harness            | ✅ |
| 1  | Core expression tree                   | ✅ |
| 2  | Assumptions                            | 🟡 minimal subset |
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

## What's in the box

- **Symbolic algebra** — Add/Mul/Pow with full canonical form,
  hash-cons cache, structural equality.
- **Number tower** — `Integer` / `Rational` (GMP) / `Float` (MPFR
  arbitrary-precision) / `ImaginaryUnit` / `Pi` / `E`.
- **30+ named functions** — sin/cos/tan/exp/log/sqrt/abs/floor/
  factorial/gamma/erf/heaviside/dirac_delta plus `Hyper` and
  `MeijerG` (proper Function classes with auto-eval) and the rest of
  the elementary + special + combinatorial canon.
- **Calculus** — `diff`, `integrate` (table + trig + parts +
  rational + heurisch), `series`, `limit` with L'Hôpital,
  `summation`, Padé, Euler-Lagrange, asymptotes.
- **Polynomials** — div/gcd/sqf, factor over ℤ, Cardano cubic,
  Ferrari quartic, `RootOf`, partial fractions, Gröbner basis.
- **Simplification** — `simplify` orchestrator chaining trigsimp
  (Pythagorean + Fu double-angle / 2sin·cos collapses),
  `expand_trig`, `fu`, powsimp, radsimp, sqrtdenest, combsimp,
  gammasimp, cse, nsimplify, `hyperexpand` (₀F₀→exp, ₁F₀, ₁F₁(1;2;z),
  ₂F₁(1,1;2;z), parameter cancellation).
- **Linear algebra** — det, inverse, eigendecomposition, LU/QR/
  Cholesky, rref / rank / nullspace, jacobian/hessian/wronskian,
  `Matrix::jordan_form()` (chains ≤ 2), `Matrix::exp(t)` matrix
  exponential via Jordan, `MatrixSymbol` expression tree.
- **Equation solvers** — `solve`, `solveset` with `_invert` chain
  (peels log/exp/sin/cos/tan/sinh/cosh/tanh/abs from the LHS;
  emits ImageSet over ℤ for periodic trig), `nsolve` / `vpasolve`
  (Newton in MPFR), inequality solver, `rsolve`, `nonlinsolve` via
  resultants and via Gröbner.
- **ODE / PDE** — `dsolve` for separable, linear, Bernoulli,
  exact, Riccati, homogeneous, autonomous, constant-coefficient,
  Cauchy-Euler, hypergeometric. Variation of parameters for
  nonhomogeneous 2nd-order linear (constant-coef and Cauchy-Euler).
  Linear systems via Jordan-form matrix exponential (handles
  defective A). PDE for first-order linear, heat, wave. IVP
  application + `checkodesol`.
- **Transforms** — Laplace, Fourier, Mellin, Z, sine/cosine —
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
  cross-checked against SymPy 1.14 via a long-lived Python
  subprocess. New `hyperexpand` op cross-validates SymPP's
  hypergeometric rewrites against SymPy's reference implementation.

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

## Build

### Prerequisites

- CMake ≥ 3.25
- C++20 compiler (gcc 11+, clang 14+, AppleClang 14+)
- GMP + GMPXX (`brew install gmp` / `apt install libgmp-dev libgmpxx4ldbl`)
- MPFR (`brew install mpfr` / `apt install libmpfr-dev`)
- Python 3.10+ with SymPy 1.14+ for the test oracle (`pip install sympy`)

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
