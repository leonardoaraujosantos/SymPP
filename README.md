# SymPP

Modern C++ symbolic mathematics library targeting feature parity with the MATLAB Symbolic Math Toolbox, validated continuously against SymPy.

## Status

Pre-development — Phase 0 (foundation + oracle harness).

See [docs/](docs/) for the full plan:

- [01 — Vision and Scope](docs/01-vision-and-scope.md)
- [02 — Architecture](docs/02-architecture.md)
- [03 — Feature Mapping (MATLAB → SymPy source)](docs/03-feature-mapping.md)
- [04 — Roadmap](docs/04-roadmap.md)
- [05 — Validation Strategy](docs/05-validation-strategy.md)
- [06 — Build and Tooling](docs/06-build-and-tooling.md)
- [07 — Coding Standards](docs/07-coding-standards.md)

## Quick start (Phase 0)

```bash
# macOS prerequisites
brew install cmake python@3.12
pip3 install sympy

# Build and test
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Layout

```
SymPP/
├── CMakeLists.txt
├── include/sympp/      — public headers
├── src/                — library implementation
├── tests/              — Catch2 test suite (oracle-validated)
├── oracle/             — SymPy oracle Python harness
├── examples/           — sample consumer projects
├── benchmarks/         — performance suite
├── cmake/              — CMake helpers
├── docs/               — design documents
└── sympy/              — SymPy reference checkout (read-only, gitignored)
```

## License

BSD 3-Clause. See [LICENSE](LICENSE).
