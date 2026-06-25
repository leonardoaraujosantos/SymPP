## Why

Day-to-day SymPP work — configure, compile, run unit vs. integration tests,
run the examples, the parity gate, benchmarks and docs — was an ad-hoc set of
`cmake` and Catch2 invocations that contributors had to remember. A single,
discoverable task runner makes the standard workflows one command each and
documents the unit/integration split (the `[oracle]` Catch2 tag) in one place.

## What Changes

- Add a `justfile` ([just](https://github.com/casey/just)) at the repository
  root exposing the common developer tasks as named recipes.
- The recipes cover: `configure`, `build`/`compile`/`rebuild`, `test`,
  `test-unit` (pure C++, `~[oracle]`), `test-integration` (SymPy-oracle,
  `[oracle]`), `test-filter`, `examples`, `example <name>`, `parity` (the
  Tier 1 + Tier 2 acceptance probe), `bench`, `docs`, `install`, `format`,
  `spec`/`spec-validate`, `ci`, and `clean`.
- Recipes are portable: defaults (`build_dir`, `build_type`) are overridable
  inline, and `parity` reads the GMP/MPFR paths from the CMake cache so it
  builds the probe on any host.
- Document the runner in `README.md` under Build.

## Capabilities

- **developer-tooling** (new): a `just` task runner standardizing build, the
  unit/integration test split, example execution, and the parity gate.

## Non-goals

- Replacing CMake — `just` only orchestrates the existing CMake build and the
  Catch2 runner; the canonical build remains `cmake`.
