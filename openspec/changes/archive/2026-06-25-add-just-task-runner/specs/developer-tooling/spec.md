# developer-tooling Specification

## ADDED Requirements

### Requirement: just task runner

The project SHALL provide a `justfile` at the repository root that exposes the
standard developer workflows as named recipes over the existing CMake build and
Catch2 test runner. `just` with no arguments SHALL list the available recipes.
Recipes SHALL accept inline overrides of `build_dir` and `build_type`.

#### Scenario: List recipes
- **WHEN** `just` is run with no arguments
- **THEN** the available recipes (build, test, test-unit, test-integration,
  examples, parity, …) are listed

#### Scenario: Build
- **WHEN** `just build` is run
- **THEN** CMake is configured (tests + examples enabled) and the library,
  tests and examples are compiled

### Requirement: Unit / integration test split

The runner SHALL separate pure-C++ unit tests from SymPy-oracle integration
tests using the `[oracle]` Catch2 tag: `just test-unit` SHALL run `~[oracle]`
(no Python needed) and `just test-integration` SHALL run `[oracle]`. `just test`
SHALL run the full suite, and `just test-filter <pattern>` SHALL forward an
arbitrary Catch2 tag/name filter.

#### Scenario: Unit tests need no oracle
- **WHEN** `just test-unit` is run
- **THEN** only tests without the `[oracle]` tag execute and they pass without a
  Python/SymPy subprocess

#### Scenario: Integration tests use the oracle
- **WHEN** `just test-integration` is run
- **THEN** the `[oracle]` tests execute, cross-checked against SymPy, and pass

### Requirement: Examples and parity recipes

The runner SHALL build and run the example programs (`just examples`, and
`just example <name>` for one) and SHALL build and run the Tier 1 + Tier 2
acceptance probe (`just parity`), reporting the probe's `REMAINING` count. The
parity recipe SHALL resolve GMP/MPFR include and library paths from the CMake
cache so it works on any configured host.

#### Scenario: Run all examples
- **WHEN** `just examples` is run
- **THEN** every built example executable is run in turn

#### Scenario: Parity gate
- **WHEN** `just parity` is run on a green tree
- **THEN** the acceptance probe is compiled and run and reports `REMAINING=0`
