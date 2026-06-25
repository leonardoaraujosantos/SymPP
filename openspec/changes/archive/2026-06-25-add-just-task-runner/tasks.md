## 1. Task runner

- [x] 1.1 Add `justfile` with build/test/examples/parity/bench/docs/ci recipes.
- [x] 1.2 Split tests via the `[oracle]` tag: `test-unit` = `~[oracle]`,
      `test-integration` = `[oracle]`.
- [x] 1.3 Make `parity` portable by reading GMP/MPFR paths from CMakeCache.
- [x] 1.4 Document the runner in `README.md`.
- [x] 1.5 Verify every recipe runs (build, test-unit, test-integration,
      test-filter, examples, example, parity → REMAINING=0).
