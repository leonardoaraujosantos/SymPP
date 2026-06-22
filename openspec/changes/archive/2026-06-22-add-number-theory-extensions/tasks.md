# Tasks

Single zero-regression increment on `main`, verified against the SymPy oracle
with regression tests and a green full suite.

## 1. Implementation (`functions/ntheory.{hpp,cpp}`)

- [x] 1.1 `crt(moduli, residues)` — pairwise CRT combination over `lcm`, with
  consistency check and `nullopt` on contradiction (CRT)
- [x] 1.2 `discrete_log(n, a, b)` — baby-step/giant-step, `nullopt` when unsolvable
  (DLOG)
- [x] 1.3 `diop_linear(a, b, c)` — extended-Euclid base solution + parametric
  family, `nullopt` when `gcd ∤ c` (DIOP-LINEAR)

## 2. Oracle

- [x] 2.1 Extend the `ntheory` oracle op with `crt` and `discrete_log`

## 3. Tests (`tests/functions/ntheory_test.cpp`)

- [x] 3.1 CRT: coprime, non-coprime-consistent, inconsistent (oracle-checked)
- [x] 3.2 discrete_log: solvable cases vs SymPy, no-solution case
- [x] 3.3 diop_linear: solution satisfies the equation for several `t`; no-solution
  case returns `nullopt`

## 4. Documentation

- [x] 4.1 README "What's in the box" number-theory bullet
- [x] 4.2 Parity roadmap: mark number-theory extensions shipped
- [x] 4.3 OpenSpec change archived into `specs/number-theory`
