## Why

SymPP's number-theory layer already ships `factorint`, `divisors`, `igcdex`,
`jacobi_symbol`, `continued_fraction`, `n_order`, `primitive_root` and
`sqrt_mod`, but three staples of `sympy.ntheory` / `sympy.solvers.diophantine`
are still missing: the **Chinese Remainder Theorem**, the **discrete
logarithm**, and **linear Diophantine** equation solving. These are the
cheapest remaining genuine feature win toward SymPy parity and underpin
common cryptography / modular-arithmetic workflows the library now supports
(it already has RSA/DH/ElGamal and elliptic curves).

## What Changes

- **`crt(moduli, residues)`** — solve a simultaneous system `x ≡ rᵢ (mod mᵢ)`
  by pairwise combination; works for non-coprime moduli, returns `(x, M)` with
  `M = lcm(mᵢ)` and `0 ≤ x < M`, or `nullopt` when the congruences are
  inconsistent.
- **`discrete_log(n, a, b)`** — least non-negative `x` with `bˣ ≡ a (mod n)`
  via baby-step/giant-step, or `nullopt` when no solution exists.
- **`diop_linear(a, b, c)`** — solve `a·x + b·y = c` over ℤ, returning a
  parametric family `x = x₀ + (b/g)·t`, `y = y₀ − (a/g)·t`, or `nullopt` when
  `gcd(a, b) ∤ c`.

All three are GMP-backed, cross-checked against SymPy through the oracle, and
documented in the README and parity roadmap.

## Capabilities

### New Capabilities
- `number-theory`: Chinese Remainder Theorem, discrete logarithm, and linear
  Diophantine equation solving over arbitrary-precision integers, alongside the
  existing residue/factor primitives.

## Non-goals

- General (non-linear / multivariate) Diophantine solving — Pell, sum-of-squares,
  ternary quadratic forms — remains future work.
- Index-calculus / Pollard-rho discrete log for cryptographically large moduli;
  baby-step/giant-step is `O(√n)` and aimed at teaching-scale inputs.
