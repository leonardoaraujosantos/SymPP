# Tier 1 + Tier 2 completion workspace

Goal: drive SymPP to SymPy parity on Tier 1 and Tier 2, looping until the
acceptance probe (`acceptance_probe.cpp`) reports REMAINING=0.

"Finished" = SymPy-parity on the concrete cases below (each oracle- or
identity-checked). Where SymPy itself returns a special function / unevaluated,
parity means matching that — not literal mathematical completeness.

Loop protocol:
- Each batch (workflow) targets the currently-FAILing acceptance cases.
- Every feature lands on main with real regression tests, green-gated.
- After each batch: rebuild + run the probe; record REMAINING here.
- No-progress guard: if an item's failing count is unchanged for 2 consecutive
  batches, mark it BLOCKED (with evidence) and stop looping on it.
- STOP when REMAINING=0 over the achievable set.

## Items & acceptance criteria

| Item | Criteria (see probe for exact cases) | Status |
|------|--------------------------------------|--------|
| T1#1 erf | definite erf/erfc ∫₀^∞ families | tracking |
| T2#4 Wang | bivariate non-monomial roots, bilinear, symmetric 3-var | tracking |
| T2#5 Risch | special-function antiderivatives (Ei/li/Si/erfi) + rational·exp | tracking |
| T2#6 hyperexpand | pFq → elementary (incl 1F1(1;2;z), atanh form) | tracking |
| T2#7 Gruntz | gamma-ratio + tower limit values | tracking |
| T2#8 group | weighted Burnside / bracelet counts | unit-tested |
| T2#9 physics | multi-mode CCR/CAR relations | unit-tested |

## Batch log
(updated each iteration)

### Batch 0 (baseline)
REMAINING=3 / TOTAL=22. Open: Wang trivariate (x^3+y^3+z^3-3xyz),
hyperexpand 2F1(1/2,1;3/2;z^2)->atanh(z)/z, Gruntz Gamma(x+1)/Gamma(x+1/2)->oo.
