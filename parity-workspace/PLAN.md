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

### Batch 5 (final) — DONE
Closed all 3: trivariate Wang (linear-form roots; x^3+y^3+z^3-3xyz),
hyperexpand 2F1(1/2,1;3/2;z^2)->atanh(z)/z, gamma-ratio limit
Gamma(x+1)/Gamma(x+1/2)->oo (engine already correct; added oracle-checked test).
Fixed a probe bug: oracle.equivalent can't certify oo (oo-oo=nan); expr_is now
matches infinite values by string equality first.

ACCEPTANCE: REMAINING=0 / TOTAL=22. Tier 1 + Tier 2 finished on the acceptance
set. Suite: 7480 assertions / 1774 cases, all green. LOOP STOPPED.

## Deep-completion pass (#4 Wang, hyperexpand, Gruntz, group theory)
Expanded probe to 38 cases. Baseline REMAINING=2:
- Wang (x+y)^2-z^2 [degree-2 multivariate, >1 other var]
- hyperexpand 0F1(;1/2;z^2/4)->cosh(z)
Gruntz towers all pass; hard bivariate quartic x^4+x^2y^2+y^4 already factors.
Group theory: needs conjugacy_classes/center/derived_subgroup API + probe section
(guarded by PROBE_GROUP).

### Deep-completion pass — DONE
Batch 6 closed all remaining: Wang degree-2 multivariate ((x+y)^2-z^2 via
linear-form deflation), hyperexpand 0F1(;1/2;z^2/4)->cosh / 0F1(;3/2;z^2/4)->
sinh(z)/z, and group theory (conjugacy_classes / num_conjugacy_classes /
group_center / derived_subgroup, validated vs S3/S4/A4/D4). Probe extended with
a PROBE_GROUP section.

ACCEPTANCE (full set incl group, build with -DPROBE_GROUP): REMAINING=0 / TOTAL=48.
Suite: 7503 assertions / 1776 cases, all green. #4 Wang, hyperexpand, Gruntz,
group theory finished on the acceptance set. LOOP STOPPED.

Probe rebuild+run:
  g++ -std=c++20 -I include -I tests -I /tmp/deps/root/usr/include \
    -I build/_deps/nlohmann_json-src/include -DPROBE_GROUP \
    -DSYMPP_TEST_ORACLE_SCRIPT="\"$PWD/oracle/oracle.py\"" \
    -DSYMPP_TEST_PYTHON_EXECUTABLE="\"python3\"" \
    parity-workspace/acceptance_probe.cpp tests/oracle/oracle.cpp -o parity-workspace/probe \
    build/src/libsympp.so -L<gmp-lib> -Wl,-rpath,<gmp-lib> -lgmpxx -lgmp -lmpfr
  LD_LIBRARY_PATH=build/src:<gmp-lib> parity-workspace/probe

## Tier 2 full-completion pass
Deepened probe to 70 cases (14 Risch integrals, Wang quartic-4-factor + biquadratic,
hyperexpand degenerate 2F1(2,1;1;z), 4 adversarial Gruntz towers). Baseline:
REMAINING=0/70 — all already pass. Remaining Tier 2 work needs NEW API:
group is_solvable/is_nilpotent, physics multi-mode operators.

### Tier 2 full-completion — DONE
Added group is_solvable/is_nilpotent (S3/S4/A4/A5/D4/C6 vs SymPy) and multi-mode
bosonic Fock operators ([a_i,a_j†]=δ_ij). Probe extended with PROBE_GROUP +
PROBE_PHYS sections.

ACCEPTANCE (full Tier 1+2, build -DPROBE_GROUP -DPROBE_PHYS): REMAINING=0 / TOTAL=79.
Suite: 7543 assertions / 1785 cases, all green. ALL OF TIER 2 finished on the
79-case acceptance set (#4 Wang, #5 Risch, #6 hyperexpand, #7 Gruntz,
#8 group theory, #9 physics). LOOP STOPPED.

## Residual-depth pass
Probe deepened to 87 cases. Baseline REMAINING=3 (all Wang 4-variable:
ab-ac+bd-cd, wx+wy+zx+zy, x^2y^2-z^2w^2). Residual integrals (x^2 e^-x^2,
1/(x log^2 x), cos/x->Ci, 0..oo x^2 e^-x^2=sqrt(pi)/4) and polylog 3F2 already pass.
New-API targets: factor(x^N-1) cyclotomic perf (x^30-1 was >30s), group
is_simple/abelian_invariants, physics normal-ordering.
BLOCKED (documented, beyond a loop / need new infra): full Risch non-elementarity
proofs (SymPy doesn't do them either), elliptic-K hyperexpand (needs elliptic
functions), arbitrary-arity Wang general Hensel (multi-week).
