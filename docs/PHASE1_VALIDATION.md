# Phase 1 Validation: Pipeline Calibration Against Known Mutations

**Date:** 2026-06-21 (Day 22, recall) / 2026-06-22 (Day 23, precision)
**Outcome:** 2 (mixed). No code bug. Pipeline is biased-conservative: poor recall, strong precision on confident negatives. Proceed to Phase 2 with documented caveats.
**Repo commit:** `<FILL IN>`

---

## 1. Purpose

Before committing Phase 2 compute (MD, docking), validate the corrected FoldX + ThermoMPNN consensus pipeline against mutations whose effect is already known — in **both** directions:

- **Recall** — given known *stabilizers*, does the pipeline find them?
- **Precision** — given known *destabilizers*, does the pipeline reject them?

Prior to this, the FoldX side had exactly **one** positive control (N233K); a residual bug (RepairPDB, BuildModel parameters, `.fxout` parsing) would have been invisible.

## 2. Method

- **Structure:** `6eqe_Repair.pdb` — IsPETase WT (PDB 6EQE is the Austin 2018 wild-type structure). Single chain A. Catalytic triad S160 / D206 / H237.
- **Numbering:** PDB position = ThermoMPNN position + 29. Confirmed by FoldX (`ASNA233` -> internal index 204) and by the N233K ThermoMPNN lookup.
- **FoldX:** BuildModel, **1 run** (`numberOfRuns` default); parse `Dif_6eqe_Repair.fxout`. Negative total energy = stabilizing.
- **ThermoMPNN:** site-saturation inference (`ThermoMPNN_inference_6eqe_Repair.csv`). Negative `ddG_pred` = stabilizing (sign convention corrected in Week 1).
- **Reproduction check:** N233K re-run reproduced the stored value **bit-for-bit** (-0.685467 kcal/mol), confirming the script is deterministic and faithful before testing unknowns.
- **Positive controls (recall) — two buckets:**
  - **Stability bucket (verdict):** S121E, D186H, R280A (ThermoPETase) and R224Q (FAST-PETase), plus N233K (FAST-PETase, already validated). Folding-stability mutations.
  - **Cleft bucket (no verdict):** W159H, S238F (Austin 2018) — active-site-narrowing *activity* mutations, not folding stabilizers. Included only to illustrate that the tools are blind to activity-driven effects.
- **Negative controls (precision):** six deeply buried (RSA 0.00, core, not near active site) mutations selected as confident destabilizers — three **cavity-creating** (W257A, F106A, L101A) and three **charge-burial** (V68D, I83R, I145D).

## 3. Results

### 3a. Recall — known stabilizers (negative = stabilizing)

| Mutation | FoldX ddG | ThermoMPNN | FoldX | ThermoMPNN | Consensus |
|---|---|---|:--:|:--:|:--:|
| N233K | -0.685 | -1.198 | YES | YES | PASS |
| S121E | -0.444 | +0.022 | YES | no | no |
| D186H | +0.526 | -0.512 | no | YES | no |
| R224Q | +0.097 | +0.378 | no | no | no |
| R280A | +0.559 | -0.069 (approx 0) | no | no | no |

**Recall on known stabilizers: FoldX 2/5 - ThermoMPNN 2/5 - Consensus 1/5 (N233K only).**

The two methods are **anti-correlated on the misses**: where FoldX recovers a stabilizer (S121E) ThermoMPNN misses it, and vice versa (D186H). The only mutation both recover is N233K — the single control with direct wet-lab validation.

**Cleft bucket (no verdict):** W159H (FoldX +3.743, ThermoMPNN +1.250) and S238F (FoldX -0.532, ThermoMPNN +0.054). Both behave as folding-non-stabilizers — expected, since these are activity mutations the tools cannot register.

### 3b. Precision — known destabilizers (positive = destabilizing = correctly rejected)

| Mutation | Type | FoldX ddG | ThermoMPNN | Verdict |
|---|---|---|---|:--:|
| W257A | cavity | +5.13 | +3.75 | rejected |
| F106A | cavity | +3.34 | +2.68 | rejected |
| L101A | cavity | +3.79 | +2.34 | rejected |
| V68D | charge burial | +4.45 | +3.26 | rejected |
| I83R | charge burial | +11.81 | +2.79 | rejected |
| I145D | charge burial | +4.26 | +2.91 | rejected |

**False-positive rate: FoldX 0/6 - ThermoMPNN 0/6 - Consensus 0/6.**

Both tools correctly flagged all six as destabilizing. FoldX penalized the bulky charge-burial I83R most heavily (+11.8). Unlike recall, the two methods **agree completely** on the negatives.

## 4. Interpretation

- **No code bug.** N233K reproduces bit-for-bit; both N233K and S121E are recovered in the correct direction; all six negatives are correctly rejected. A systematic error would corrupt all cases. **Outcome 2 (mixed), not Outcome 3 (debug).**
- **The pipeline is biased-conservative — and recall and precision are two faces of the same bias.** FoldX leans toward predicting destabilization. That single tendency *lowers recall* (it under-calls real stabilizers, 2/5) and *protects precision* (a tool reluctant to say "stabilizing" rarely says it about a bad mutation, 0/6 false positives). Bad recall and good precision are the same coin.
- **Why the methods disagree on recall but agree on precision.** Strong destabilization (cavities, buried charges) produces large, unambiguous energy penalties that both a physics model and an ML model catch every time. Real stabilization is small, distributed, and subtle — genuinely hard, and where the two methods diverge. The tools are **reliable at rejecting clearly-bad mutations, unreliable at finding good ones.**
- **Consensus is a strict, precision-oriented filter.** It rejected 4/5 known stabilizers (high false-negative rate) *and* 6/6 confident destabilizers (zero false positives). For a candidate-selection filter this is the favorable error profile: it misses good mutations but does not wave through bad ones, so a mutation that *passes* consensus carries real (if not absolute) weight.

## 5. Limitations

- **Negatives are *easy* cases.** Deeply buried cavity/charge-burial mutations carry large, obvious penalties. 0/6 here is necessary but **not sufficient** — it does not establish precision on *subtle* or *neutral* mutations, which remain untested. Failing even these would have been damning; passing them is reassuring but bounded.
- **Precision is about direction, not magnitude.** Correctly signing destabilizers says nothing about the accuracy of a stabilizing *value*. The FoldX ddG magnitude (e.g. -6.57) remains unreliable.
- **Recall controls are combinatorial-design singles.** ThermoPETase/FAST-PETase mutations were characterized mainly within multi-mutation combinations; no published *per-mutation* experimental ddG is in hand, so the magnitude of each miss cannot be quantified and some individual effects may genuinely be small.
- **FoldX at 1 run — no error bars.** Consistent with how the original singles were generated, but conflicts with the handoff's "numberOfRuns=3" convention. A 3-run rerun is pending.

## 6. Implications for the lead candidate

The 5-mutant lead (N233K / T77I / A179V / R260F / Q119D) was selected substantially on FoldX/ThermoMPNN consensus. Given both validation halves:

- **The headline FoldX ddG (-6.57 kcal/mol) is not a claimable result.** FoldX mis-signs 3 of 5 known stabilizers on this protein; its ddG magnitude on novel mutations is unreliable.
- **But "passed consensus" now carries weight it did not two days ago.** The filter demonstrably rejects obvious destabilizers (0/6) and is strict (rejects 4/5 real stabilizers). A strict filter that does not accept bad mutations makes its rare positive verdicts — including the 5 candidates — meaningfully more credible. Confidence moved from *uninformative* to *weak-but-real positive evidence*.
- **The lead is a candidate, not a validated design.** Honest description: *"Generated by FoldX/ThermoMPNN consensus, anchored on the experimentally-validated N233K mutation; the consensus filter rejects confident destabilizers (0/6) but its discrimination on subtle cases is untested. Pending MD and experimental validation."*
- **Confidence rests on:** (1) the N233K experimental anchor, (2) consensus as a strict filter shown to reject easy negatives, and (3) forthcoming Phase 2 MD — **not** on the FoldX ddG magnitude.

## 7. Decision and next steps

**Decision:** Outcome 2. Proceed to Phase 2 **with documented caveats** — not full-confidence production, not pipeline debugging.

**Next steps (priority order):**
1. **3-run FoldX rerun** of the recall + precision sets for error bars (defensibility for the writeup / Furst).
2. **Subtle/neutral negatives** — extend the precision test beyond easy buried cases (surface mutations, conservative substitutions) to probe fine-grained discrimination.
3. **Stronger positive controls** — identify >=1 mutation with published *individual* experimental ddG to anchor future calibration.

---

*Working draft. Verify every value against source files and edit into the author's own voice before committing.*