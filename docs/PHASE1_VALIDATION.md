# Phase 1 Validation: Pipeline Calibration Against Known Stabilizers

**Date:** 2026-06-21 (Week 3, Day 22)
**Outcome:** 2 (mixed). No code bug; tool reliability limited on IsPETase. Proceed to Phase 2 with documented caveats.

---

## 1. Purpose

Before committing Phase 2 compute (MD, docking), validate the corrected FoldX + ThermoMPNN consensus pipeline against mutations whose stabilizing effect is already known from the literature.

Prior to this step the FoldX side of the pipeline had exactly **one** positive control (N233K). A residual bug in FoldX-side code (RepairPDB use, BuildModel parameters, `.fxout` parsing) would have been invisible. This step adds calibration points and asks a direct question: **does the pipeline recover known stabilizers?**

## 2. Method

- **Structure:** `6eqe_Repair.pdb` — IsPETase WT (PDB 6EQE is the Austin 2018 wild-type structure). Single chain A. Catalytic triad S160 / D206 / H237.
- **Numbering:** PDB position = ThermoMPNN position + 29. Confirmed by FoldX (`ASNA233` → internal index 204) and by the N233K ThermoMPNN lookup.
- **FoldX:** BuildModel, **1 run** (`numberOfRuns` default); parse `Dif_6eqe_Repair.fxout`. Negative total energy = stabilizing.
- **ThermoMPNN:** site-saturation inference (`ThermoMPNN_inference_6eqe_Repair.csv`). Negative `ddG_pred` = stabilizing (sign convention corrected in Week 1).
- **Reproduction check:** N233K re-run through the current script reproduced the stored value **bit-for-bit** (−0.685467 kcal/mol), confirming the script is deterministic and faithful before testing unknowns.
- **Controls — two buckets:**
  - **Stability bucket (verdict):** S121E, D186H, R280A (ThermoPETase) and R224Q (FAST-PETase), plus N233K (FAST-PETase, already validated). All are folding-stability mutations.
  - **Cleft bucket (no verdict):** W159H, S238F (Austin 2018) — active-site-narrowing *activity* mutations, not folding stabilizers. FoldX/ThermoMPNN measure folding ΔΔG and are not expected to score these as stabilizing; included only to illustrate that the tools are blind to activity-driven effects.

## 3. Results

### Stability bucket (negative = stabilizing)

| Mutation | FoldX ΔΔG | ThermoMPNN | FoldX | ThermoMPNN | Consensus |
|---|---|---|:--:|:--:|:--:|
| N233K | −0.685 | −1.198 | ✅ | ✅ | ✅ pass |
| S121E | −0.444 | +0.022 | ✅ | ✗ | ✗ |
| D186H | +0.526 | −0.512 | ✗ | ✅ | ✗ |
| R224Q | +0.097 | +0.378 | ✗ | ✗ | ✗ |
| R280A | +0.559 | −0.069 (≈0) | ✗ | ✗ | ✗ |

**Recall on known stabilizers: FoldX 2/5 · ThermoMPNN 2/5 · Consensus 1/5 (N233K only).**

The two methods are **anti-correlated on the misses**: where FoldX recovers a stabilizer (S121E) ThermoMPNN misses it, and vice versa (D186H). The only mutation both recover is N233K — the single control with direct wet-lab validation.

### Cleft bucket (no verdict)

| Mutation | FoldX ΔΔG | ThermoMPNN | Note |
|---|---|---|---|
| W159H | +3.743 | +1.250 | Both score folding-destabilizing — expected; activity mutation, folding-costly (Trp→His). |
| S238F | −0.532 | +0.054 | FoldX-favorable (added aromatic packing); ThermoMPNN ≈ neutral. Not a folding-stability signal. |

Behaves as predicted: the tools do not register the activity benefit these mutations were designed for.

## 4. Interpretation

- **No code bug.** N233K reproduces bit-for-bit, and both N233K and S121E are recovered in the correct direction. A systematic error (sign flip, wrong WT reference) would corrupt *all* cases, not 2 of 5. **This is Outcome 2 (mixed), not Outcome 3 (debug).**
- **The tools are unreliable on IsPETase folding stability.** Each recovers only ~2/5 known stabilizers, and they disagree with each other. FoldX mis-signs 3 of 5 cases. The misses concentrate on **charge-altering surface mutations**, including two arginine removals (R224Q, R280A) — FoldX's documented weak spot (it over-penalizes loss of charged/salt-bridging residues and undercounts diffuse electrostatic / loop-rigidity stabilization).
- **Recall was measured; precision was not.** Every control was a known stabilizer, so this experiment only tested whether the pipeline *finds* real stabilizers (recall: poor). It did **not** test whether mutations the pipeline *passes* are real (precision). Precision is the property relevant to the design's selected candidates, and it remains unmeasured.
- **Consensus is conservative, not validating.** Requiring both tools to agree would have rejected **4 of 5 known stabilizers** (high false-negative rate). For a design filter this is the safe failure direction — but it means **"passed consensus" is weak positive evidence, not validation.**

## 5. Limitations

- **Controls are combinatorial-design singles.** S121E/D186H/R280A (ThermoPETase) and R224Q/N233K (FAST-PETase) were characterized primarily as parts of multi-mutation combinations, where the large Tm gains arise from synergy. Published *per-mutation* experimental ΔΔG is not in hand, so the magnitude of each miss cannot be quantified, and some individual effects may genuinely be small.
- **FoldX at 1 run — no error bars.** Consistent with how the original single mutants were generated, but this conflicts with the "numberOfRuns=3" convention recorded in the handoff. A 3-run rerun is pending.
- **Precision untested.** No negative controls (known destabilizers/neutrals) were run, so the false-positive rate of the pipeline is unknown.

## 6. Implications for the lead candidate

The 5-mutant lead (N233K / T77I / A179V / R260F / Q119D) was selected substantially on FoldX/ThermoMPNN consensus. Given the above:

- **The headline FoldX ΔΔG (−6.57 kcal/mol) is not a claimable result.** FoldX mis-signs 3 of 5 known cases on this protein; its ΔΔG on novel mutations is correspondingly unreliable as a magnitude or a guarantee.
- **The lead is a candidate, not a validated design.** Honest description: *"Generated by FoldX/ThermoMPNN consensus, anchored on the experimentally-validated N233K mutation, pending MD and experimental validation."*
- **Confidence rests on:** (1) the N233K experimental anchor, (2) consensus as a precision-oriented but *unquantified* filter, and (3) forthcoming Phase 2 MD. It does **not** rest on the FoldX ΔΔG magnitude.

## 7. Decision and next steps

**Decision:** Outcome 2. Proceed to Phase 2 **with documented caveats** — not full-confidence production, not pipeline debugging.

**Next steps (priority order):**
1. **Precision / false-positive test** — run known destabilizing/neutral mutations through the consensus pipeline; measure how many are wrongly passed. This is the experiment that earns or refutes the claim that selected candidates are likely stabilizing.
2. **3-run FoldX rerun** of the controls for error bars (defensibility).
3. **Stronger controls** — identify ≥1 mutation with published *individual* experimental ΔΔG to anchor future calibration.