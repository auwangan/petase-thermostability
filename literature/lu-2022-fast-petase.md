Lu et al. 2022 — FAST-PETase (machine learning-aided PETase engineering)
Citation
Lu, H. et al. (2022). Machine learning-aided engineering of hydrolases for PET depolymerization. Nature, 604(7907), 662-667.
DOI: 10.1038/s41586-022-04599-z
What was discovered
This paper used a structure-based, self-supervised 3D convolutional neural network (MutCompute) to predict stabilizing/activity-improving mutations on IsPETase, producing an engineered variant called FAST-PETase (Functional, Active, Stable, Tolerant). The model was trained on ~19,000 protein structures of varied function. From an initial pool of predicted variants, the authors narrowed to a final design carrying five mutations, which fully degraded 51 different post-consumer PET products at 50C within a week (some within 24h), and enabled closed-loop recycling (degraded monomers re-polymerised into virgin PET).
Why it matters for this project

Closest methodological analog to my pipeline: ML-predicted mutations on IsPETase, narrowed to a small validated set. This is my blueprint.
It demonstrates the target outcome — a small number of mutations (5) producing a robust, more thermostable, more active enzyme.
CRITICAL CONTRAST with my current candidate list (see below): FAST-PETase's mutations are surface, charge/polar rearrangements — NOT buried-core prolines or charges. My ThermoMPNN top-30 was dominated by exactly the mutation types FAST-PETase avoided, suggesting my top-30 is mostly false positives that FoldX should reject.

Key facts to remember

Model used: MutCompute — 3D self-supervised CNN; predicts which native residues are a poor structural fit in their local environment, flags better-fitting substitutions.
Training set: ~19,000 protein structures (diverse functions).
Starting scaffolds: wild-type IsPETase plus two prior engineered variants (ThermoPETase, DuraPETase) — i.e. they did NOT start only from wild-type.
Final FAST-PETase = 5 mutations vs wild-type:

From ML prediction: N233K, R224Q, S121E
From scaffold: D186H, R280A


Location of the five mutations: none are mutations TO proline (confirmed: N->K, R->Q, S->E, D->H, R->A). Burial state (core vs surface) of each NOT yet verified per-residue — [VERIFY exact RSA of each from the structure / SI before citing "all surface" as a hard claim].
Active site / R280A caveat: R280A is NOT a pure stability mutation. Companion paper states it "alters the substrate binding site conformation to allow stronger substrate interactions" — i.e. partly an activity/substrate-binding effect. So the clean story "all five are surface stability mutations" is too simple; R280A touches the binding site. [Confirmed from companion paper.]
Tm / thermostability gain: [STILL UNVERIFIED — do NOT state a number. I loosely recalled "+12C" from a summary; this may be a confusion with the Furst CHMO paper's +13C. Confirm from Lu 2022 directly before using.]
Activity result: up to ~38-fold higher activity at 40-50C vs prior PETases (from secondary coverage; verify exact figure/conditions in the paper itself — UNVERIFIED).

Mechanism (from open-access MD companion paper) — VERIFIED
Companion: "Computational analysis reveals temperature-induced stabilization of FAST-PETase", Comput Struct Biotechnol J 2025 (PMC11946493 / doi:10.1016/j.csbj.2025.03.006, open access).

Comparative Constraint Network Analysis (CNAnalysis) + MD of WT vs FAST-PETase at 30C and 50C. [Confirmed]
Identified thermolabile (flexible) sequence stretches in the wild-type enzyme. [Confirmed]
All five FAST-PETase mutations affect these thermolabile/critical regions — stabilization works by rigidifying flexible regions. [Confirmed, direct quote]
Specific examples: rigidity increased in region 110-155 (around position 121 -> S121E) and region 200-230 (around position 224 -> R224Q). Region 55-60 actually weakened. [Confirmed]
R280A region also gained rigidity AND alters substrate binding site conformation for stronger substrate interactions. [Confirmed]
Counterintuitive: flexibility of these regions DECREASED at 50C in FAST-PETase (opposite of usual thermodynamic expectation). [Confirmed]

Implications for my methodology

Distrust raw ThermoMPNN ranking. My top-30 (core prolines, buried charges) contradicts the gold-standard result. Use FoldX as the physics-based cross-check; require agreement.
Consider chemical-sanity filters in filter_candidates.py: down-weight/flag mutations TO proline and charged residues into low-RSA (buried) positions. Justified by the FAST-PETase pattern.
Site selection for Phase 2: use MD-identified flexible regions as the mutation-site selector, not just raw ddG ranking.
Combinations may be non-additive (epistasis): likely cannot assume top-5 individual mutations sum to +5C; should test combinations explicitly. [UNVERIFIED — confirm the specific epistasis claim against Lu 2022 before citing.]

Open questions / follow-ups

How exactly did they rank/narrow from the full prediction pool to the final set? (Free Supplementary Methods on the Nature page — read directly, do not rely on summary.)
Does the scaffold-first strategy (building on ThermoPETase/DuraPETase) apply to me, or do I deliberately start from wild-type to keep the story clean?
Email Furst (already agreed to give feedback) the specific question of whether requiring FoldX/ThermoMPNN agreement is a sound way to suppress the proline false positives.

My takeaway
FAST-PETase reached a robust, more thermostable enzyme with five mutations that all target thermolabile/flexible regions of the protein (confirmed by the MD companion paper), and none are core prolines. My ThermoMPNN top-30 looks the opposite (core prolines, buried charges), which is suggestive evidence those high-ranked predictions are mostly artifacts. The paper validates my two-tool design (predict, then physics-check) and points Phase 2 toward flexibility-guided site selection rather than blind ddG ranking. Note: R280A is partly a substrate-binding/activity mutation, not pure stability — the "all five are clean surface stability mutations" framing is too simple.