# Austin et al. 2018 — IsPETase crystal structure and engineering

## Citation
Austin, H.P. et al. (2018). Characterization and engineering of a plastic-degrading aromatic polyesterase. *PNAS*, 115(19), E4350-E4357.
DOI: 10.1073/pnas.1718804115

## What was discovered
This paper solved the first crystal structure of IsPETase (PDB 6EQE, 0.92 Å resolution) — the enzyme Yoshida 2016 had discovered. The structure revealed three things that distinguish IsPETase from other cutinases: a wider active site cleft, an unusually polarized surface charge distribution, and two disulfide bonds providing rigidity near the active site. They also performed the first engineering attempts on IsPETase, including the S238F/W159H double mutant which surprisingly improved activity.

## Why it matters for this project
- This paper produced PDB 6EQE — the structure I parsed in notebook 02 (catalytic triad: Ser-160, Asp-206, His-237)
- It identifies which structural features make IsPETase PET-specific — these features must be preserved when engineering thermostable variants
- The disulfide bonds (Cys-203/Cys-239 and Cys-273/Cys-289) contribute to the rigidity of the active site region; they should NOT be mutated
- It set the precedent that even subtle mutations near the active site can change activity — a warning that random mutations are risky

## Key facts to remember
- **PDB ID solved:** 6EQE (also 5XJH from another group around the same time)
- **Resolution:** 0.92 Å (very high — every atom positioned precisely)
- **Catalytic triad confirmed:** Ser-160, Asp-206, His-237 (matches Yoshida)
- **Disulfide bonds:** two of them, located near the active site (Cys-203/Cys-239 and Cys-273/Cys-289)
- **Structurally unique vs. other cutinases:**
  - Wider, more open active site cleft (allows rigid PET chain to fit; narrower cleft of LCC/TfH cannot accommodate PET)
  - Trp-159 sits at the cleft edge, extending the hydrophobic surface adjacent to the catalytic center
  - Highly polarized surface charge distribution (one face is more charged than the other)
- **Mutations they tested:** S238F/W159H double mutant; W185A single mutant
- **Most successful mutation:** S238F/W159H — modestly improved PET-degrading activity
- **Mechanism of improvement:** primarily better packing near the active site, with a small secondary thermostability gain (~2-3°C Tm increase)
- **What S238F/W159H does structurally:** narrows the cleft slightly to be more like other cutinases — counterintuitively, this helped activity at moderate temperatures

## Open questions / follow-ups
- Austin only tested a handful of mutations by intuition. With ThermoMPNN + AlphaFold + FoldX, the scale of computational exploration is now thousands of mutations. What did Austin miss?
- The wide active site cleft is what makes IsPETase PET-specific. Are there mutations that raise thermostability without narrowing the cleft? (This is the core question for my project)
- Lu 2022 (FAST-PETase) used machine learning to find thermostable variants and got +12°C Tm. Read that paper next.

## My takeaway
IsPETase's PET-specificity comes from three structural features: a wide active site cleft, hydrophobic residues like Trp-159 lining it, and two disulfide bonds maintaining its geometry. Engineering thermostable variants requires preserving all three. Austin's 2018 attempt made 2 mutations by visual inspection of the structure; I'll use modern computational tools to screen thousands while respecting these constraints.