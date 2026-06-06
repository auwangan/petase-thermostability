```markdown
# Week 0 Reflection — Setup and Foundations

**Dates:** May 27 – June 3, 2026
**Hours invested:** ~24-26 (estimate)
**Status:** Week 0 complete. Ready for Phase 0.

## What I built

- Public GitHub repo with 13 commits at github.com/auwangan/petase-thermostability
- Conda environment with Biopython, pandas, numpy, py3dmol
- Solved Rosalind problems 1, 2, 3, 5, 7 (skipped 4 originally)
- Literature notes on Yoshida 2016 and Austin 2018
- Notebook 01: pandas intro with synthetic thermostability dataset
- Notebook 02: PDB 6EQE parsing + catalytic triad verification (2.94 Å, 2.66 Å H-bonds confirmed)
- Script: src/find_disulfide_candidates.py — recovered both IsPETase disulfide bonds

## What I struggled with most

-I struggled with rosalind the most, they took the most time
-I struggled with prograstinations
-I struggled with the fundementals in cs as well as I have not properly learnt all the jargons
-fibonacci sequence in rosalind was really bad
-i almost skipped yoshida paper cuz i was tired
-most importantly i don't have the basic intuition for the basic use of fundementals python like loops

## Skills gained (1-3 honest scale)

| Skill | Day 1 | Day 7 |
|---|---|---|
| Linux / WSL / CLI | 0 | 1.5 |
| Python | 0 | 2 |
| Git / GitHub | 0 | 2 |
| Pandas | 0 | 1.5 |
| Biopython | 0 | 1.5 |
| PDB / structural biology tools | 0 | 1.5 |
| Literature reading | 1 | 1.5 |
| Project planning / workflow | 1 | 2 |

## What I'm not yet good at

-Writing functions from scratch without copying examples from past work
-comceptual understanding of the PETase not there yet
-understanding conceptual knowledge for coding

## Key conceptual wins

- PET thermostability matters because of semicrystalline → amorphous transition above 70°C
- Catalytic triad geometry verified computationally — these residues must NOT be mutated
- Mutation tradeoffs aren't strictly zero-sum — scaffold mutations can improve stability without losing specificity

## Honest accountability

- I did manage to follow the plan everyday which i did expect myself to do
- but still gets lazy somehow with coding and use copilot for some of the work
- I googled 2 answers for rosalind on day 7
- Used Copilot autocomplete a few times despite the no-AI rule for Rosalind. 

## Plan for Week 1

[Brief — 3-5 bullets]
- FoldX install (Days 8-9)
- ThermoMPNN setup (Day 9-10)
- Read Lu 2022 FAST-PETase paper (Day 13)
- Friend joins Day 11
- End of week: finalize Phase 1 scope