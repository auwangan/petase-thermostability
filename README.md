# PETase Thermostability — Computational Engineering

Computational engineering of thermostable PETase variants for whole-cell
biocatalysis in *E. coli*. Dry-lab only, summer 2026.

## Project goal

Predict 5–10 PETase mutants with ΔTm ≥ +5°C vs IsPETase wild-type, validated
through a combined pipeline of ThermoMPNN, FoldX, AlphaFold confidence
scoring, and short MD simulations. A Phase 3 COBRApy/FBA model evaluates
whole-cell production feasibility in engineered *E. coli*.

## Pipeline

1. **Mutation generation** — ThermoMPNN + saturation scan on stability-relevant residues
2. **Energetic validation** — FoldX ΔΔG predictions
3. **Structural validation** — AlphaFold2 / AlphaFold3 confidence scoring (pLDDT, PAE)
4. **Whole-cell context** — COBRApy FBA on iML1515 with engineered enzyme

## Repository structure

- data/ — input PDB files, sequences, datasets (mostly .gitignored)
- notebooks/ — exploratory Jupyter notebooks
- src/ — reusable Python scripts and modules
- figures/ — generated plots and structural images
- literature/ — paper summaries and references
- learning/ — practice code and tutorials
- rosalind/ — Rosalind problem solutions (Python practice, no AI)
- docs/ — environment recipes, project notes

## Environment setup

conda env create -f docs/environment.yml
conda activate petase

## Author

Auwangan, Year 12 student, Athens, Greece. Summer 2026 independent project.
