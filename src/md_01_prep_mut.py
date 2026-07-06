from pdbfixer import PDBFixer
from openmm.app import PDBFile

inp = "data/alphafold/ispetase_5mut_cfbb9_0/ispetase_5mut_cfbb9_0_relaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb"
out = "data/md/mut_fixed.pdb"

fixer = PDBFixer(filename=inp)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()          # fill any incomplete sidechains
fixer.addMissingHydrogens(7.0)   # add H atoms at pH 7.0

import os
os.makedirs("data/md", exist_ok=True)
with open(out, "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

n_atoms = fixer.topology.getNumAtoms()
print(f"Prepared structure written to {out}")
print(f"Atom count after prep: {n_atoms}")