from Bio.PDB import PDBParser, Superimposer
import numpy as np

parser = PDBParser(QUIET=True)

af_path = "data/alphafold/ispetase_wt_val_924d7_2/ispetase_wt_val_924d7_2_relaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb"
xtal_path = "data/6eqe_Repair.pdb"

af = parser.get_structure("af", af_path)[0]
xtal = parser.get_structure("xtal", xtal_path)[0]

# AlphaFold numbers 1..265; crystal is PDB-numbered (~29..). Offset = +28.
OFFSET = 28

def ca_by_resnum(chain):
    d = {}
    for res in chain:
        if "CA" in res:
            d[res.id[1]] = res["CA"]
    return d

af_ca   = ca_by_resnum(next(af.get_chains()))
xtal_ca = ca_by_resnum(next(xtal.get_chains()))

# match AF residue i to crystal residue i+OFFSET
af_atoms, xtal_atoms = [], []
for i, atom in af_ca.items():
    j = i + OFFSET
    if j in xtal_ca:
        af_atoms.append(atom)
        xtal_atoms.append(xtal_ca[j])

print(f"AF CA: {len(af_ca)} | crystal CA: {len(xtal_ca)} | matched: {len(af_atoms)}")

sup = Superimposer()
sup.set_atoms(xtal_atoms, af_atoms)   # (reference, moving)
print(f"Backbone Cα RMSD: {sup.rms:.3f} Å")