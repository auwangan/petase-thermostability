#Step 1 — Imports and parse
#Step 2 — Define your reference data (constants at top of file
#Step 3 — Compute SASA for the whole structure
#Step 4 — Identify the triad atoms
#Step 5 — Main loop: one pass over chain A residues
from pandas import DataFrame
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley   
import numpy as np 
parser = PDBParser(QUIET=True)
structure = parser.get_structure("6QED", 'data/6eqe_Repair.pdb')
model = structure[0]
chain_a = model['A']
triad_ids = [160, 206, 237]
sasa_calculator = ShrakeRupley()
sasa_calculator.compute(structure, level="R")
triad_atoms =[]
rows = []
for residue in chain_a:
    if residue.get_id()[1] in triad_ids:
        triad_atoms.extend(residue.get_atoms())
MAX_SASA = {
    'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0,
    'GLU': 223.0, 'GLN': 225.0, 'GLY': 104.0, 'HIS': 224.0, 'ILE': 197.0,
    'LEU': 201.0, 'LYS': 236.0, 'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0,
    'SER': 155.0, 'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0,
}
KD_HYDROPHOBICITY = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
    'GLU': -3.5, 'GLN': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
    'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
    'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2,
}
for residue in chain_a:
    position = residue.get_id()[1]
    resname = residue.get_resname()

    # distance (the nested loop from Step 5)
    dists = []
    for atom in residue.get_atoms():
        for triad_atom in triad_atoms:
            dists.append(atom - triad_atom)
    min_dist = min(dists)

    # rsa
    rsa = residue.sasa / MAX_SASA[resname]
    location = "core" if rsa < 0.20 else "surface"

    # hydrophobicity
    hydrophobicity = KD_HYDROPHOBICITY[resname]

    # active site flag
    near_active_site = "YES" if min_dist < 8.0 else "NO"

    rows.append({
    "position": position,
    "resname": resname,
    "min_dist_to_triad": round(min_dist, 2),
    "rsa": round(rsa, 2),
    "location": location,
    "hydrophobicity": hydrophobicity,
    "near_active_site": near_active_site,
    })
df = DataFrame(rows)
df.to_csv("data/residue_features.csv", index=False)
print(df)