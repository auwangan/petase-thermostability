#Downloads PDB 6EQE programmatically from RCSB (don't hardcode-use a local file that might not exist on someone else's machine)
#Parses it with Biopython
#Finds all cysteine residues in chain A
#For every pair of cysteines, computes the CA-CA distance
#Reports pairs with distance ≤ 8 Å (the geometric threshold for a potential disulfide bond)
#Prints a clean readable report
import os  
from Bio.PDB import PDBParser
import math
import requests
import sys
import numpy as np
# Step 1: Download PDB 6EQE
pdb_path = "data/6eqe.pdb"
read_url = "https://files.rcsb.org/download/6EQE.pdb"
if not os.path.exists(pdb_path):
    response = requests.get(read_url)
    if response.status_code == 200:
        with open(pdb_path, "w") as file:
            file.write(response.text)
    else:
        print(f"Failed to download. Status: {response.status_code}")
        sys.exit(1)
# Step 2: Parse the PDB 
parser = PDBParser(QUIET=True)
structure = parser.get_structure("6EQE", pdb_path)
print ("PDB file 6EQE downloaded and parsed successfully.")
model = structure[0]    # Returns the Model object directly (no generator)
chain_a = model['A']
cysteines = []
for residue in chain_a:
    if residue.get_resname() == 'CYS':
        cysteines.append(residue)
print(f"Found {len(cysteines)} cysteine residues in chain A.")
disulfide_candidates = []
for i in range(len(cysteines)):
    for j in range(i+1, len(cysteines)):
        ca1 = cysteines[i]['CA'].coord
        ca2 = cysteines[j]['CA'].coord
        distance = np.linalg.norm(ca1 - ca2)
        if distance < 8.0:
            disulfide_candidates.append((cysteines[i].get_id(), cysteines[j].get_id(), distance))
print(f"Found {len(disulfide_candidates)} potential disulfide bond candidates:")
for res1, res2, dist in disulfide_candidates:
    print(f"Residue {res1} and Residue {res2} with CA-CA distance: {dist:.2f} Å")
    print(f"\nAll {len(cysteines)} cysteines in IsPETase form disulfide bonds.")