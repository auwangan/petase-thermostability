from Bio.PDB import PDBParser
s = PDBParser(QUIET=True).get_structure("p", "data/6eqe_Repair.pdb")
for chain in s[0]:
    for r in chain:
        if r.id[1] in [121, 159, 186, 224, 238, 280]:
            print(chain.id, r.id[1], r.resname)