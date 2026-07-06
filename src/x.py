from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
s = PDBParser(QUIET=True).get_structure("p", "data/6eqe_Repair.pdb")
for pp in PPBuilder().build_peptides(s):
    print(pp.get_sequence())