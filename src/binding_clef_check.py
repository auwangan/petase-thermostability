# src/check_binding_cleft.py
from Bio.PDB import PDBParser

# IsPETase substrate-binding cleft (6EQE numbering): triad + pocket-lining residues
CLEFT = [87, 159, 160, 161, 185, 206, 208, 237, 238, 241, 280]
candidates = [(260,'R','F'), (119,'Q','D'), (127,'Q','L'), (233,'N','A'),
              (225,'N','C'), (77,'T','I'), (179,'A','V')]

chain = PDBParser(QUIET=True).get_structure("x", "data/6eqe_Repair.pdb")[0]["A"]
cleft_atoms = [a for res in chain if res.get_id()[1] in CLEFT for a in res.get_atoms()]

print(f"{'mutation':9}{'min_dist_to_cleft':>18}   flag")
for pos, wt, mt in candidates:
    d = min(a - c for a in chain[pos].get_atoms() for c in cleft_atoms)
    flag = "CHECK (near cleft)" if d < 8.0 else "ok"
    print(f"{wt}{pos}{mt:<5}{round(d,2):>18}   {flag}")