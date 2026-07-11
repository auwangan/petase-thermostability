from openmm.app import PDBFile, ForceField, Modeller
from pdbfixer import PDBFixer
import openmm.app as app
 
IN = "data/md/mut_v2_ss.pdb"
 
fixer = PDBFixer(filename=IN)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)   # protonate at pH 7; bonded Cys should get NO HG
 
top = fixer.topology
 
# 1) list all disulfide (SG-SG) bonds OpenMM sees
ss = []
for bond in top.bonds():
    a1, a2 = bond[0], bond[1]
    if a1.name == "SG" and a2.name == "SG":
        ss.append((a1.residue.id, a2.residue.id))
print("SG-SG bonds in topology:", ss)
 
want = {"79", "153"}
found = any(set(p) == want for p in ss)
print(f"79-153 disulfide present: {'YES' if found else 'NO  <-- PROBLEM'}")
 
# 2) confirm the two cysteines are deprotonated (no HG hydrogen)
for res in top.residues():
    if res.id in ("79", "153") and res.name in ("CYS", "CYX"):
        hs = [a.name for a in res.atoms() if a.name in ("HG", "HG1")]
        print(f"  resi {res.id} {res.name}: thiol H = {hs if hs else 'none (bonded, correct)'}")
 
# 3) also confirm native disulfides still intact (sanity)
natives = [p for p in ss if set(p) in ({"203","239"}, {"273","289"})]
print(f"native disulfides intact: {len(natives)}/2")