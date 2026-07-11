import pymol
from pymol import cmd
import math
 
IN  = "data/md/mut_v2_renumbered.pdb"      # was the 6eqe FoldX file
OUT = "data/md/mut_v2_ss.pdb"              # new output name
CHAIN = "A"
IDEAL = 2.05
 
def atom_xyz(resi, name):
    m = cmd.get_model(f"prot and chain {CHAIN} and resi {resi} and name {name}").atom
    return m[0].coord if m else None
 
def sg_sg():
    a, b = atom_xyz(79, "SG"), atom_xyz(153, "SG")
    return math.dist(a, b) if (a and b) else None
 
pymol.finish_launching(["pymol", "-qc"])
cmd.load(IN, "prot")
cmd.remove("solvent")
print(f"SG-SG before: {sg_sg():.2f} A")
 
best = None
# scan chi1 of BOTH cysteines over full rotation, 5-degree steps
for a79 in range(0, 360, 5):
    cmd.set_dihedral(f"prot and resi 79 and name N",
                     f"prot and resi 79 and name CA",
                     f"prot and resi 79 and name CB",
                     f"prot and resi 79 and name SG", a79)
    for a153 in range(0, 360, 5):
        cmd.set_dihedral(f"prot and resi 153 and name N",
                         f"prot and resi 153 and name CA",
                         f"prot and resi 153 and name CB",
                         f"prot and resi 153 and name SG", a153)
        d = sg_sg()
        if d is not None and (best is None or abs(d - IDEAL) < abs(best[0] - IDEAL)):
            best = (d, a79, a153)
 
d, a79, a153 = best
# re-apply the winning dihedrals
cmd.set_dihedral("prot and resi 79 and name N", "prot and resi 79 and name CA",
                 "prot and resi 79 and name CB", "prot and resi 79 and name SG", a79)
cmd.set_dihedral("prot and resi 153 and name N", "prot and resi 153 and name CA",
                 "prot and resi 153 and name CB", "prot and resi 153 and name SG", a153)
 
print(f"best SG-SG   : {d:.2f} A  (chi1_79={a79}, chi1_153={a153})")
if d <= 2.5:
    print("  OK -> bonding geometry set.")
    cmd.save(OUT, "prot")
    print(f"  wrote {OUT}")
else:
    print("  best reachable is still >2.5 A by chi1 alone; may need chi2 or the "
          "FoldX backbone differs from the crystal we validated on.")