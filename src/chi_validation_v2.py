import pymol
from pymol import cmd
import itertools
 
PDB = "/home/auwangan/petase-thermostability/data/6eqe_Repair.pdb"  # absolute path
CHAIN = "A"
PAIRS = [(135, 151), (79, 153), (69, 102)]
 
 
def rotamer_sg_positions(resi, chain=CHAIN):
    """Mutate one residue to CYS, keep all rotamers as states, return list of SG xyz."""
    cmd.wizard("mutagenesis")
    cmd.refresh_wizard()
    cmd.get_wizard().set_mode("CYS")
    cmd.get_wizard().do_select(f"prot and chain {chain} and resi {resi}")
    cmd.refresh_wizard()
    # the wizard preview object is called "mutation"; it holds N rotamer states
    n = cmd.count_states("mutation")
    positions = []
    for s in range(1, max(n, 1) + 1):
        model = cmd.get_model("mutation and name SG", state=s)
        if model.atom:
            positions.append(tuple(model.atom[0].coord))
    cmd.set_wizard()  # close without applying
    return positions
 
 
def best_pair(resi_a, resi_b):
    cmd.reinitialize()
    cmd.load(PDB, "prot")
    cmd.remove("solvent")
    sgs_a = rotamer_sg_positions(resi_a)
    # reload fresh so residue A's preview doesn't interfere with B
    cmd.reinitialize()
    cmd.load(PDB, "prot")
    cmd.remove("solvent")
    sgs_b = rotamer_sg_positions(resi_b)
 
    if not sgs_a or not sgs_b:
        return None
    best = min(
        sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5
        for a, b in itertools.product(sgs_a, sgs_b)
    )
    return best
 
 
def verdict(d):
    if d is None:            return "NO ROTAMERS"
    if 1.8 <= d <= 2.5:      return "PASS"
    if d <= 3.5:             return "BORDERLINE"
    return "FAIL"
 
 
if __name__ == "__main__":
    pymol.finish_launching(["pymol", "-qc"])
    print(f"{'pair':<14}{'best Sg-Sg':<14}{'verdict'}")
    for a, b in PAIRS:
        d = best_pair(a, b)
        ds = f"{d:.2f} A" if d is not None else "n/a"
        print(f"{a}C/{b}C{'':<6}{ds:<14}{verdict(d)}")
 