import argparse, os, sys, math, urllib.request
from pdbfixer import PDBFixer
from openmm.app import PDBFile, ForceField, Modeller, Simulation, PME, HBonds
from openmm import LangevinMiddleIntegrator, unit, Platform
 
# --- Stockinger 5.2 parameters ---
FF        = ("amber14-all.xml", "amber14/tip4pew.xml")   # GAFF/Amber14 + tip4p-Ew
WATER     = "tip4pew"
PAD       = 1.0 * unit.nanometer                          # "Box padding of 1 nm ... cubic box"
IONIC     = 0.1 * unit.molar                              # "ionic strength of 0.1 M NaCl"
PH        = 8.0                                           # "utilizing a pH 8"
MIN_TOL   = 10.0 * unit.kilojoule_per_mole                # "minimized until 10 kJ/mole tolerance"
 
NATIVE_SS = [(203, 239), (273, 289)]   # our numbering; verify against crystal numbering
 
 
def fetch_pdb(pdb_id, path):
    if os.path.exists(path):
        print(f"[fetch] {path} already present"); return
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"[fetch] {url}")
    urllib.request.urlretrieve(url, path)
 
 
def ss_bonds(topology):
    out = []
    for a1, a2 in topology.bonds():
        if a1.name == "SG" and a2.name == "SG":
            out.append(tuple(sorted((int(a1.residue.id), int(a2.residue.id)))))
    return sorted(out)
 
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-id", required=True, help="e.g. 5XJH or 7SH6")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--outdir", default="data/v3/structures")
    ap.add_argument("--platform", default="CPU")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
 
    raw = f"{args.outdir}/{args.pdb_id.lower()}_raw.pdb"
    fetch_pdb(args.pdb_id, raw)
 
    print("=" * 60)
    print(f"  V3 PREP [{args.tag}]  <- {args.pdb_id} (Stockinger protocol)")
    print(f"  water {WATER} | ionic {IONIC} | pH {PH} | padding {PAD}")
    print("=" * 60)
 
    # --- clean crystal: drop waters/ligands, fill missing atoms, protonate at pH 8 ---
    fixer = PDBFixer(filename=raw)
    fixer.removeHeterogens(keepWater=False)      # strip crystallographic waters + ligands
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(PH)                # approximates their PROPKA/PDB2PQR pH 8
 
    ss = ss_bonds(fixer.topology)
    print(f"\n[fix] SG-SG bonds found: {ss}")
    n_res = sum(1 for _ in fixer.topology.residues())
    print(f"[fix] residues: {n_res}")
 
    # --- solvate (tip4p-Ew, cubic, 1 nm pad, 0.1 M NaCl) ---
    ff = ForceField(*FF)
    modeller = Modeller(fixer.topology, fixer.positions)
    print(f"\n[solvate] {WATER}, padding {PAD}, ionicStrength {IONIC}")
    modeller.addSolvent(ff, model=WATER, padding=PAD,
                        ionicStrength=IONIC, neutralize=True)
    n_tot = modeller.topology.getNumAtoms()
    n_prot = sum(1 for a in modeller.topology.atoms()
                 if a.residue.name not in ("HOH","WAT","NA","CL"))
    print(f"      protein {n_prot:,} / total {n_tot:,} particles")
    print(f"      (tip4p-Ew has a virtual site per water -> ~33% more particles than tip3p)")
 
    # --- minimise to 10 kJ/mol tolerance ---
    system = ff.createSystem(modeller.topology, nonbondedMethod=PME,
                             nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
    integ = LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond,
                                     0.002*unit.picoseconds)
    sim = Simulation(modeller.topology, system, integ,
                     Platform.getPlatformByName(args.platform))
    sim.context.setPositions(modeller.positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"\n[minimise] before {e0}")
    sim.minimizeEnergy(tolerance=MIN_TOL)
    st = sim.context.getState(getEnergy=True, getPositions=True)
    print(f"           after  {st.getPotentialEnergy()}")
 
    out = f"{args.outdir}/{args.tag}_v3_min.pdb"
    with open(out, "w") as fh:
        PDBFile.writeFile(modeller.topology, st.getPositions(), fh)
    print(f"\nWrote {out}")
    print(f"  -> next: md_v3_equilibrate.py --structure {out} --tag {args.tag} --temp 303.15")
 
 
if __name__ == "__main__":
    main()