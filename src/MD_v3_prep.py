import argparse, os, sys, urllib.request
import numpy as np
from pdbfixer import PDBFixer
from openmm.app import PDBFile, ForceField, Modeller, Simulation, PME, HBonds
from openmm import LangevinMiddleIntegrator, unit, Platform
 
FF       = ("amber14-all.xml", "amber14/tip4pew.xml")
WATER    = "tip4pew"
PAD      = 1.0 * unit.nanometer
IONIC    = 0.1 * unit.molar
PH       = 8.0
MIN_TOL  = 10.0 * unit.kilojoule_per_mole / unit.nanometer   # FORCE units (OpenMM 8)
 
# common range across 5XJH (30-292) and 7SH6 (29-289)
RES_LO, RES_HI = 30, 289
MAX_PARTICLES  = 70000        # sane upper bound for a 265-res protein + 1nm pad
 
 
def fetch(pdb_id, path):
    if not os.path.exists(path):
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        print(f"[fetch] {url}")
        urllib.request.urlretrieve(url, path)
    else:
        print(f"[fetch] {path} present")
 
 
def drop_terminal_missing(fixer):
    """Keep internal gap-filling; discard missing residues at chain termini."""
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    removed = 0
    for key in list(fixer.missingResidues.keys()):
        chain = chains[key[0]]
        n_res = len(list(chain.residues()))
        if key[1] == 0 or key[1] == n_res:        # start or end of chain
            del fixer.missingResidues[key]
            removed += 1
    print(f"[fix] dropped {removed} terminal missing-residue block(s); "
          f"internal gaps to fill: {len(fixer.missingResidues)}")
 
 
def trim_range(fixer, lo, hi):
    """Delete residues outside [lo, hi] and any non-protein chains."""
    to_delete = []
    for chain in fixer.topology.chains():
        for res in chain.residues():
            try:
                rid = int(res.id)
            except ValueError:
                to_delete.append(res); continue
            if rid < lo or rid > hi:
                to_delete.append(res)
    if to_delete:
        fixer.removeChains([])          # no-op, keeps API happy
        modeller = Modeller(fixer.topology, fixer.positions)
        modeller.delete(to_delete)
        fixer.topology, fixer.positions = modeller.topology, modeller.positions
    print(f"[fix] trimmed to residues {lo}-{hi}: "
          f"{sum(1 for _ in fixer.topology.residues())} residues remain")
 
 
def ss_bonds(top):
    return sorted(tuple(sorted((int(a.residue.id), int(b.residue.id))))
                  for a, b in top.bonds() if a.name == "SG" and b.name == "SG")
 
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-id", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--outdir", default="data/v3/structures")
    ap.add_argument("--platform", default="CPU")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
 
    raw = f"{a.outdir}/{a.pdb_id.lower()}_raw.pdb"
    fetch(a.pdb_id, raw)
 
    print("=" * 62)
    print(f"  V3 PREP [{a.tag}] <- {a.pdb_id}   (Stockinger protocol, fixed)")
    print(f"  {WATER} | {IONIC} NaCl | pH {PH} | pad {PAD} | residues {RES_LO}-{RES_HI}")
    print("=" * 62)

    fixer = PDBFixer(filename=raw)
    fixer.removeHeterogens(keepWater=False)      # drop HOH, SO4, ligands
    drop_terminal_missing(fixer)                 # FIX 1
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(PH)
    # trim_range(fixer, RES_LO, RES_HI)  # disabled: breaks terminal capping            # FIX 3

    n_res = sum(1 for _ in fixer.topology.residues())
    ids = [int(r.id) for r in fixer.topology.residues() if r.id.isdigit()]
    print(f"[fix] residues {min(ids)}..{max(ids)}  (n={n_res})")
    print(f"[fix] SG-SG bonds: {ss_bonds(fixer.topology)}")
 
    pos = np.array(fixer.positions.value_in_unit(unit.nanometer))
    ext = pos.max(0) - pos.min(0)
    print(f"[fix] extent (nm): {ext.round(2)}")
    if ext.max() > 8.0:
        sys.exit(f"ABORT: extent {ext.max():.1f} nm too large - stray atoms remain.")
 
    ff = ForceField(*FF)
    modeller = Modeller(fixer.topology, fixer.positions)
    print(f"\n[solvate] {WATER}, pad {PAD}, {IONIC}")
    modeller.addSolvent(ff, model=WATER, padding=PAD,
                        ionicStrength=IONIC, neutralize=True)
    n_tot = modeller.topology.getNumAtoms()
    n_prot = sum(1 for at in modeller.topology.atoms()
                 if at.residue.name not in ("HOH","WAT","NA","CL"))
    print(f"      protein {n_prot:,} / total {n_tot:,} particles")
    if n_tot > MAX_PARTICLES:                    # FIX 4
        sys.exit(f"ABORT: {n_tot:,} particles > {MAX_PARTICLES:,}. "
                 f"Box blow-up - do not send this to the GPU.")
    system = ff.createSystem(modeller.topology, nonbondedMethod=PME,
                             nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
    sim = Simulation(modeller.topology, system,
                     LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond,
                                              0.002*unit.picoseconds),
                     Platform.getPlatformByName(a.platform))
    sim.context.setPositions(modeller.positions)
    print(f"\n[minimise] before {sim.context.getState(getEnergy=True).getPotentialEnergy()}")
    sim.minimizeEnergy(tolerance=MIN_TOL)        # FIX 2
    st = sim.context.getState(getEnergy=True, getPositions=True)
    print(f"           after  {st.getPotentialEnergy()}")
 
    out = f"{a.outdir}/{a.tag}_v3_min.pdb"
    with open(out, "w") as fh:
        PDBFile.writeFile(modeller.topology, st.getPositions(), fh)
    print(f"\nWrote {out}")

if __name__ == "__main__":
    main()