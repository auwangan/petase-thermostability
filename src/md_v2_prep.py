#!/usr/bin/env python3
"""V2 prep: ONE script for both wt and mutant so preparation can never differ."""
import argparse, os, sys, math, time
from pdbfixer import PDBFixer
from openmm.app import PDBFile, ForceField, Modeller, Simulation, PME, HBonds
from openmm import LangevinMiddleIntegrator, unit, Platform

ENGINEERED_SS = (79, 153)
NATIVE_SS = [(203, 239), (273, 289)]
FF = ("amber14-all.xml", "amber14/tip3pfb.xml")
PAD = 1.0 * unit.nanometer
IONIC = 0.15 * unit.molar

def ss_bonds(topology):
    out = []
    for a1, a2 in topology.bonds():
        if a1.name == "SG" and a2.name == "SG":
            out.append(tuple(sorted((int(a1.residue.id), int(a2.residue.id)))))
    return sorted(out)

def sg_dist(positions, topology, ra, rb):
    idx = {}
    for atom in topology.atoms():
        if atom.name == "SG" and int(atom.residue.id) in (ra, rb):
            idx[int(atom.residue.id)] = atom.index
    if ra not in idx or rb not in idx:
        return None
    p1 = positions[idx[ra]].value_in_unit(unit.angstrom)
    p2 = positions[idx[rb]].value_in_unit(unit.angstrom)
    return math.dist(p1, p2)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--tag", required=True, choices=["wt", "mut"])
    ap.add_argument("--outdir", default="data/v2/structures")
    ap.add_argument("--platform", default="CPU")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    expect = list(NATIVE_SS) + ([tuple(sorted(ENGINEERED_SS))] if args.tag == "mut" else [])
    print("=" * 58)
    print(f"  V2 PREP  [{args.tag}]   input: {args.inp}")
    print(f"  expecting disulfides: {expect}")
    print("=" * 58)

    fixer = PDBFixer(filename=args.inp)
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)

    ss = ss_bonds(fixer.topology)
    print(f"\n[fix] SG-SG bonds: {ss}")
    missing = [w for w in expect if w not in ss]
    if missing:
        sys.exit(f"ABORT after fixing: disulfide(s) {missing} not formed.")
    for (ra, rb) in expect:
        print(f"      {ra}-{rb}: SG-SG = {sg_dist(fixer.positions, fixer.topology, ra, rb):.2f} A")

    ff = ForceField(*FF)
    modeller = Modeller(fixer.topology, fixer.positions)
    print(f"\n[solvate] padding {PAD}, ionicStrength {IONIC}, TIP3P")
    modeller.addSolvent(ff, model="tip3p", padding=PAD, ionicStrength=IONIC, neutralize=True)
    n_total = modeller.topology.getNumAtoms()
    n_prot = sum(1 for a in modeller.topology.atoms()
                 if a.residue.name not in ("HOH","WAT","NA","CL"))
    print(f"      protein atoms {n_prot:,} / total {n_total:,}")

    system = ff.createSystem(modeller.topology, nonbondedMethod=PME,
                             nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
    sim = Simulation(modeller.topology, system, integrator,
                     Platform.getPlatformByName(args.platform))
    sim.context.setPositions(modeller.positions)

    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"\n[minimize] energy before: {e0}")
    t = time.time(); sim.minimizeEnergy()
    state = sim.context.getState(getEnergy=True, getPositions=True)
    print(f"      energy after : {state.getPotentialEnergy()}   ({time.time()-t:.0f} s)")

    print("\n[check] disulfides after minimization:")
    bad = False
    for (ra, rb) in expect:
        d = sg_dist(state.getPositions(), modeller.topology, ra, rb)
        ok = d is not None and d < 2.5
        bad |= not ok
        print(f"      {ra}-{rb}: SG-SG = {d:.2f} A   {'intact' if ok else '*** STRAINED/BROKEN ***'}")
    if bad:
        sys.exit("ABORT: a disulfide is strained or broken after minimization.")

    out = f"{args.outdir}/{args.tag}_v2_minimized.pdb"
    with open(out, "w") as f:
        PDBFile.writeFile(modeller.topology, state.getPositions(), f)
    print(f"\nWrote {out}")

if __name__ == "__main__":
    main()
