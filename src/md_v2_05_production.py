import sys, time, os, argparse, math
from openmm.app import (PDBFile, ForceField, Simulation, PME, HBonds,
                        StateDataReporter, DCDReporter, CheckpointReporter)
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit
 
# engineered disulfide (mutant only); natives present in both
ENGINEERED_SS = (79, 153)
NATIVE_SS = [(203, 239), (273, 289)]
 
 
def find_ss_bonds(topology):
    """Return list of (resid_a, resid_b) for every SG-SG bond in the topology."""
    out = []
    for a1, a2 in topology.bonds():
        if a1.name == "SG" and a2.name == "SG":
            out.append(tuple(sorted((int(a1.residue.id), int(a2.residue.id)))))
    return sorted(out)
 
 
def sg_distance(positions, topology, ra, rb):
    """SG-SG distance in Angstrom, or None if either SG is missing."""
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
    p = argparse.ArgumentParser()
    p.add_argument("--structure", required=True)
    p.add_argument("--temp", type=float, required=True)
    p.add_argument("--ns", type=float, required=True)
    p.add_argument("--rep", type=int, required=True)
    p.add_argument("--tag", required=True, choices=["wt", "mut"])
    p.add_argument("--platform", default="CUDA")
    p.add_argument("--report-every", type=int, default=10000,
                   help="steps between frames (10000 * 2fs = 20 ps)")
    p.add_argument("--equil-ns", type=float, default=1.0,
                   help="equilibration at target T before production")
    args = p.parse_args()
 
    T = args.temp * unit.kelvin
    out = f"data/md/v2_prod_{args.tag}_{int(args.temp)}K_rep{args.rep}"
    os.makedirs("data/md", exist_ok=True)
 
    pdb = PDBFile(args.structure)
    ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    system = ff.createSystem(pdb.topology, nonbondedMethod=PME,
                             nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
    system.addForce(MonteCarloBarostat(1.0*unit.atmosphere, T))
 
    # ---------- GUARD 1: protein-only trajectory ----------
    protein_atoms = [a.index for a in pdb.topology.atoms() if a.residue.name not in
                     ("HOH", "WAT", "NA", "CL", "K", "MG")]
    n_frames = int(args.ns * 1000 / 0.002 / args.report_every)
    est_mb = n_frames * len(protein_atoms) * 12 / 1e6
 
    # ---------- GUARD 3a: disulfides BEFORE the run ----------
    ss = find_ss_bonds(pdb.topology)
    want = list(NATIVE_SS) + ([tuple(sorted(ENGINEERED_SS))] if args.tag == "mut" else [])
    d_eng = (sg_distance(pdb.positions, pdb.topology, *ENGINEERED_SS)
             if args.tag == "mut" else None)
 
    # ---------- GUARD 2: announce everything BEFORE doing work ----------
    n_steps = int(args.ns * 1000 / 0.002)
    print("=" * 56)
    print(f"  V2 PRODUCTION  [{args.tag}]  rep {args.rep}")
    print("=" * 56)
    print(f"  structure      : {args.structure}")
    print(f"  temperature    : {args.temp} K")
    print(f"  production     : {args.ns} ns  ({n_steps:,} steps @ 2 fs)")
    print(f"  equilibration  : {args.equil_ns} ns at target T")
    print(f"  report_every   : {args.report_every} steps "
          f"({args.report_every*0.002:.0f} ps/frame, {n_frames:,} frames)")
    print(f"  trajectory     : PROTEIN ONLY "
          f"({len(protein_atoms):,} of {pdb.topology.getNumAtoms():,} atoms)")
    print(f"  est traj size  : ~{est_mb:.0f} MB")
    print(f"  SG-SG bonds    : {ss}")
    for w in want:
        ok = "YES" if w in ss else "*** MISSING ***"
        print(f"    disulfide {w[0]}-{w[1]:<4}: {ok}")
    if d_eng is not None:
        print(f"  engineered SG-SG distance (start): {d_eng:.2f} A")
    print("=" * 56)
 
    missing = [w for w in want if w not in ss]
    if missing:
        sys.exit(f"ABORT: expected disulfide(s) {missing} not in topology. "
                 f"Do not spend GPU time on a broken system.")
 
    platform = Platform.getPlatformByName(args.platform)
    print(f"  platform: {platform.getName()}")
    sim = Simulation(pdb.topology, system, 
                     LangevinMiddleIntegrator(T, 1.0/unit.picosecond, 0.002*unit.picoseconds),
                     platform)
    sim.context.setPositions(pdb.positions)
 
    # per-replicate randomness lives here
    sim.context.setVelocitiesToTemperature(T)
 
    print(f"\n[{args.tag} {int(args.temp)}K rep{args.rep}] equilibrating {args.equil_ns} ns...")
    sim.step(int(args.equil_ns * 1000 / 0.002))
 
    sim.reporters.append(StateDataReporter(sys.stdout, args.report_every, step=True,
        potentialEnergy=True, temperature=True, progress=True,
        totalSteps=n_steps, speed=True, remainingTime=True))
    sim.reporters.append(DCDReporter(out + ".dcd", args.report_every,
                                     atomSubset=protein_atoms))     # <-- GUARD 1
    sim.reporters.append(CheckpointReporter(out + ".chk", 100000))
 
    print(f"[{args.tag} {int(args.temp)}K rep{args.rep}] production {args.ns} ns...")
    t0 = time.time()
    sim.step(n_steps)
    print(f"  done in {(time.time()-t0)/3600:.2f} h")
 
    state = sim.context.getState(getPositions=True)
    with open(out + "_final.pdb", "w") as f:
        PDBFile.writeFile(pdb.topology, state.getPositions(), f)
 
    # ---------- GUARD 3b: disulfides AFTER the run ----------
    print("\n--- post-run disulfide check ---")
    for (ra, rb) in want:
        d = sg_distance(state.getPositions(), pdb.topology, ra, rb)
        verdict = "intact" if d and d < 2.5 else "*** BROKEN / STRAINED ***"
        print(f"  {ra}-{rb}: SG-SG = {d:.2f} A  {verdict}")
    print(f"\nWrote {out}.dcd ({est_mb:.0f} MB est) and {out}_final.pdb")
 
 
if __name__ == "__main__":
    main()