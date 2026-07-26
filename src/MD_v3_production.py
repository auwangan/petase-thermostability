import argparse, os, sys, time
from openmm.app import (PDBFile, ForceField, Simulation, PME, HBonds,
                        StateDataReporter, DCDReporter, CheckpointReporter)
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit
 
FF = ("amber14-all.xml", "amber14/tip4pew.xml")
SOLVENT = ("HOH", "WAT", "NA", "CL", "K", "MG")
 
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structure", required=True, help="equilibrated PDB")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--temp", type=float, required=True)
    ap.add_argument("--rep", type=int, required=True)
    ap.add_argument("--ns", type=float, default=100.0)
    ap.add_argument("--frame-ps", type=float, default=5.0,
                    help="ps between saved frames (paper: 2.0)")
    ap.add_argument("--outdir", default="data/v3/md")
    ap.add_argument("--platform", default="CUDA")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
 
    T = a.temp * unit.kelvin
    report_every = int(a.frame_ps / 0.002)
    n_steps = int(a.ns * 1000 / 0.002)
    out = f"{a.outdir}/v3_{a.tag}_{int(round(a.temp))}K_rep{a.rep}"
 
    pdb = PDBFile(a.structure)
    ff = ForceField(*FF)
    system = ff.createSystem(pdb.topology, nonbondedMethod=PME,
                             nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
    system.addForce(MonteCarloBarostat(1.0*unit.atmosphere, T))
 
    protein = [at.index for at in pdb.topology.atoms() if at.residue.name not in SOLVENT]
    n_frames = n_steps // report_every
    est_mb = n_frames * len(protein) * 12 / 1e6
 
    print("=" * 60)
    print(f"  V3 PRODUCTION [{a.tag}] {a.temp} K  rep {a.rep}")
    print(f"  structure    : {a.structure}")
    print(f"  production   : {a.ns} ns ({n_steps:,} steps @ 2 fs)")
    print(f"  frames       : every {a.frame_ps} ps -> {n_frames:,} frames")
    print(f"  trajectory   : PROTEIN ONLY ({len(protein):,} of {pdb.topology.getNumAtoms():,})")
    print(f"  est size     : ~{est_mb:.0f} MB")
    print("=" * 60)
 
    integ = LangevinMiddleIntegrator(T, 1.0/unit.picosecond, 0.002*unit.picoseconds)
    sim = Simulation(pdb.topology, system, integ,
                     Platform.getPlatformByName(a.platform))
    sim.context.setPositions(pdb.positions)
    sim.context.setVelocitiesToTemperature(T)      # per-replicate randomness
 
    sim.reporters.append(StateDataReporter(sys.stdout, report_every, step=True,
        potentialEnergy=True, temperature=True, progress=True,
        totalSteps=n_steps, speed=True, remainingTime=True))
    sim.reporters.append(DCDReporter(out + ".dcd", report_every, atomSubset=protein))
    sim.reporters.append(CheckpointReporter(out + ".chk", 250000))
 
    t0 = time.time()
    sim.step(n_steps)
    print(f"  done in {(time.time()-t0)/3600:.2f} h")
 
    with open(out + "_final.pdb", "w") as fh:
        PDBFile.writeFile(pdb.topology, sim.context.getState(getPositions=True).getPositions(), fh)
    print(f"Wrote {out}.dcd")
 
 
if __name__ == "__main__":
    main()