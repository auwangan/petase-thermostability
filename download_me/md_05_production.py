import sys, time, os
from openmm.app import (PDBFile, ForceField, Simulation, PME, HBonds,
                        StateDataReporter, DCDReporter, CheckpointReporter)
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit
import argparse

p = argparse.ArgumentParser()
p.add_argument("--structure", required=True, help="minimized PDB (e.g. data/md/wt_minimized.pdb)")
p.add_argument("--temp", type=float, required=True, help="temperature in K (300 or 400)")
p.add_argument("--ns", type=float, required=True, help="production length in ns")
p.add_argument("--rep", type=int, required=True, help="replicate number (1,2,3)")
p.add_argument("--tag", required=True, help="label, e.g. wt or mut")
p.add_argument("--platform", default="CUDA", help="CUDA or CPU")
args = p.parse_args()

T = args.temp * unit.kelvin
out = f"data/md/prod_{args.tag}_{int(args.temp)}K_rep{args.rep}"
os.makedirs("data/md", exist_ok=True)

pdb = PDBFile(args.structure)
ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
system = ff.createSystem(pdb.topology, nonbondedMethod=PME,
                         nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
system.addForce(MonteCarloBarostat(1.0*unit.atmosphere, T))

integrator = LangevinMiddleIntegrator(T, 1.0/unit.picosecond, 0.002*unit.picoseconds)
platform = Platform.getPlatformByName(args.platform)
sim = Simulation(pdb.topology, system, integrator, platform)
sim.context.setPositions(pdb.positions)

# fresh random velocities at target T -> this is also the per-replicate randomness
sim.context.setVelocitiesToTemperature(T)

# --- equilibrate at target T (100 ps) ---
print(f"[{args.tag} {int(args.temp)}K rep{args.rep}] equilibrating 100 ps at {args.temp} K...")
sim.minimizeEnergy()                     # quick re-min in case T differs from prep
sim.context.setVelocitiesToTemperature(T)
sim.step(50000)                          # 50,000 * 2 fs = 100 ps

# --- production ---
n_steps = int(args.ns * 1000 / 0.002)    # ns -> steps (2 fs each)
report_every = 50000                      # log/trajectory frame every 10 ps
sim.reporters.append(StateDataReporter(sys.stdout, report_every, step=True,
    potentialEnergy=True, temperature=True, progress=True,
    totalSteps=n_steps + 50000, speed=True, remainingTime=True))
sim.reporters.append(DCDReporter(out + ".dcd", report_every))
sim.reporters.append(CheckpointReporter(out + ".chk", 50000))   # restart safety

print(f"[{args.tag} {int(args.temp)}K rep{args.rep}] production {args.ns} ns ({n_steps} steps)...")
t0 = time.time()
sim.step(n_steps)
print(f"  done in {(time.time()-t0)/3600:.2f} h")

with open(out + "_final.pdb", "w") as f:
    PDBFile.writeFile(pdb.topology, sim.context.getState(getPositions=True).getPositions(), f)
print(f"Wrote {out}.dcd and {out}_final.pdb")