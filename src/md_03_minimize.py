from openmm.app import PDBFile, ForceField, Modeller, Simulation, PME, HBonds
from openmm import LangevinMiddleIntegrator, unit
import time

# rebuild the system from the solvated structure
pdb = PDBFile("data/md/wt_solvated.pdb")
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

# create the System: the full physics model OpenMM will simulate
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,                 # how long-range electrostatics are computed
    nonbondedCutoff=1.0*unit.nanometer,  # direct-space cutoff for short-range forces
    constraints=HBonds,                  # freeze H-bond lengths -> allows 2 fs timesteps
)

# integrator = the rule for stepping atoms forward in time
integrator = LangevinMiddleIntegrator(
    300*unit.kelvin,                     # target temperature
    1.0/unit.picosecond,                 # friction (thermostat coupling)
    0.002*unit.picoseconds,              # 2 fs timestep
)

# Simulation ties topology + system + integrator together (runs on CPU here)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# energy BEFORE minimization
state0 = simulation.context.getState(getEnergy=True)
e0 = state0.getPotentialEnergy()
print(f"Energy before: {e0}")

# the minimization itself
print("Minimizing...")
t = time.time()
simulation.minimizeEnergy()
print(f"  done in {time.time()-t:.1f} s")

# energy AFTER
state1 = simulation.context.getState(getEnergy=True, getPositions=True)
e1 = state1.getPotentialEnergy()
print(f"Energy after:  {e1}")

# save the minimized structure
with open("data/md/wt_minimized.pdb", "w") as f:
    PDBFile.writeFile(pdb.topology, state1.getPositions(), f)
print("Wrote data/md/wt_minimized.pdb")