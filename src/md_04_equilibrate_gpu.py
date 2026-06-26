from openmm.app import PDBFile, ForceField, Simulation, PME, HBonds, StateDataReporter, DCDReporter
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, unit
import sys, time
from openmm import Platform
platform = Platform.getPlatformByName("CUDA")

pdb = PDBFile("data/md/wt_minimized.pdb")
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

system = forcefield.createSystem(
    pdb.topology, nonbondedMethod=PME,
    nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds,
)

# barostat: keeps pressure at 1 atm by letting the box volume breathe (NPT ensemble)
system.addForce(MonteCarloBarostat(1.0*unit.atmosphere, 300*unit.kelvin))

integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

# assign starting velocities drawn from the 300 K distribution
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)

# reporters: write progress to screen + trajectory to disk
simulation.reporters.append(StateDataReporter(
    sys.stdout, 500, step=True, potentialEnergy=True,
    temperature=True, progress=True, totalSteps=25000, speed=True,
))
simulation.reporters.append(DCDReporter("data/md/wt_equil.dcd", 500))

print("Equilibrating (50 ps)...")
t = time.time()
simulation.step(25000)          # 25,000 steps x 2 fs = 50 ps
print(f"  done in {time.time()-t:.1f} s")

with open("data/md/wt_equilibrated.pdb", "w") as f:
    PDBFile.writeFile(pdb.topology, simulation.context.getState(getPositions=True).getPositions(), f)
print("Wrote data/md/wt_equilibrated.pdb")