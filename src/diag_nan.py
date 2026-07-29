
import sys
from openmm.app import PDBFile, ForceField, Simulation, PME, HBonds
from openmm import LangevinMiddleIntegrator, Platform, unit

path = sys.argv[1]
pdb = PDBFile(path)
ff = ForceField("amber14-all.xml", "amber14/tip4pew.xml")

bv = pdb.topology.getPeriodicBoxVectors()
print(f"file        : {path}")
print(f"particles   : {pdb.topology.getNumAtoms():,}")
if bv is None:
    print("BOX VECTORS : *** NONE (no CRYST1 record) *** <- this alone causes NaN with PME")
else:
    d = [bv[i][i].value_in_unit(unit.nanometer) for i in range(3)]
    vol = d[0]*d[1]*d[2]
    print(f"box (nm)    : {[round(x,3) for x in d]}   volume {vol:.1f} nm^3")
    n_wat = sum(1 for r in pdb.topology.residues() if r.name in ("HOH","WAT"))
    print(f"waters      : {n_wat:,}  -> density {n_wat/vol:.1f} /nm^3 (bulk ~33.4)")

system = ff.createSystem(pdb.topology, nonbondedMethod=PME,
                         nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
sim = Simulation(pdb.topology, system,
                 LangevinMiddleIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond,
                                          0.002*unit.picoseconds),
                 Platform.getPlatformByName("CUDA"))
sim.context.setPositions(pdb.positions)
sim.context.computeVirtualSites()
e = sim.context.getState(getEnergy=True).getPotentialEnergy()
print(f"\nenergy as loaded : {e}")
ev = e.value_in_unit(unit.kilojoule_per_mole)
if ev != ev:
    print("  -> NaN ALREADY AT LOAD. Coordinates or box are broken.")
elif ev > -300000:
    print("  -> far above the ~-515,000 seen during equilibration: CLASH present.")
else:
    print("  -> looks sane; NaN must arise during dynamics.")

print("\ntrying minimisation...")
sim.minimizeEnergy(maxIterations=200)
e2 = sim.context.getState(getEnergy=True).getPotentialEnergy()
print(f"after minimise   : {e2}")