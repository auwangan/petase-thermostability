from openmm.app import PDBFile, ForceField, Modeller, PME, PDBxFile
from openmm import unit

# load the prepped (hydrogen-added) structure from step 1
pdb = PDBFile("data/md/mut_fixed.pdb")

# force field = the physics rulebook: how every atom-type interacts
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

modeller = Modeller(pdb.topology, pdb.positions)

# add water box + neutralising ions
modeller.addSolvent(
    forcefield,
    model="tip3p",                       # water model
    padding=1.0*unit.nanometer,          # 1 nm water shell around protein
    ionicStrength=0.15*unit.molar,       # ~physiological salt (0.15 M NaCl)
    neutralize=True,                     # add ions to cancel net charge
)

# save the solvated system
with open("data/md/mut_solvated.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

n = modeller.topology.getNumAtoms()
print(f"Solvated system written to data/md/mut_solvated.pdb")
print(f"Total atoms now: {n}")