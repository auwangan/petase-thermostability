import argparse, os, time, sys
from openmm.app import (PDBFile, ForceField, Simulation, PME, HBonds,
                        StateDataReporter)
from openmm import (LangevinMiddleIntegrator, MonteCarloBarostat,
                    CustomExternalForce, Platform, unit)
 
FF        = ("amber14-all.xml", "amber14/tip4pew.xml")
K_RESTR   = 100.0 * unit.kilojoule_per_mole / unit.angstrom**2   # their value
BACKBONE  = ("N", "CA", "C", "O")
NS_NVT    = 5.0
NS_NPT    = 5.0
NS_FREE   = 5.0
 
 
def add_backbone_restraint(system, topology, positions, k):
    """Harmonic position restraint on backbone atoms (periodic-safe)."""
    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", k)
    for p in ("x0", "y0", "z0"):
        force.addPerParticleParameter(p)
    n = 0
    for atom in topology.atoms():
        if atom.name in BACKBONE and atom.residue.name not in ("HOH","WAT","NA","CL"):
            force.addParticle(atom.index, positions[atom.index].value_in_unit(unit.nanometer))
            n += 1
    idx = system.addForce(force)
    return idx, n
 
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structure", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--temp", type=float, required=True, help="303.15 or 323.15")
    ap.add_argument("--outdir", default="data/v3/equil")
    ap.add_argument("--platform", default="CUDA")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    T = args.temp * unit.kelvin
    tag = f"{args.tag}_{int(round(args.temp))}K"
 
    pdb = PDBFile(args.structure)
    ff = ForceField(*FF)
    system = ff.createSystem(pdb.topology, nonbondedMethod=PME,
                             nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
    r_idx, n_restr = add_backbone_restraint(system, pdb.topology, pdb.positions, K_RESTR)
 
    print("=" * 60)
    print(f"  V3 EQUILIBRATION [{tag}]")
    print(f"  structure   : {args.structure}")
    print(f"  temperature : {args.temp} K")
    print(f"  restraint   : {K_RESTR} on {n_restr} backbone atoms")
    print(f"  stages      : {NS_NVT} ns NVT(restr) -> {NS_NPT} ns NPT(restr) -> {NS_FREE} ns free")
    print(f"  total       : {NS_NVT+NS_NPT+NS_FREE} ns")
    print("=" * 60)
 
    integ = LangevinMiddleIntegrator(T, 1.0/unit.picosecond, 0.002*unit.picoseconds)
    sim = Simulation(pdb.topology, system, integ,
                     Platform.getPlatformByName(args.platform))
    sim.context.setPositions(pdb.positions)
    sim.context.setVelocitiesToTemperature(T)
    sim.reporters.append(StateDataReporter(sys.stdout, 25000, step=True,
        potentialEnergy=True, temperature=True, speed=True))
 
    t0 = time.time()
    # --- stage 1: NVT, restrained ---
    print(f"\n[1/3] NVT restrained, {NS_NVT} ns ...")
    sim.step(int(NS_NVT * 1000 / 0.002))
 
    # --- stage 2: NPT, restrained ---
    print(f"\n[2/3] NPT restrained (MC barostat 1 atm), {NS_NPT} ns ...")
    system.addForce(MonteCarloBarostat(1.0*unit.atmosphere, T))
    sim.context.reinitialize(preserveState=True)
    sim.step(int(NS_NPT * 1000 / 0.002))
 
    # --- stage 3: free NPT ---
    print(f"\n[3/3] free NPT (restraints removed), {NS_FREE} ns ...")
    sim.context.setParameter("k", 0.0)      # switch restraint off
    sim.step(int(NS_FREE * 1000 / 0.002))
 
    print(f"\n  equilibration done in {(time.time()-t0)/3600:.2f} h")
 
    out = f"{args.outdir}/{tag}_equil.pdb"
    st = sim.context.getState(getPositions=True)
    with open(out, "w") as fh:
        PDBFile.writeFile(pdb.topology, st.getPositions(), fh)
    print(f"Wrote {out}")
    print(f"  -> next: md_v3_production.py --structure {out} --tag {args.tag} "
          f"--temp {args.temp} --rep 1..5")
 
 
if __name__ == "__main__":
    main()