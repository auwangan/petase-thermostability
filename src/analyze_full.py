import mdtraj as md
import numpy as np
import glob, os

os.makedirs("data/md/analysis", exist_ok=True)
summary = []
for dcd in sorted(glob.glob("data/md/prod_*.dcd")):
    name = os.path.basename(dcd).replace("prod_", "").replace(".dcd", "")
    this_top = "data/md/mut_minimized.pdb" if name.startswith("mut") else "data/md/wt_minimized.pdb"
    if os.path.getsize(dcd) == 0:
        print("skip (empty):", name); continue
    print("analyzing", name)
    t = md.load(dcd, top=this_top)
    ca = t.topology.select("name CA")
    t.superpose(t, 0, atom_indices=ca)

    rmsd = md.rmsd(t, t, 0, atom_indices=ca) * 10        # per-frame, Angstrom
    rmsf = md.rmsf(t, t, 0, atom_indices=ca) * 10        # per-residue, Angstrom
    rg   = md.compute_rg(t) * 10                          # per-frame, Angstrom

    # residue numbers for the CA atoms (so RMSF maps to real residue IDs)
    resids = [t.topology.atom(i).residue.resSeq for i in ca]

    np.savetxt(f"data/md/analysis/{name}_rmsd.csv", rmsd, header="rmsd_A", comments="")
    np.savetxt(f"data/md/analysis/{name}_rg.csv", rg, header="rg_A", comments="")
    np.savetxt(f"data/md/analysis/{name}_rmsf.csv",
               np.column_stack([resids, rmsf]),
               delimiter=",", header="resid,rmsf_A", comments="")

    plateau = float(rmsd[-len(rmsd)//4:].mean())
    summary.append((name, float(rmsd.mean()), plateau, float(rg.mean())))
    print(f"  {name}: meanRMSD {rmsd.mean():.2f} plateau {plateau:.2f} Rg {rg.mean():.2f} (RMSF saved per-residue)")

with open("data/md/analysis/summary.csv", "w") as f:
    f.write("run,mean_rmsd_A,plateau_rmsd_A,mean_rg_A\n")
    for r in summary:
        f.write("%s,%.3f,%.3f,%.3f\n" % r)
print("DONE - RMSD, RMSF (per-residue), Rg all saved")