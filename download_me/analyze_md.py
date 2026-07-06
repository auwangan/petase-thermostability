import mdtraj as md
import numpy as np
import glob, os
os.makedirs("data/md/analysis", exist_ok=True)
results = []
for dcd in sorted(glob.glob("data/md/prod_*.dcd")):
    name = os.path.basename(dcd).replace("prod_","").replace(".dcd","")
    this_top = "data/md/mut_minimized.pdb" if name.startswith("mut") else "data/md/wt_minimized.pdb"
    if os.path.getsize(dcd) == 0:
        print("skip "+name); continue
    print("analyzing "+name)
    t = md.load(dcd, top=this_top)
    ca = t.topology.select("name CA")
    t.superpose(t, 0, atom_indices=ca)
    rmsd = md.rmsd(t, t, 0, atom_indices=ca) * 10
    rg = md.compute_rg(t) * 10
    plateau = float(rmsd[-len(rmsd)//4:].mean())
    results.append((name, float(rmsd.mean()), plateau, float(rg.mean())))
    print("  "+name+" meanRMSD %.2f plateau %.2f Rg %.2f" % (rmsd.mean(), plateau, rg.mean()))
with open("data/md/analysis/summary.csv","w") as f:
    f.write("run,mean_rmsd_A,plateau_rmsd_A,mean_rg_A\n")
    for r in results:
        f.write("%s,%.3f,%.3f,%.3f\n" % r)
print("DONE")
