import sys, os, shutil
 
TARGET = "src/analyse_v2.py"
 
SASA_BLOCK = '''
    # --- SASA (Stockinger 5.4: MDTraj, every 100th frame) ---
    print("  computing SASA (stride 100)...")
    stride = 100 if traj.n_frames > 200 else 1
    sasa_res = md.shrake_rupley(traj[::stride], mode="residue")   # nm^2
    sasa_total = sasa_res.sum(axis=1)          # total SASA per analysed frame
    sasa_mean_res = sasa_res.mean(axis=0)      # mean SASA per residue
'''
 
SASA_WRITE = '''    np.savetxt(f"{o}/{args.tag}_sasa.csv",
               np.c_[np.arange(len(sasa_total))*stride, sasa_total],
               delimiter=",", header="frame,sasa_total_nm2", comments="")
    with open(f"{o}/{args.tag}_sasa_perres.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["resid", "sasa_nm2"])
        for r, v in zip(resids, sasa_mean_res): w.writerow([r, f"{v:.4f}"])
'''
 
def main():
    if not os.path.exists(TARGET):
        sys.exit(f"ERROR: {TARGET} not found. Run from the repo root.")
    src = open(TARGET).read()
 
    if "shrake_rupley" in src:
        sys.exit("Already patched (shrake_rupley present). Nothing to do.")
 
    # anchor 1: insert the computation just before the write section
    a1 = "    # --- write everything ---"
    if a1 not in src:
        sys.exit(f"ERROR: anchor not found: {a1!r}\n"
                 f"Your analyse_v2.py differs from expected - patch manually.")
    src = src.replace(a1, SASA_BLOCK + "\n" + a1, 1)
 
    # anchor 2: add the CSV writes right after the Q csv write
    a2 = '    np.savetxt(f"{o}/{args.tag}_rmsd.csv"'
    if a2 not in src:
        sys.exit(f"ERROR: anchor not found: {a2!r}")
    src = src.replace(a2, SASA_WRITE + a2, 1)
 
    # anchor 3: summary header
    a3 = '"mean_rg_A", "n_contacts", "n_frames"'
    if a3 in src:
        src = src.replace(a3, '"mean_rg_A", "mean_sasa_nm2", "n_contacts", "n_frames"', 1)
    else:
        print("  warn: summary header anchor not found, skipping header update")
 
    # anchor 4: summary row
    a4 = 'f"{rg.mean():.3f}", len(pairs), traj.n_frames])'
    if a4 in src:
        src = src.replace(a4,
            'f"{rg.mean():.3f}", f"{sasa_total.mean():.2f}", len(pairs), traj.n_frames])', 1)
    else:
        print("  warn: summary row anchor not found, skipping row update")
 
    shutil.copy(TARGET, TARGET + ".bak")
    open(TARGET, "w").write(src)
    print(f"patched {TARGET}  (backup: {TARGET}.bak)")
    print("  + md.shrake_rupley, stride 100")
    print("  + writes <tag>_sasa.csv and <tag>_sasa_perres.csv")
    print("  + mean_sasa_nm2 column in summary_v2.csv")
 
if __name__ == "__main__":
    main()