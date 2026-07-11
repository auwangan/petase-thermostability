import argparse, os, csv
import numpy as np
import mdtraj as md

CONTACT_CUTOFF = 0.45     # nm (4.5 A) -- native contact definition
SEQ_SEP        = 3        # |i-j| >= 3
BETA           = 1.2      # a contact is "kept" if r < BETA * r_native
AVG_LAST_NS    = 30.0     # average Q over final 30 ns


def native_contacts(ref):
    """Return (pairs, r_native) for heavy-atom contacts in the reference frame."""
    heavy = ref.topology.select("protein and not element H")
    ref_h = ref.atom_slice(heavy)

    # residue index of each heavy atom
    res_of = np.array([a.residue.index for a in ref_h.topology.atoms])

    # all heavy-atom pairs within cutoff in the reference
    pairs = md.compute_neighborlist(ref_h, CONTACT_CUTOFF, 0)
    contacts = set()
    for i, neigh in enumerate(pairs):
        for j in neigh:
            if j <= i:
                continue
            ri, rj = res_of[i], res_of[j]
            if abs(ri - rj) >= SEQ_SEP:
                contacts.add((min(i, j), max(i, j)))   # atom-level pair

    if not contacts:
        raise RuntimeError("no native contacts found - check topology/selection")

    atom_pairs = np.array(sorted(contacts))
    # map back to the FULL-trajectory atom indices
    atom_pairs = heavy[atom_pairs]
    r_native = md.compute_distances(ref, atom_pairs)[0]
    return atom_pairs, r_native


def compute_Q(traj, atom_pairs, r_native):
    """Q(t) = fraction of native contacts retained (hard cutoff at BETA*r_native)."""
    d = md.compute_distances(traj, atom_pairs)          # (n_frames, n_pairs)
    kept = d < (BETA * r_native)                        # broadcast per pair
    return kept.mean(axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True, help="topology / reference PDB (the MINIMIZED structure)")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--outdir", default="analysis_v2")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"[{args.tag}] loading...")
    traj = md.load(args.traj, top=args.top)
    ref = md.load(args.top)
    print(f"  {traj.n_frames} frames, {traj.n_atoms} atoms, {traj.n_residues} residues")

    # --- align on protein backbone before RMSD/RMSF ---
    prot = traj.topology.select("protein")
    traj = traj.atom_slice(prot)
    ref = ref.atom_slice(ref.topology.select("protein"))
    traj.superpose(ref)

    # --- Q ---
    print("  computing Q (native contacts)...")
    pairs, r_nat = native_contacts(ref)
    q = compute_Q(traj, pairs, r_nat)
    print(f"    {len(pairs)} native contacts")

    # frames covering the final AVG_LAST_NS
    if traj.n_frames > 1:
        dt_ns = (traj.time[1] - traj.time[0]) / 1000.0
        n_last = max(1, int(AVG_LAST_NS / dt_ns)) if dt_ns > 0 else traj.n_frames
    else:
        n_last = 1
    q_final = float(np.mean(q[-n_last:]))
    print(f"    mean Q (final {AVG_LAST_NS} ns, {n_last} frames) = {q_final:.4f}")

    # --- RMSD / Rg ---
    rmsd = md.rmsd(traj, ref) * 10.0        # nm -> Angstrom
    rg = md.compute_rg(traj) * 10.0

    # --- RMSF (per residue, CA) ---
    ca = traj.topology.select("name CA")
    rmsf = md.rmsf(traj.atom_slice(ca), ref.atom_slice(ref.topology.select("name CA"))) * 10.0
    resids = [traj.topology.atom(i).residue.resSeq for i in ca]

    # --- DSSP ---
    print("  computing DSSP...")
    dssp = md.compute_dssp(traj, simplified=True)   # (n_frames, n_residues)
    # fraction of frames each residue is helix/sheet
    frac_H = (dssp == "H").mean(axis=0)
    frac_E = (dssp == "E").mean(axis=0)

    # --- write everything ---
    o = args.outdir
    np.savetxt(f"{o}/{args.tag}_Q.csv", np.c_[traj.time/1000.0, q],
               delimiter=",", header="time_ns,Q", comments="")
    np.savetxt(f"{o}/{args.tag}_rmsd.csv", np.c_[traj.time/1000.0, rmsd],
               delimiter=",", header="time_ns,rmsd_A", comments="")
    np.savetxt(f"{o}/{args.tag}_rg.csv", np.c_[traj.time/1000.0, rg],
               delimiter=",", header="time_ns,rg_A", comments="")
    with open(f"{o}/{args.tag}_rmsf.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["resid", "rmsf_A"])
        for r, v in zip(resids, rmsf): w.writerow([r, f"{v:.4f}"])
    with open(f"{o}/{args.tag}_dssp.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["resid", "frac_helix", "frac_sheet"])
        for r, h, e in zip(resids, frac_H, frac_E): w.writerow([r, f"{h:.3f}", f"{e:.3f}"])

    # --- summary line ---
    summary = f"{o}/summary_v2.csv"
    new = not os.path.exists(summary)
    with open(summary, "a", newline="") as fh:
        w = csv.writer(fh)
        if new: w.writerow(["run", "mean_Q_final", "mean_rmsd_A", "mean_rg_A", "n_contacts", "n_frames"])
        w.writerow([args.tag, f"{q_final:.4f}", f"{rmsd.mean():.3f}",
                    f"{rg.mean():.3f}", len(pairs), traj.n_frames])

    print(f"  wrote Q / RMSD / Rg / RMSF / DSSP to {o}/  (appended {summary})")
    print(f"\n  >>> mean Q = {q_final:.4f}  <<<")
    if q_final > 0.9:
        print("      (high Q = fold well retained; expected for a stable 300K run)")
    elif q_final < 0.6:
        print("      (low Q = substantially unfolded)")


if __name__ == "__main__":
    main()