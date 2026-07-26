import argparse, os, glob
import mdtraj as md
 
ap = argparse.ArgumentParser()
ap.add_argument("--tag", required=True, choices=["wt", "mut", "fast"])
ap.add_argument("--nframes", type=int, default=10, help="total frames across all reps")
ap.add_argument("--last-ns", type=float, default=30.0, help="sample from final N ns")
ap.add_argument("--outdir", default="data/v2/ensemble")
args = ap.parse_args()
 
TOP = f"data/v2/structures/{args.tag}_v2_top.pdb"
trajs = sorted(glob.glob(f"data/md/v2_prod_{args.tag}_350K_rep[1-5].dcd"))
if not trajs:
    raise SystemExit(f"no trajectories found for {args.tag}")
os.makedirs(args.outdir, exist_ok=True)
 
per_rep = max(1, args.nframes // len(trajs))
print(f"[{args.tag}] {len(trajs)} reps, taking {per_rep} frame(s) each "
      f"from the last {args.last_ns} ns")
 
n = 0
for ti, tpath in enumerate(trajs, start=1):
    t = md.load(tpath, top=TOP)
    dt_ns = (t.time[1] - t.time[0]) / 1000.0 if t.n_frames > 1 else 0.02
    n_last = min(t.n_frames, max(1, int(args.last_ns / dt_ns)))
    tail = t[-n_last:]
    # evenly spaced picks within the equilibrated tail
    idxs = [int(i * (tail.n_frames - 1) / max(1, per_rep - 1)) for i in range(per_rep)] \
           if per_rep > 1 else [tail.n_frames - 1]
    for k, fi in enumerate(idxs):
        out = f"{args.outdir}/{args.tag}_rep{ti}_f{k}.pdb"
        tail[fi].save_pdb(out)
        n += 1
        print(f"  wrote {out}  (frame {fi} of tail, t={tail.time[fi]/1000:.1f} ns)")
print(f"[{args.tag}] {n} frames written to {args.outdir}")