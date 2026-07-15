import argparse, os, glob
import numpy as np
import csv
 
IS_REGIONS = {
    "IS1": (54, 70),  "IS2": (78, 102), "IS3": (110, 134),
    "IS4": (150, 182),"IS5": (186, 198),"IS6": (206, 222),"IS7": (230, 246),
}
DISULFIDE = (79, 153)
TRIAD = (160, 206, 237)
 
 
def load_rmsf(path):
    resid, rmsf = [], []
    with open(path) as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            resid.append(int(row[0])); rmsf.append(float(row[1]))
    return np.array(resid), np.array(rmsf)
 
 
def load_q(path):
    t, q = [], []
    with open(path) as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            t.append(float(row[0])); q.append(float(row[1]))
    return np.array(t), np.array(q)
 
 
def mean_rmsf(analysis, tag):
    files = sorted(glob.glob(f"{analysis}/{tag}_350K_rep*_rmsf.csv"))
    if not files:
        raise SystemExit(f"no RMSF files for {tag} in {analysis}")
    resid = None; stack = []
    for f in files:
        r, v = load_rmsf(f)
        resid = r if resid is None else resid
        stack.append(v)
    return resid, np.mean(stack, axis=0), np.std(stack, axis=0), len(files)
 
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis", default="data/v2/analysis")
    ap.add_argument("--out", default="figures")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
 
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
 
    # ---- 1. per-residue dRMSF (mut - wt) ----
    resid, wt_m, wt_s, nwt = mean_rmsf(args.analysis, "wt")
    _,     mut_m, mut_s, nmut = mean_rmsf(args.analysis, "mut")
    d = mut_m - wt_m                       # negative = mutant MORE rigid
 
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.axhline(0, color="k", lw=0.6)
    ax.plot(resid, d, lw=1.1, color="crimson")
    for name, (a, b) in IS_REGIONS.items():
        ax.axvspan(a, b, color="grey", alpha=0.12)
        ax.text((a+b)/2, ax.get_ylim()[1]*0.85, name, ha="center", fontsize=7, color="grey")
    for r in DISULFIDE:
        ax.axvline(r, color="navy", ls="--", lw=1)
    ax.text(DISULFIDE[0], ax.get_ylim()[0]*0.9, "C79", color="navy", fontsize=8)
    ax.text(DISULFIDE[1], ax.get_ylim()[0]*0.9, "C153", color="navy", fontsize=8)
    for r in TRIAD:
        ax.axvline(r, color="green", ls=":", lw=0.8)
    ax.set_xlabel("residue"); ax.set_ylabel("dRMSF  mut - wt  (A)")
    ax.set_title(f"Per-residue flexibility change (mut - wt), {nwt}+{nmut} reps\n"
                 "negative = disulfide made it MORE rigid")
    fig.tight_layout(); fig.savefig(f"{args.out}/drmsf_mut_minus_wt.png", dpi=150)
    print(f"wrote {args.out}/drmsf_mut_minus_wt.png")
 
    # ---- 2. Q(t) overlay ----
    fig, ax = plt.subplots(figsize=(10, 5))
    for tag, col in (("wt", "steelblue"), ("mut", "crimson")):
        for f in sorted(glob.glob(f"{args.analysis}/{tag}_350K_rep*_Q.csv")):
            t, q = load_q(f)
            ax.plot(t, q, color=col, alpha=0.5, lw=0.8,
                    label=tag if f.endswith("rep1_Q.csv") else None)
    ax.set_xlabel("time (ns)"); ax.set_ylabel("Q (fraction native contacts)")
    ax.set_title("Q(t), all reps  (blue=WT, red=mutant)")
    ax.legend(); ax.set_ylim(0.6, 1.0)
    fig.tight_layout(); fig.savefig(f"{args.out}/Q_vs_time.png", dpi=150)
    print(f"wrote {args.out}/Q_vs_time.png")
 
    # ---- 3. numbers that explain the null ----
    def near(res, window=3):
        mask = np.zeros_like(resid, dtype=bool)
        for r in res:
            mask |= np.abs(resid - r) <= window
        return mask
 
    ds_mask = near(DISULFIDE)
    print("\n=== DIAGNOSIS ===")
    print(f"WT  global mean RMSF: {wt_m.mean():.3f} A")
    print(f"mut global mean RMSF: {mut_m.mean():.3f} A")
    print(f"  -> global dRMSF: {d.mean():+.3f} A (negative = mutant more rigid overall)")
    print(f"\nAt/near the disulfide (res {DISULFIDE[0]},{DISULFIDE[1]} +/-3):")
    print(f"  WT  RMSF: {wt_m[ds_mask].mean():.3f} A")
    print(f"  mut RMSF: {mut_m[ds_mask].mean():.3f} A")
    print(f"  -> local dRMSF: {d[ds_mask].mean():+.3f} A")
    if d[ds_mask].mean() < -0.05:
        print("  VERDICT: the disulfide DID rigidify its own region locally.")
    elif d[ds_mask].mean() > 0.05:
        print("  VERDICT: the disulfide region got MORE flexible (bond straining fold?).")
    else:
        print("  VERDICT: no local rigidification -- the staple isn't doing its job.")
 
    print("\nMost-flexible WT regions (where unfolding actually happens):")
    order = np.argsort(wt_m)[::-1][:8]
    for i in order:
        reg = next((n for n,(a,b) in IS_REGIONS.items() if a <= resid[i] <= b), "-")
        print(f"  res {resid[i]} ({reg}): WT RMSF {wt_m[i]:.2f} A")
    print("\n-> if the disulfide is NOT near these, you rigidified the wrong place.")
 
 
if __name__ == "__main__":
    main()