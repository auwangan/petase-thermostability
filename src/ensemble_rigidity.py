import argparse, os, glob, subprocess
import numpy as np
import xml.etree.ElementTree as ET
 
NS = "{tram:components}"
ENS = "data/v2/ensemble"
THRESHOLDS = [0.75, 0.50, 0.25]
GRID = np.arange(0.0, -6.0, -0.02)   # common energy axis
 
# ---------- stage 1: extract ----------
def extract(tag, per_rep, tail_frac):
    import mdtraj as md
    top = f"data/v2/structures/{tag}_v2_top.pdb"
    trajs = sorted(glob.glob(f"data/md/v2_prod_{tag}_350K_rep[1-5].dcd"))
    os.makedirs(ENS, exist_ok=True)
    n = 0
    for ti, tp in enumerate(trajs, 1):
        t = md.load(tp, top=top)
        start = int(t.n_frames * (1 - tail_frac))          # index-based, not time
        idxs = np.linspace(start, t.n_frames - 1, per_rep).astype(int)
        for k, fi in enumerate(idxs):
            out = f"{ENS}/{tag}_r{ti}_f{k}.pdb"
            t[int(fi)].save_pdb(out)
            print(f"  {out}  (frame {fi}/{t.n_frames})")
            n += 1
    print(f"[{tag}] {n} frames")
 
# ---------- stage 2: dilution ----------
def dilute(tag):
    frames = sorted(glob.glob(f"{ENS}/{tag}_r*_f*.pdb"))
    for i, f in enumerate(frames, 1):
        name = os.path.basename(f).replace(".pdb", "")
        if os.path.exists(f"{ENS}/{name}_components.xml"):
            print(f"  [{i}/{len(frames)}] {name} exists, skip"); continue
        print(f"  [{i}/{len(frames)}] {name} ...", flush=True)
        subprocess.run(["tram-pdb", "-p", f, "-o", ENS, "-n", name,
                        "-t", "-0.1", "-l", "ERROR"], check=False)
 
# ---------- P_inf from one components.xml ----------
def pinf_curve(path):
    root = ET.parse(path).getroot()
    id_size = {int(c.get("id")): int(c.get("size"))
               for c in root.find(f"{NS}components").findall(f"{NS}component")}
    total = sum(id_size.values())
    pts = []
    for st in root.find(f"{NS}states").findall(f"{NS}state"):
        k = st.get("key")
        if k == "-INF":
            continue
        E = float(k)
        largest = 0
        for body in st.findall(f"{NS}component"):
            inner = body.find(f"{NS}components")
            ids = [int(x.get("id")) for x in inner.findall(f"{NS}component")] if inner is not None else []
            s = sum(id_size.get(i, 0) for i in ids)
            halo = body.find(f"{NS}halo")
            if halo is not None:
                s += len(halo.findall(f"{NS}node"))
            largest = max(largest, s)
        pts.append((E, largest / total))
    pts.sort(key=lambda x: -x[0])           # folded end first
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
 
def on_grid(E, P):
    """interpolate P_inf onto the common (descending) energy grid"""
    order = np.argsort(E)                    # ascending for np.interp
    return np.interp(GRID, E[order], P[order])
 
# ---------- stage 3: analyse ----------
def analyse():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 5.5))
    summary = {}
    for tag, col in (("wt", "steelblue"), ("fast", "crimson"), ("mut", "seagreen")):
        files = sorted(glob.glob(f"{ENS}/{tag}_r*_f*_components.xml"))
        if not files:
            print(f"no ensemble files for {tag}"); continue
        curves = []
        for f in files:
            try:
                E, P = pinf_curve(f)
                curves.append(on_grid(E, P))
            except Exception as e:
                print(f"  skip {os.path.basename(f)}: {e}")
        M = np.vstack(curves)
        mean, sd = M.mean(0), M.std(0)
        ax.plot(GRID, mean, color=col, lw=2, label=f"{tag} (n={len(curves)})")
        ax.fill_between(GRID, mean - sd, mean + sd, color=col, alpha=0.2)
        # threshold crossings on the MEAN curve, plus per-frame scatter
        row = {}
        for th in THRESHOLDS:
            idx = np.argmax(mean < th) if (mean < th).any() else -1
            e_mean = GRID[idx] if idx >= 0 else np.nan
            per_frame = []
            for c in M:
                j = np.argmax(c < th) if (c < th).any() else -1
                per_frame.append(GRID[j] if j >= 0 else np.nan)
            row[th] = (e_mean, np.nanmean(per_frame), np.nanstd(per_frame))
        summary[tag] = row
        print(f"\n{tag} (n={len(curves)} frames):")
        for th in THRESHOLDS:
            em, pm, ps = row[th]
            print(f"  P<{th:.2f}: mean-curve E={em:6.2f} | per-frame {pm:6.2f} +/- {ps:.2f} kcal/mol"
                  f"  (~{-20*pm+300:.0f} K)")
    ax.set_xlabel("E_cut,hb (kcal/mol)"); ax.set_ylabel(r"$P_\infty$")
    ax.set_title("Ensemble rigidity (mean +/- SD over MD frames)")
    for th in THRESHOLDS: ax.axhline(th, color="k", lw=0.4, ls=":")
    ax.legend(); fig.tight_layout()
    os.makedirs("figures", exist_ok=True)
    fig.savefig("figures/rigidity_ensemble.png", dpi=150)
    print("\nwrote figures/rigidity_ensemble.png")
 
    print("\n=== VALIDATION: FAST vs WT (experiment: +18 K) ===")
    if "wt" in summary and "fast" in summary:
        for th in THRESHOLDS:
            _, wm, ws = summary["wt"][th]
            _, fm, fs = summary["fast"][th]
            dT = (-20*fm) - (-20*wm)
            sep = abs(fm - wm) / max(1e-9, np.hypot(ws, fs))
            print(f"  P<{th:.2f}: WT {wm:.2f}+/-{ws:.2f} | FAST {fm:.2f}+/-{fs:.2f} | "
                  f"shift {dT:+.0f} K | separation {sep:.1f} sigma")
    if "wt" in summary and "mut" in summary:
        print("\n=== mutant vs WT ===")
        for th in THRESHOLDS:
            _, wm, ws = summary["wt"][th]
            _, mm, ms = summary["mut"][th]
            dT = (-20*mm) - (-20*wm)
            sep = abs(mm - wm) / max(1e-9, np.hypot(ws, ms))
            print(f"  P<{th:.2f}: WT {wm:.2f}+/-{ws:.2f} | mut {mm:.2f}+/-{ms:.2f} | "
                  f"shift {dT:+.0f} K | separation {sep:.1f} sigma")
 
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["extract", "dilute", "analyse"])
    ap.add_argument("--tag", choices=["wt", "mut", "fast"])
    ap.add_argument("--per-rep", type=int, default=2)
    ap.add_argument("--tail-frac", type=float, default=0.3, help="sample from last 30% of each rep")
    a = ap.parse_args()
    if a.stage == "extract": extract(a.tag, a.per_rep, a.tail_frac)
    elif a.stage == "dilute": dilute(a.tag)
    else: analyse()