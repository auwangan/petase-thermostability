#!/usr/bin/env python3
"""
delta_rmsf.py -- Stockinger et al. 2025 Eq. 1

    dRMSF = mean(RMSF over 5 reps at 50 C) - mean(RMSF over 5 reps at 30 C)

Reads the per-replicate *_rmsf.csv files written by analyse_v2.py.

CRITICAL: WT (5XJH, res 30-292) and FAST (7SH6, res 29-289) cover different
ranges, so residues are matched by NUMBER, not array index.

Expected signature (their Fig. 3): in IS1-IS6, WT dRMSF is POSITIVE
(heat -> more flexible) while FAST is NEGATIVE (heat -> more rigid).
IS7 was the exception, only slightly reduced vs WT.

Usage:
  python delta_rmsf.py --analysis data/v3/analysis --out figures
"""
import argparse, glob, os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IS_REGIONS = {"IS1": (54, 70),  "IS2": (78, 102), "IS3": (110, 134),
              "IS4": (150, 182),"IS5": (186, 198),"IS6": (206, 222),
              "IS7": (230, 246)}
FAST_MUTS = {121: "S121E", 186: "D186H", 224: "R224Q", 233: "N233K", 280: "R280A"}
# Stockinger excluded termini 30-50 and 290+ from the dRMSF plots
PLOT_LO, PLOT_HI = 51, 289


def load_rmsf(path):
    """-> dict {resid: rmsf}"""
    d = {}
    with open(path) as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            d[int(row[0])] = float(row[1])
    return d


def mean_over_reps(analysis, tag, temp):
    files = sorted(glob.glob(f"{analysis}/{tag}_{temp}K_rep*_rmsf.csv"))
    if not files:
        raise SystemExit(f"no RMSF csv for {tag} at {temp}K in {analysis}")
    dicts = [load_rmsf(f) for f in files]
    common = set(dicts[0])
    for d in dicts[1:]:
        common &= set(d)
    out = {r: float(np.mean([d[r] for d in dicts])) for r in common}
    sd  = {r: float(np.std ([d[r] for d in dicts])) for r in common}
    return out, sd, len(files)


def delta(analysis, tag, t_lo, t_hi):
    lo, lo_sd, n_lo = mean_over_reps(analysis, tag, t_lo)
    hi, hi_sd, n_hi = mean_over_reps(analysis, tag, t_hi)
    res = sorted(set(lo) & set(hi))
    d   = np.array([hi[r] - lo[r] for r in res])
    err = np.array([np.hypot(hi_sd[r], lo_sd[r]) for r in res])
    print(f"  {tag}: {n_lo} reps @{t_lo}K, {n_hi} reps @{t_hi}K, {len(res)} shared residues")
    return np.array(res), d, err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis", default="data/v3/analysis")
    ap.add_argument("--t-lo", default="303")
    ap.add_argument("--t-hi", default="323")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    print("computing dRMSF = mean(RMSF_hi) - mean(RMSF_lo)")
    curves = {}
    for tag in ("wt", "fast", "mut"):
        try:
            curves[tag] = delta(a.analysis, tag, a.t_lo, a.t_hi)
        except SystemExit as e:
            print(f"  (skip {tag}: {e})")

    fig, ax = plt.subplots(figsize=(13, 4.5))
    ax.axhline(0, color="k", lw=0.7)
    for name, (lo_r, hi_r) in IS_REGIONS.items():
        ax.axvspan(lo_r, hi_r, color="grey", alpha=0.13)
        ax.text((lo_r+hi_r)/2, 0.9, name, ha="center", fontsize=7,
                transform=ax.get_xaxis_transform(), color="dimgrey")
    for tag, col in (("wt","steelblue"), ("fast","crimson"), ("mut","seagreen")):
        if tag not in curves: continue
        res, d, err = curves[tag]
        m = (res >= PLOT_LO) & (res <= PLOT_HI)
        ax.plot(res[m], d[m], color=col, lw=1.4, label=tag)
        ax.fill_between(res[m], (d-err)[m], (d+err)[m], color=col, alpha=0.18)
    for pos, lab in FAST_MUTS.items():
        ax.axvline(pos, color="k", ls=":", lw=0.6)
    ax.set_xlabel("residue"); ax.set_ylabel(r"$\Delta$RMSF  (nm)")
    ax.set_title(r"$\Delta$RMSF = mean(RMSF$_{50^\circ C}$) $-$ mean(RMSF$_{30^\circ C}$)"
                 "\npositive = heat increases flexibility;  negative = heat rigidifies")
    ax.set_xlim(PLOT_LO, PLOT_HI); ax.legend()
    fig.tight_layout(); fig.savefig(f"{a.out}/delta_rmsf.png", dpi=150)
    print(f"\nwrote {a.out}/delta_rmsf.png")

    # --- per-IS-region summary: the actual test ---
    print("\n=== mean dRMSF per instability region ===")
    hdr = f"{'region':6s} " + " ".join(f"{t:>18s}" for t in curves)
    print(hdr); print("-" * len(hdr))
    for name, (lo_r, hi_r) in IS_REGIONS.items():
        cells = []
        for tag in curves:
            res, d, err = curves[tag]
            m = (res >= lo_r) & (res <= hi_r)
            cells.append(f"{d[m].mean():+.4f}+/-{err[m].mean():.4f}" if m.any() else "    n/a")
        print(f"{name:6s} " + " ".join(f"{c:>18s}" for c in cells))

    if "wt" in curves and "fast" in curves:
        print("\n=== VALIDATION ===")
        ok = 0
        for name, (lo_r, hi_r) in IS_REGIONS.items():
            rw, dw, _ = curves["wt"]; rf, df, _ = curves["fast"]
            mw = (rw >= lo_r) & (rw <= hi_r); mf = (rf >= lo_r) & (rf <= hi_r)
            w, f = dw[mw].mean(), df[mf].mean()
            sign_flip = (w > 0) and (f < 0)
            ok += sign_flip
            print(f"  {name}: WT {w:+.4f} | FAST {f:+.4f} | "
                  f"{'SIGN FLIP (as published)' if sign_flip else 'no flip'}")
        print(f"\n  {ok}/7 regions reproduce the published WT-positive / FAST-negative pattern")
        print("  (Stockinger: IS1-IS6 flip, IS7 does not -> 6/7 expected)")


if __name__ == "__main__":
    main()