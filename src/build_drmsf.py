import csv
import argparse
 
OFFSET = 28   # resid 1  ->  PDB residue 29
 
# Terminal residues that flail freely at high T and are NOT real unfolding
# nuclei. Stockinger et al. excluded these exact ranges from their dRMSF map.
# Given in PDB numbering.
TERMINI_EXCLUDE = set(range(29, 51)) | set(range(290, 294))
 
def load(path):
    """Read an analyze_full.py rmsf csv -> {resid: rmsf}."""
    d = {}
    with open(path) as fh:
        r = csv.reader(fh)
        header = next(r)                      # skip 'resid,rmsf_A'
        for row in r:
            if len(row) < 2:
                continue
            d[int(float(row[0]))] = float(row[1])
    return d
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hot", nargs="+", required=True,
                    help="one or more 400K rmsf csvs (averaged)")
    ap.add_argument("--cold", nargs="+", required=True,
                    help="one or more 300K rmsf csvs (averaged)")
    ap.add_argument("--out", default="data/wt_drmsf.txt")
    args = ap.parse_args()
 
    hot = [load(p) for p in args.hot]
    cold = [load(p) for p in args.cold]
 
    # residues present in every file
    resids = set(hot[0])
    for d in hot[1:] + cold:
        resids &= set(d)
    resids = sorted(resids)
 
    rows = []
    skipped = 0
    for rid in resids:
        resseq = rid + OFFSET                   # <-- numbering fix
        if resseq in TERMINI_EXCLUDE:           # <-- drop flailing termini
            skipped += 1
            continue
        hot_mean = sum(d[rid] for d in hot) / len(hot)
        cold_mean = sum(d[rid] for d in cold) / len(cold)
        drmsf = hot_mean - cold_mean
        rows.append((resseq, drmsf))
 
    with open(args.out, "w") as fh:
        for resseq, drmsf in rows:
            fh.write(f"{resseq} {drmsf:.4f}\n")
 
    print(f"Wrote {len(rows)} residues to {args.out} (skipped {skipped} terminal)")
    print(f"  numbering: resid {resids[0]}..{resids[-1]}  ->  PDB {resids[0]+OFFSET}..{resids[-1]+OFFSET}")
    hi = sorted(rows, key=lambda x: x[1], reverse=True)[:8]
    print("  most heat-destabilised (highest dRMSF), PDB numbering:")
    for resseq, d in hi:
        print(f"    res {resseq}: dRMSF={d:.3f}")
 
if __name__ == "__main__":
    main()