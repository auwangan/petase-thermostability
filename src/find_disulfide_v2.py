import argparse
import csv
import math
 
# ----------------------------------------------------------------------
# CONFIG -- EDIT THESE
# ----------------------------------------------------------------------

# Native disulfide Cys + already-published engineered Cys -> never touch.
#   203,239 / 273,289 : native IsPETase disulfides
#   233,282           : Zhong-Johnson 2021 engineered disulfide (+14C) -> NOT NOVEL
#   241               : HotPETase third introduced Cys
#   233 also carries your N233K anchor, so it must stay non-Cys regardless.
BLACKLIST = {203, 233, 239, 241, 273, 282, 289}
 
# Catalytic triad (PDB numbering). Pairs with either residue closer than
# TRIAD_EXCLUSION_ANGSTROM to a triad Cb are dropped as activity-risky.
TRIAD = {160, 206, 237}
TRIAD_EXCLUSION_ANGSTROM = 10.0
 
# Instability regions IS1-IS7 (Stockinger et al., CSBJ 2025, Fig. 3a).
# Read off the plot to nearest tick (~4 res) -> APPROXIMATE, fine for targeting.
IS_REGIONS = {
    "IS1": range(54, 71),
    "IS2": range(78, 103),
    "IS3": range(110, 135),
    "IS4": range(150, 183),   # LARGEST peak; contains catalytic S160 -> activity risk
    "IS5": range(186, 199),
    "IS6": range(206, 223),   # contains D206
    "IS7": range(230, 247),   # N233 anchor + published disulfide -> avoid
}
# Regions we will accept a disulfide in (novel + safe-ish). IS3/IS7 deprioritised.
TARGET_IS = {"IS1", "IS2", "IS4", "IS5", "IS6"}
 
# Geometry windows for a modellable disulfide (Angstrom).
CA_CA_MIN, CA_CA_MAX = 4.5, 7.0
CB_CB_MIN, CB_CB_MAX = 3.0, 4.5
MIN_SEQ_SEP = 3           # partners must be >= this many residues apart in sequence
 
# ----------------------------------------------------------------------
 
 
def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def scale(a, s): return (a[0]*s, a[1]*s, a[2]*s)
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
 
 
def virtual_cb(n, ca, c):
    """Idealised Cb from backbone N, CA, C (for GLY or missing Cb). Approximate."""
    b = sub(ca, n)
    c_ = sub(c, ca)
    a = cross(b, c_)
    return add(scale(a, -0.58273431),
               add(scale(b, 0.56802827),
                   add(scale(c_, -0.54067466), ca)))
 
 
def parse_pdb(path, chain_pref=None):
    """Return {resSeq: (resname, CA, CB)} from ATOM records. Picks one chain."""
    backbone = {}           # (chain,resSeq) -> {atomname: xyz}
    resname = {}            # (chain,resSeq) -> name
    chains_seen = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            name = line[12:16].strip()
            if name not in ("N", "CA", "C", "CB"):
                continue
            chain = line[21]
            if chain not in chains_seen:
                chains_seen.append(chain)
            try:
                resseq = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            key = (chain, resseq)
            backbone.setdefault(key, {})[name] = (x, y, z)
            resname[key] = line[17:20].strip()
 
    chain = chain_pref or ("A" if "A" in chains_seen else chains_seen[0])
    if len(chains_seen) > 1:
        print(f"[note] chains present {chains_seen}; using chain '{chain}'. "
              f"Override with --chain.")
 
    out = {}
    for (ch, rs), atoms in backbone.items():
        if ch != chain or "CA" not in atoms:
            continue
        ca = atoms["CA"]
        if "CB" in atoms:
            cb = atoms["CB"]
        elif "N" in atoms and "C" in atoms:
            cb = virtual_cb(atoms["N"], ca, atoms["C"])
        else:
            continue
        out[rs] = (resname[(ch, rs)], ca, cb)
    return out
 
 
def is_region_of(resseq):
    for nmreg, rng in IS_REGIONS.items():
        if resseq in rng:
            return nmreg
    return None
 
 
def load_rmsf(path):
    if not path:
        return None
    d = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", " ").split()
            try:
                d[int(float(parts[0]))] = float(parts[1])
            except (ValueError, IndexError):
                continue
    return d
 
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("--rmsf", default=None,
                    help="per-residue RMSF file: 'resSeq rmsf' per line")
    ap.add_argument("--chain", default=None)
    ap.add_argument("--out", default="disulfide_candidates_v2.csv")
    args = ap.parse_args()
 
    residues = parse_pdb(args.pdb, args.chain)
    rmsf = load_rmsf(args.rmsf)
    print(f"Parsed {len(residues)} residues. RMSF file: {'yes' if rmsf else 'no (geometry+IS only)'}")
 
    triad_cb = [residues[r][2] for r in TRIAD if r in residues]
    resseqs = sorted(residues)
    hits = []
 
    for i, ri in enumerate(resseqs):
        if ri in BLACKLIST:
            continue
        nm_i, ca_i, cb_i = residues[ri]
        for rj in resseqs[i + 1:]:
            if rj in BLACKLIST or abs(ri - rj) < MIN_SEQ_SEP:
                continue
            nm_j, ca_j, cb_j = residues[rj]
 
            d_caca = math.dist(ca_i, ca_j)
            if not (CA_CA_MIN <= d_caca <= CA_CA_MAX):
                continue
            d_cbcb = math.dist(cb_i, cb_j)
            if not (CB_CB_MIN <= d_cbcb <= CB_CB_MAX):
                continue
 
            if triad_cb:
                min_triad = min([math.dist(cb_i, t) for t in triad_cb] +
                                [math.dist(cb_j, t) for t in triad_cb])
            else:
                min_triad = 999.0
            if min_triad < TRIAD_EXCLUSION_ANGSTROM:
                continue
 
            is_i, is_j = is_region_of(ri), is_region_of(rj)
            if (230 <= ri <= 246) or (230 <= rj <= 246):
                continue
            if not ((is_i in TARGET_IS) or (is_j in TARGET_IS)):
                continue
 
            score, note = 0.0, []
            if rmsf and ri in rmsf and rj in rmsf:
                diff = abs(rmsf[ri] - rmsf[rj])
                score += diff
                note.append(f"dRMSF={diff:.2f}")
            score += max(0.0, 1.0 - abs(d_cbcb - 3.8))
            score += min(min_triad, 25.0) / 25.0
 
            hits.append({
                "mutation": f"{nm_i}{ri}C/{nm_j}{rj}C",
                "resi": ri, "aa_i": nm_i, "IS_i": is_i or "-",
                "resj": rj, "aa_j": nm_j, "IS_j": is_j or "-",
                "CaCa": round(d_caca, 2), "CbCb": round(d_cbcb, 2),
                "min_triad_dist": round(min_triad, 1),
                "score": round(score, 3), "notes": ";".join(note),
            })
            
    hits.sort(key=lambda h: h["score"], reverse=True)
    cols = ["mutation", "resi", "aa_i", "IS_i", "resj", "aa_j", "IS_j",
            "CaCa", "CbCb", "min_triad_dist", "score", "notes"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for h in hits:
            w.writerow(h)
 
    print(f"\n{len(hits)} candidate pairs (target IS regions, blacklist + activity filter).")
    print("Top hits (VALIDATE chi-angles in Disulfide by Design 2 / Rosetta before designing):")
    for h in hits[:15]:
        print(f"  {h['mutation']:16s} IS={h['IS_i']}/{h['IS_j']:4s} "
              f"CaCa={h['CaCa']} CbCb={h['CbCb']} triad={h['min_triad_dist']}A "
              f"score={h['score']} {h['notes']}")
    if not hits:
        print("  EMPTY -> widen CA_CA_MAX/CB_CB_MAX a little, re-check IS ranges, "
              "or drop to the salt-bridge fallback.")
 
 
if __name__ == "__main__":
    main()