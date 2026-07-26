import xml.etree.ElementTree as ET
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NS = "{tram:components}"
RIGID = "data/v2/rigidity"
THRESHOLDS = [0.75, 0.50, 0.25]

def pinf_curve(path):
    root = ET.parse(path).getroot()
    id_size = {int(c.get("id")): int(c.get("size"))
               for c in root.find(f"{NS}components").findall(f"{NS}component")}
    total = sum(id_size.values())
    pts = []
    for st in root.find(f"{NS}states").findall(f"{NS}state"):
        key = st.get("key")
        E = None if key == "-INF" else float(key)
        largest = 0
        for body in st.findall(f"{NS}component"):
            inner = body.find(f"{NS}components")
            ids = [int(x.get("id")) for x in inner.findall(f"{NS}component")] if inner is not None else []
            s = sum(id_size.get(i, 0) for i in ids)
            halo = body.find(f"{NS}halo")
            if halo is not None:
                s += len(halo.findall(f"{NS}node"))
            largest = max(largest, s)
        pts.append((E, largest/total if total else 0))
    reals = [e for e, _ in pts if e is not None]
    floor = (min(reals) - 0.5) if reals else 0
    pts = [((floor if e is None else e), p) for e, p in pts]
    # sort DESCENDING by energy: start at folded end (E near 0) -> unfolding
    pts.sort(key=lambda x: -x[0])
    return [e for e, _ in pts], [p for _, p in pts], total

def crossings(E, P, thresholds):
    """First energy (scanning from folded end) where P_inf drops below each threshold."""
    out = {}
    for th in thresholds:
        hit = None
        for e, p in zip(E, P):
            if p < th:
                hit = e
                break
        out[th] = hit
    return out

fig, ax = plt.subplots(figsize=(10,5.5))
res = {}
for tag, col in (("wt","steelblue"),("fast","crimson"),("mut","seagreen")):
    f = f"{RIGID}/{tag}_components.xml"
    if not os.path.exists(f):
        print(f"MISSING {f}"); continue
    E, P, total = pinf_curve(f)
    cr = crossings(E, P, THRESHOLDS)
    res[tag] = cr
    ax.plot(E, P, color=col, lw=1.9, label=tag)
    e50 = cr[0.50]
    if e50 is not None:
        ax.plot([e50],[0.5], 'o', color=col, ms=7)
    print(f"{tag:5s} (total {total} atoms, Pmax {max(P):.3f}):")
    for th in THRESHOLDS:
        e = cr[th]
        T = -20*e + 300 if e is not None else float('nan')
        print(f"    P_inf < {th:.2f}  at E = {e:7.3f} kcal/mol   (~{T:.0f} K / {T-273:.0f} C)")

for th in THRESHOLDS:
    ax.axhline(th, color="k", lw=0.4, ls=":")
ax.set_xlabel("H-bond energy cutoff E_cut,hb (kcal/mol)   [rightward = more bonds removed = hotter]")
ax.set_ylabel(r"$P_\infty$  (giant rigid cluster fraction)")
ax.set_title("Rigidity order parameter vs dilution (Radestock-Gohlke)")
ax.legend()
fig.tight_layout()
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/rigidity_pinf_fixed.png", dpi=150)
print("\nwrote figures/rigidity_pinf_fixed.png")

print("\n=== VALIDATION: does rigidity detect FAST (+18 C experimentally)? ===")
if "wt" in res and "fast" in res:
    for th in THRESHOLDS:
        ew, ef = res["wt"][th], res["fast"][th]
        if ew is None or ef is None: continue
        dT = (-20*ef + 300) - (-20*ew + 300)
        verdict = "FAST more stable" if ef < ew - 0.05 else ("WT more stable" if ef > ew + 0.05 else "~equal")
        print(f"  at P_inf={th:.2f}:  WT {ew:7.3f} | FAST {ef:7.3f} | shift {dT:+.0f} K  -> {verdict}")
if "mut" in res and "wt" in res:
    print("\n=== disulfide mutant vs WT ===")
    for th in THRESHOLDS:
        ew, em = res["wt"][th], res["mut"][th]
        if ew is None or em is None: continue
        dT = (-20*em + 300) - (-20*ew + 300)
        verdict = "mut more stable" if em < ew - 0.05 else ("mut LESS stable" if em > ew + 0.05 else "~equal")
        print(f"  at P_inf={th:.2f}:  WT {ew:7.3f} | mut {em:7.3f} | shift {dT:+.0f} K  -> {verdict}")