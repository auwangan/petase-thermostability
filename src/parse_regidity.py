import xml.etree.ElementTree as ET
import glob, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
 
RIGID = "data/v2/rigidity"
NS = "{tram:components}"
 
def curve(path):
    """Return (energies, largest_cluster_sizes) sorted by energy ascending."""
    tree = ET.parse(path)
    root = tree.getroot()
    states = root.find(f"{NS}states")
    pts = []
    for st in states.findall(f"{NS}state"):
        key = st.get("key")
        if key == "-INF":
            energy = None   # handle separately (fully folded start)
        else:
            energy = float(key)
        comps = st.find(f"{NS}components")
        if comps is None:
            sizes = [int(c.get("size")) for c in st.iter(f"{NS}component")]
        else:
            sizes = [int(c.get("size")) for c in comps.findall(f"{NS}component")]
        largest = max(sizes) if sizes else 0
        pts.append((energy, largest))
    # -INF is the folded reference -> put it at the most-negative real energy - 0.5
    real = [e for e, _ in pts if e is not None]
    if real:
        floor = min(real) - 0.5
        pts = [(floor if e is None else e, s) for e, s in pts]
    pts.sort()
    return [e for e, _ in pts], [s for _, s in pts]
 
def melting_point(energies, sizes):
    """Energy at which largest cluster drops below 50% of its max (the collapse)."""
    smax = max(sizes)
    half = smax / 2
    for e, s in zip(energies, sizes):
        if s < half:
            return e
    return energies[-1]
 
fig, ax = plt.subplots(figsize=(9,5))
results = {}
for tag, col in (("wt","steelblue"),("fast","crimson"),("mut","seagreen")):
    f = f"{RIGID}/{tag}_components.xml"
    if not os.path.exists(f):
        print(f"MISSING: {f}"); continue
    e, s = curve(f)
    mp = melting_point(e, s)
    results[tag] = mp
    T = -20*mp + 300
    ax.plot(e, s, color=col, lw=1.8, label=f"{tag} (collapse @ {mp:.2f} kcal/mol ~ {T:.0f}K)")
    ax.axvline(mp, color=col, ls="--", lw=0.8, alpha=0.6)
    print(f"{tag:5s}: largest cluster {max(s)} atoms, collapses at E={mp:.3f} kcal/mol  (~{T:.0f} K)")
 
ax.set_xlabel("H-bond energy cutoff (kcal/mol)  [more negative = higher temperature ->]")
ax.set_ylabel("largest rigid cluster (atoms)")
ax.set_title("Rigidity dilution / computational melting")
ax.legend()
ax.invert_xaxis()   # so heating goes left->right
fig.tight_layout()
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/rigidity_melting.png", dpi=150)
print(f"\nwrote figures/rigidity_melting.png")
 
print("\n=== VERDICT ===")
if "wt" in results and "fast" in results:
    dw, df = results["wt"], results["fast"]
    print(f"WT collapses at   {dw:.3f} kcal/mol")
    print(f"FAST collapses at {df:.3f} kcal/mol")
    if df < dw - 0.05:
        print("-> FAST holds rigidity to MORE NEGATIVE energy = MORE STABLE. "
              "METHOD DETECTS THE +18C DIFFERENCE. Instrument validated.")
    elif df > dw + 0.05:
        print("-> FAST collapses EARLIER than WT = wrong direction. Method fails on these structures.")
    else:
        print("-> WT and FAST collapse at ~same energy. Method cannot separate them (like MD).")
if "mut" in results and "wt" in results:
    print(f"\nmut collapses at  {results['mut']:.3f} kcal/mol (vs WT {results['wt']:.3f})")
    print("  (only meaningful if the WT/FAST validation above passed)")