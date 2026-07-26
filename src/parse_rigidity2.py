import xml.etree.ElementTree as ET
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
 
RIGID = "data/v2/rigidity"
NS = "{tram:components}"
 
def q(tag): return f"{NS}{tag}"
 
def parse(path):
    root = ET.parse(path).getroot()
    states_el = root.find(q("states"))
    states = states_el.findall(q("state"))
 
    # --- step 1: sizes from the -INF (definition) state ---
    id_size = {}
    first = states[0]
    # in the -INF state, components are defined with size attr + node list
    for comp in first.iter(q("component")):
        cid = comp.get("id")
        size_attr = comp.get("size")
        if cid is None:
            continue
        if size_attr is not None:
            id_size[cid] = int(size_attr)
        else:
            # count nodes as fallback
            nodes = comp.find(q("nodes"))
            if nodes is not None:
                id_size[cid] = len(nodes.findall(q("node")))
 
    curve = []
    for st in states:
        key = st.get("key")
        energy = None if key == "-INF" else float(key)
 
        # for each top-level <component> group in this state, sum sizes of the
        # cluster-ids it references; track the largest group
        largest = 0
        top_comps = [c for c in st.findall(q("component"))]
        if not top_comps:
            # -INF state style: components directly sized
            sizes = [int(c.get("size")) for c in st.iter(q("component"))
                     if c.get("size") is not None]
            largest = max(sizes) if sizes else 0
        else:
            for grp in top_comps:
                inner = grp.find(q("components"))
                ids = []
                if inner is not None:
                    ids = [c.get("id") for c in inner.findall(q("component"))
                           if c.get("id") is not None]
                grp_size = sum(id_size.get(i, 0) for i in ids)
                # also add halo atoms (they're part of this rigid body's boundary)
                halo = grp.find(q("halo"))
                if halo is not None:
                    grp_size += len(halo.findall(q("node")))
                largest = max(largest, grp_size)
        curve.append((energy, largest))
 
    # place -INF at (min_real_energy - 0.5) so it plots at the cold end
    reals = [e for e, _ in curve if e is not None]
    if reals:
        floor = min(reals) - 0.5
        curve = [(floor if e is None else e, s) for e, s in curve]
    curve.sort()
    return [e for e, _ in curve], [s for _, s in curve]
 
def collapse_energy(E, S):
    smax = max(S) if S else 0
    for e, s in zip(E, S):
        if s < 0.5 * smax:
            return e
    return E[-1] if E else 0
 
fig, ax = plt.subplots(figsize=(9,5))
res = {}
for tag, col in (("wt","steelblue"),("fast","crimson"),("mut","seagreen")):
    p = f"{RIGID}/{tag}_components.xml"
    if not os.path.exists(p):
        print(f"MISSING {p}"); continue
    E, S = parse(p)
    ce = collapse_energy(E, S)
    T = -20*ce + 300
    res[tag] = (ce, max(S))
    ax.plot(E, S, color=col, lw=1.8, label=f"{tag}: collapse @ {ce:.2f} (~{T:.0f}K), max {max(S)}")
    ax.axvline(ce, color=col, ls="--", lw=0.7, alpha=0.6)
    print(f"{tag:5s}: max rigid cluster {max(S):5d} atoms | collapse E={ce:.3f} kcal/mol (~{T:.0f}K)")
 
ax.set_xlabel("H-bond energy cutoff (kcal/mol)")
ax.set_ylabel("largest rigid cluster (atoms)")
ax.set_title("Rigidity dilution: computational melting (leftward = hotter)")
ax.invert_xaxis()
ax.legend(fontsize=8)
fig.tight_layout()
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/rigidity_melting.png", dpi=150)
print("\nwrote figures/rigidity_melting.png")
 
print("\n=== VERDICT ===")
if "wt" in res and "fast" in res:
    we, fe = res["wt"][0], res["fast"][0]
    print(f"WT   collapse: {we:.3f} kcal/mol")
    print(f"FAST collapse: {fe:.3f} kcal/mol")
    if fe < we - 0.05:
        print("-> FAST stays rigid to more negative energy = MORE STABLE. "
              "METHOD DETECTS +18C. Instrument validated.")
    elif fe > we + 0.05:
        print("-> FAST collapses earlier (wrong direction). Method fails here.")
    else:
        print("-> WT/FAST collapse ~equal. Method can't separate them on these structures.")
if "mut" in res and "wt" in res:
    print(f"\nmut  collapse: {res['mut'][0]:.3f} kcal/mol (WT {res['wt'][0]:.3f})")
    print("  interpret only if WT/FAST validation passed above")