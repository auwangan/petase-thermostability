import xml.etree.ElementTree as ET
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NS = "{tram:components}"
RIGID = "data/v2/rigidity"

def pinf_curve(path):
    root = ET.parse(path).getroot()
    # 1) dictionary: id -> size
    id_size = {int(c.get("id")): int(c.get("size"))
               for c in root.find(f"{NS}components").findall(f"{NS}component")}
    total = sum(id_size.values())

    # 2) each state: largest rigid body = max over grouped-id-sum
    out = []
    for st in root.find(f"{NS}states").findall(f"{NS}state"):
        key = st.get("key")
        E = None if key == "-INF" else float(key)
        largest = 0
        for body in st.findall(f"{NS}component"):
            inner = body.find(f"{NS}components")
            ids = [int(x.get("id")) for x in inner.findall(f"{NS}component")] if inner is not None else []
            s = sum(id_size.get(i, 0) for i in ids)
            # add halo atoms of this body
            halo = body.find(f"{NS}halo")
            if halo is not None:
                s += len(halo.findall(f"{NS}node"))
            largest = max(largest, s)
        out.append((E, largest / total if total else 0))
    # place -INF at cold end
    reals = [e for e, _ in out if e is not None]
    floor = (min(reals) - 0.5) if reals else 0
    out = [((floor if e is None else e), p) for e, p in out]
    out.sort()
    return [e for e, _ in out], [p for _, p in out], total

def transition_E(E, P):
    """Late transition: E where P_inf drops below 0.5 (giant cluster loses majority)."""
    for e, p in zip(E, P):
        if p < 0.5:
            return e
    return E[-1]

fig, ax = plt.subplots(figsize=(9,5))
res = {}
for tag, col in (("wt","steelblue"),("fast","crimson"),("mut","seagreen")):
    p = f"{RIGID}/{tag}_components.xml"
    if not os.path.exists(p):
        print(f"MISSING {p}"); continue
    E, P, total = pinf_curve(p)
    Et = transition_E(E, P)
    T = -20*Et + 300
    res[tag] = (Et, T, max(P), total)
    ax.plot(E, P, color=col, lw=1.9, marker='.', ms=3,
            label=f"{tag}: transition E={Et:.2f} (~{T:.0f}K), Pmax={max(P):.2f}")
    ax.axvline(Et, color=col, ls="--", lw=0.7, alpha=0.5)
    print(f"{tag:5s}: total {total} atoms | Pinf_max {max(P):.3f} | "
          f"transition E={Et:.3f} kcal/mol (~{T:.0f}K)")

ax.axhline(0.5, color="k", lw=0.5, ls=":")
ax.set_xlabel("H-bond energy cutoff E_cut,hb (kcal/mol)")
ax.set_ylabel(r"$P_\infty$  (giant rigid cluster fraction)")
ax.set_title("Rigidity order parameter vs dilution (Radestock-Gohlke)\nleftward = higher temperature")
ax.invert_xaxis()
ax.legend(fontsize=8)
fig.tight_layout()
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/rigidity_pinf.png", dpi=150)
print("\nwrote figures/rigidity_pinf.png")

print("\n=== VERDICT (validation on FAST positive control) ===")
if "wt" in res and "fast" in res:
    ew, tw = res["wt"][0], res["wt"][1]
    ef, tf = res["fast"][0], res["fast"][1]
    print(f"WT   transition: E={ew:.3f} kcal/mol  (~{tw:.0f} K)")
    print(f"FAST transition: E={ef:.3f} kcal/mol  (~{tf:.0f} K)")
    print(f"shift: {tf-tw:+.0f} K  (experiment: FAST is +18 K)")
    if ef < ew - 0.1:
        print("-> FAST transition shifted to higher T = MORE STABLE. "
              "Rigidity analysis DETECTS the effect. (Radestock-Gohlke signature)")
    elif ef > ew + 0.1:
        print("-> FAST transition LOWER than WT (wrong direction).")
    else:
        print("-> WT/FAST transitions ~equal. Single-structure rigidity cannot "
              "separate them; literature says use ENSEMBLE (tram-xtc on trajectories).")
if "mut" in res and "wt" in res:
    print(f"\nmut transition: E={res['mut'][0]:.3f} (~{res['mut'][1]:.0f}K) vs WT {res['wt'][1]:.0f}K")
    print("  (interpret only if FAST validation passed)")