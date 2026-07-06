import numpy as np
import matplotlib.pyplot as plt

def load(f):
    d = np.loadtxt(f, delimiter=",", skiprows=1)
    return d[:,0], d[:,1]  # resid, rmsf

# average the 2 reps for each
def avg(f1, f2):
    r, a = load(f1); _, b = load(f2)
    return r, (a + b) / 2

r_wt, wt = avg("analysis/wt_400K_rep1_rmsf.csv", "analysis/wt_400K_rep2_rmsf.csv")
r_mut, mut = avg("analysis/mut_400K_rep1_rmsf.csv", "analysis/mut_400K_rep2_rmsf.csv")

fig, ax = plt.subplots(figsize=(14,6))
ax.plot(r_wt, wt, color="blue", label="WT (400K, 2-rep avg)", lw=1.2)
ax.plot(r_mut, mut, color="red", label="Mutant (400K, 2-rep avg)", lw=1.2)
ax.fill_between(r_mut, wt, mut, where=(mut>wt), alpha=0.2, color="red")

# mark the 5 mutations
muts = {77:"T77I", 119:"Q119D", 179:"A179V", 233:"N233K", 260:"R260F"}
for pos, name in muts.items():
    ax.axvline(pos, color="gray", ls="--", alpha=0.5)
    ax.text(pos, ax.get_ylim()[1]*0.9, name, rotation=90, fontsize=8, ha="right")

# mark catalytic triad
for pos in [160, 206, 237]:
    ax.axvline(pos, color="green", ls=":", alpha=0.4)

ax.set_xlabel("Residue"); ax.set_ylabel("RMSF (Å)")
ax.set_title("Per-residue flexibility at 400K: WT vs 5-mutant")
ax.legend(); ax.set_ylim(0, 6)  # cap y so termini don't dominate
plt.tight_layout()
plt.savefig("rmsf_comparison.png", dpi=130)
print("saved rmsf_comparison.png")

# quantify: RMSF difference near each mutation site (+-3 residues)
print("\nRMSF change (mut - wt) in +-3 residue window around each mutation:")
for pos, name in muts.items():
    mask = (r_wt >= pos-3) & (r_wt <= pos+3)
    diff = (mut[mask] - wt[mask]).mean()
    print(f"  {name}: {diff:+.2f} A  ({'MORE flexible' if diff>0 else 'less flexible'} in mutant)")