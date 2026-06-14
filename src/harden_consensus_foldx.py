# src/harden_consensus_foldx.py
# Block 1: re-run FoldX 3x on the consensus hits and average,
# to separate real consensus stabilizers from single-run noise.
import os, subprocess
import pandas as pd

# 1. consensus hits = the ones FoldX also called stabilizing yesterday
v3 = pd.read_csv("data/foldx_vs_thermompnn_v3.csv")
consensus = v3[v3["foldx_ddg"] < 0].copy()
print(f"consensus hits to harden: {len(consensus)}")

# 2. build the FoldX mutation list (same format as run_foldx_batch.py)
foldx_mutations = [f"{r.wildtype}A{int(r.position)}{r.mutation};" for r in consensus.itertuples()]
print(foldx_mutations)
with open("data/individual_list.txt", "w") as f:
    f.write("\n".join(foldx_mutations) + "\n")

# 3. clear old outputs so we don't read stale data
for fp in ["data/Dif_6eqe_Repair.fxout", "data/Average_6eqe_Repair.fxout"]:
    if os.path.exists(fp):
        os.remove(fp)

# 4. run FoldX — the ONLY change vs your batch script is the last line
subprocess.run([
    "foldx",
    "--command=BuildModel",
    "--pdb=6eqe_Repair.pdb",
    "--pdb-dir=data/",
    "--mutant-file=data/individual_list.txt",
    "--output-dir=data/",
    "--numberOfRuns=3",
])

# 5. INSPECT the averaged output before trusting any parse
for name in ["Average_6eqe_Repair.fxout", "Dif_6eqe_Repair.fxout"]:
    print(f"\n=== {name} (first 20 lines) ===")
    with open(f"data/{name}") as f:
        for i, line in enumerate(f):
            if i >= 20: break
            print(line.rstrip())