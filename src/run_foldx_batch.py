import pandas as pd
import subprocess
import os

# input: the 6 negative controls
df = pd.read_csv("data/negatives.csv")

# clear stale FoldX output so we never parse an old run
if os.path.exists("data/Dif_6eqe_Repair.fxout"):
    os.remove("data/Dif_6eqe_Repair.fxout")

# build FoldX mutation strings (chain A, PDB numbering)
foldx_mutations = []
for _, row in df.iterrows():
    foldx_mutations.append(f"{row['wildtype']}A{int(row['position'])}{row['mutation']};")
print(foldx_mutations)

with open("data/individual_list.txt", "w") as f:
    for m in foldx_mutations:
        f.write(m + "\n")

# run FoldX BuildModel at 1 run (same as controls + N233K)
subprocess.run([
    "foldx",
    "--command=BuildModel",
    "--pdb=6eqe_Repair.pdb",
    "--pdb-dir=data/",
    "--mutant-file=data/individual_list.txt",
    "--output-dir=data/",
])

# parse the ddG (Dif_ = mutant - WT), label each row
foldx_df = pd.read_csv("data/Dif_6eqe_Repair.fxout", sep="\t", skiprows=8)
foldx_df["mutation"] = foldx_mutations
print(foldx_df[["mutation", "total energy"]])     
