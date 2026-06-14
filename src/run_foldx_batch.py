import pandas as pd
import subprocess
#Load the candidate mutations from data/phase1_candidates.csv
candidates_df = pd.read_csv("data/phase1_candidates_v3.csv")
foldx_mutations = []
import os
for fpath in ["data/Dif_6eqe_Repair.fxout"]:
    if os.path.exists(fpath):
        os.remove(fpath)
for index, row in candidates_df.iterrows():
    position = row["position"]
    wildtype = row["wildtype"]
    mutation = row["mutation"]
    mutation_str = f"{wildtype}A{int(position)}{mutation};"
    foldx_mutations.append(mutation_str)
print(foldx_mutations)
with open("data/individual_list.txt", "w") as f:
    for m in foldx_mutations:
        f.write(m + "\n")
subprocess.run([
  "foldx",
  "--command=BuildModel",
    "--pdb=6eqe_Repair.pdb",
   "--pdb-dir=data/",
   "--mutant-file=data/individual_list.txt",
   "--output-dir=data/"
])
#save it as csv
foldx_df = pd.read_csv("data/Dif_6eqe_Repair.fxout", sep="\t", skiprows=8)
print(foldx_df.columns.tolist())
print(foldx_df.head(10))
foldx_df = foldx_df.reset_index(drop=True)
print(foldx_df.shape)
foldx_df["mutation_str"] = foldx_mutations
# attach FoldX ddG alongside the candidate info (same order)
candidates_df["foldx_ddg"] = foldx_df["total energy"].values
candidates_df["foldx_mutation"] = foldx_df["mutation_str"].values

# build a clean comparison table
comparison = candidates_df[["position", "wildtype", "mutation", "ddG_pred", "foldx_ddg", "location", "rsa"]]
comparison = comparison.rename(columns={"ddG_pred": "thermompnn_ddg"})

# sort by FoldX (most stabilizing = most negative first)
comparison = comparison.sort_values("foldx_ddg")

print(comparison.to_string())
comparison.to_csv("data/foldx_vs_thermompnn_v3.csv", index=False)
print(comparison.shape)          
