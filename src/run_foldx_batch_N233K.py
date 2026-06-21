import pandas as pd
import subprocess
#Load the candidate mutations from data/phase1_candidates.csv
candidates_df = pd.read_csv("data/n233k_repro.csv")
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
print(foldx_df[["total energy"]])
