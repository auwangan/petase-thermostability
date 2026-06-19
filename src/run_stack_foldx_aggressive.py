# src/run_stack_foldx.py
import os, subprocess
import pandas as pd

stack = ["NA233K", "TA77I", "AA179V", "RA260F", "QA119D"]
individual = ",".join(stack) + ";"
print("stack:", individual)

with open("data/individual_list.txt", "w") as f:
    f.write(individual + "\n")
for fp in ["data/Dif_6eqe_Repair.fxout", "data/Average_6eqe_Repair.fxout"]:
    if os.path.exists(fp): os.remove(fp)

subprocess.run([
    "foldx", "--command=BuildModel", "--pdb=6eqe_Repair.pdb",
    "--pdb-dir=data/", "--mutant-file=data/individual_list.txt",
    "--output-dir=data/", "--numberOfRuns=3",
])

avg = pd.read_csv("data/Average_6eqe_Repair.fxout", sep="\t", skiprows=8)
measured, sd = avg["total energy"].iloc[0], avg["SD"].iloc[0]

singles = pd.read_csv("data/consensus_hardened.csv")
def single_mean(m):
    wt, pos, mt = m[0], int(m[2:-1]), m[-1]
    r = singles[(singles.position==pos)&(singles.wildtype==wt)&(singles.mutation==mt)]
    return float(r["foldx_mean"].iloc[0])
additive = sum(single_mean(m) for m in stack)

print(f"\nmeasured stack ddG : {measured:.2f} ± {sd:.2f}")
print(f"additive (sum)     : {additive:.2f}")
print(f"epistasis          : {measured - additive:.2f}")