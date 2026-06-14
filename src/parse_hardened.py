# src/parse_hardened.py
import pandas as pd

# the 15 consensus mutations, in the SAME order FoldX received them
v3 = pd.read_csv("data/foldx_vs_thermompnn_v3.csv")
consensus = v3[v3["foldx_ddg"] < 0].reset_index(drop=True)

# Average_ file: header on line 9 -> skiprows=8; cols = Pdb, SD, total energy, ...
avg = pd.read_csv("data/Average_6eqe_Repair.fxout", sep="\t", skiprows=8).reset_index(drop=True)

consensus["foldx_single"] = consensus["foldx_ddg"]        # yesterday's 1-run value
consensus["foldx_mean"]   = avg["total energy"].values    # 3-run average
consensus["foldx_sd"]     = avg["SD"].values
consensus["upper"]        = consensus["foldx_mean"] + consensus["foldx_sd"]
consensus["robust"]       = consensus["upper"] < 0        # stabilizing even at +1 SD

out = consensus[["position","wildtype","mutation","thermompnn_ddg",
                 "foldx_single","foldx_mean","foldx_sd","robust","location","rsa"]
               ].sort_values("foldx_mean")
print(out.to_string(index=False))
out.to_csv("data/consensus_hardened.csv", index=False)
print("\nrobust consensus:", int(out["robust"].sum()), "of", len(out))