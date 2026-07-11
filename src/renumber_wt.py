IN  = ("data/alphafold/ispetase_wt_val_924d7_2/"
       "ispetase_wt_val_924d7_2_relaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb")
OUT = "data/v2/structures/wt_v2.pdb"
OFFSET = 28
import os
os.makedirs("data/v2/structures", exist_ok=True)
n = 0
with open(IN) as fin, open(OUT, "w") as fout:
    for line in fin:
        if line.startswith(("ATOM", "HETATM", "TER", "ANISOU")):
            line = line[:22] + f"{int(line[22:26]) + OFFSET:>4}" + line[26:]
            n += 1
        fout.write(line)
print(f"renumbered {n} atoms -> {OUT}")
cys, first, last = [], None, None
with open(OUT) as fh:
    for line in fh:
        if line.startswith("ATOM"):
            rs = int(line[22:26])
            first = rs if first is None else first; last = rs
            if line[17:20].strip() == "CYS" and line[12:16].strip() == "CA":
                cys.append(rs)
print(f"range: {first}..{last}  (expect 29..293)")
print(f"cysteines: {cys}  (expect 203,239,273,289 -- FOUR, not six)")
