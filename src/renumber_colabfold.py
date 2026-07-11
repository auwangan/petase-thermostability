IN  = ("data/alphafold/ispetase_G79C_R153C/ispetase_G79C_R153C_af61a/"
       "ispetase_G79C_R153C_af61a_relaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb")
OUT = "data/md/mut_v2_renumbered.pdb"
OFFSET = 28
import os
os.makedirs("data/md", exist_ok=True)
n = 0
with open(IN) as fin, open(OUT, "w") as fout:
    for line in fin:
        if line.startswith(("ATOM", "HETATM", "TER", "ANISOU")):
            resseq = int(line[22:26]) + OFFSET
            line = line[:22] + f"{resseq:>4}" + line[26:]
            n += 1
        fout.write(line)
print(f"renumbered {n} atom records (+{OFFSET}) -> {OUT}")
cys, first, last = [], None, None
with open(OUT) as fh:
    for line in fh:
        if line.startswith("ATOM"):
            rs = int(line[22:26])
            first = rs if first is None else first
            last = rs
            if line[17:20].strip() == "CYS" and line[12:16].strip() == "CA":
                cys.append(rs)
print(f"new residue range: {first}..{last}  (expect 29..293)")
print(f"cysteines now at:  {cys}  (expect 79,153,203,239,273,289)")
