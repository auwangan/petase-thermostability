import sys
 
PDB = sys.argv[1] if len(sys.argv) > 1 else (
    "data/alphafold/ispetase_G79C_R153C/ispetase_G79C_R153C_af61a/"
    "ispetase_G79C_R153C_af61a_relaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb")
 
first = last = None
cys = []
with open(PDB) as fh:
    for line in fh:
        if not line.startswith("ATOM"):
            continue
        resseq = int(line[22:26])
        resn = line[17:20].strip()
        if first is None:
            first = resseq
        last = resseq
        if resn == "CYS" and line[12:16].strip() == "CA":
            cys.append(resseq)
 
print(f"residue range: {first} .. {last}")
print(f"cysteine positions: {cys}")
print()
if first == 1:
    print("-> ColabFold numbering (starts at 1). OFFSET = +28 to reach PDB numbering.")
    print(f"   engineered Cys expected at 51,125 ; natives at 175,211,245,261")
elif first == 29:
    print("-> already in PDB (29-293) numbering. no offset needed.")
else:
    print(f"-> unexpected start ({first}); inspect before proceeding.")
print("\nexpect SIX cysteines total (2 engineered + 4 native). got", len(cys))