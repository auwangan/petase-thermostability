import subprocess, os, math
 
PDB = "6eqe_Repair.pdb"
PDB_DIR = "data/"
MUT_LINE = "GA79C,RA153C;"      # <-- both in ONE variant (comma-separated)
 
# ---- 1. write the mutant file (single line = single double-variant) ----
os.makedirs(PDB_DIR, exist_ok=True)
for stale in ("data/Dif_6eqe_Repair.fxout", "data/6eqe_Repair_1.pdb"):
    if os.path.exists(stale):
        os.remove(stale)
 
with open("data/individual_list.txt", "w") as f:
    f.write(MUT_LINE + "\n")
print(f"individual_list.txt -> {MUT_LINE}")
 
# ---- 2. run FoldX BuildModel ----
subprocess.run([
    "foldx",
    "--command=BuildModel",
    f"--pdb={PDB}",
    f"--pdb-dir={PDB_DIR}",
    "--mutant-file=data/individual_list.txt",
    f"--output-dir={PDB_DIR}",
], check=True)
 
# FoldX writes the mutant as <pdbbase>_1.pdb
MUT_PDB = "data/6eqe_Repair_1.pdb"
print(f"\nmutant structure: {MUT_PDB}")
 
# ---- 3. sanity: report ddG ----
try:
    with open("data/Dif_6eqe_Repair.fxout") as fh:
        lines = [l for l in fh if l.strip()]
    print("ddG line:", lines[-1].split("\t")[1] if len(lines) > 8 else "(check file)")
except Exception as e:
    print("ddG parse skipped:", e)
 
# ---- 4. CRITICAL CHECK: did the two cysteines land at bonding distance? ----
def sg_coord(pdb, resi, chain="A"):
    with open(pdb) as fh:
        for line in fh:
            if (line.startswith("ATOM") and line[12:16].strip() == "SG"
                    and line[21] == chain and int(line[22:26]) == resi):
                return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return None
 
a = sg_coord(MUT_PDB, 79)
b = sg_coord(MUT_PDB, 153)
if a and b:
    d = math.dist(a, b)
    print(f"\nSG(79)-SG(153) distance = {d:.2f} A")
    if d <= 2.5:
        print("  OK -> bonding distance; OpenMM/PDBFixer should form the disulfide.")
    elif d <= 3.5:
        print("  BORDERLINE -> may need to set the disulfide rotamer before MD.")
    else:
        print("  TOO FAR -> FoldX repacked to a non-bonding rotamer. "
              "Fix the rotamer (PyMOL) before building the MD system, or the "
              "'disulfide' will be two free thiols doing nothing.")
else:
    print("\n[warn] could not find both SG atoms - check residue numbering in the mutant PDB")