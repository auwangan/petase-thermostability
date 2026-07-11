"""OpenMM Modeller re-indexes residues from 1. Restore 29-293 PDB numbering
on the minimized structures, but ONLY for protein residues (not water/ions)."""
import sys
IN, OUT = sys.argv[1], sys.argv[2]
OFFSET = 28
SOLVENT = {"HOH","WAT","NA","CL","K","MG"}
n = 0
with open(IN) as fin, open(OUT, "w") as fout:
    for line in fin:
        if line.startswith(("ATOM","HETATM","TER","ANISOU")) and len(line) > 26:
            resn = line[17:20].strip()
            if resn not in SOLVENT:
                line = line[:22] + f"{int(line[22:26]) + OFFSET:>4}" + line[26:]
                n += 1
        fout.write(line)
print(f"renumbered {n} protein atom records -> {OUT}")
cys = sorted({int(l[22:26]) for l in open(OUT)
              if l.startswith("ATOM") and l[17:20].strip()=="CYS" and l[12:16].strip()=="CA"})
print(f"cysteines now at: {cys}")
