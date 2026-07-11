"""Write a protein-only topology matching the protein-only DCD atom subset.
Must use the SAME selection as md_v2_05_production.py's atomSubset."""
import sys
IN, OUT = sys.argv[1], sys.argv[2]
SOLVENT = {"HOH","WAT","NA","CL","K","MG"}
n = 0
with open(IN) as fin, open(OUT, "w") as fout:
    for line in fin:
        if line.startswith(("ATOM","HETATM")):
            if line[17:20].strip() in SOLVENT:
                continue
            n += 1
            fout.write(line)
        elif line.startswith(("TER","END")):
            fout.write(line)
print(f"wrote {n} protein atoms -> {OUT}")
