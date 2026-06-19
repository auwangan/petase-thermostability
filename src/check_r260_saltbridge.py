from Bio.PDB import PDBParser

chain = PDBParser(QUIET=True).get_structure("x", "data/6eqe_Repair.pdb")[0]["A"]

arg = chain[260]                                  # R260
pos_atoms = [arg[a] for a in ("NE", "NH1", "NH2") if a in arg]   # guanidinium N+

neg = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}             # carboxylate O-
hits = []
for res in chain:
    rn = res.get_resname()
    if rn in neg:
        for oname in neg[rn]:
            if oname in res:
                for n in pos_atoms:
                    d = n - res[oname]
                    hits.append((round(d, 2), f"{rn}{res.get_id()[1]}", oname, n.get_name()))

hits.sort()
print(f"{'dist(A)':>8}  partner       contact")
for d, partner, oname, nname in hits[:8]:
    flag = "  <-- SALT BRIDGE" if d <= 4.0 else ("  (long-range)" if d <= 6.0 else "")
    print(f"{d:>8}  {partner:<11} R260:{nname}–{partner}:{oname}{flag}")

closest = hits[0][0] if hits else None
print(f"\nclosest R260(+) to any Asp/Glu(-): {closest} A")
if closest is None or closest > 6.0:
    print("=> no carboxylate nearby. R260F looks structurally clean — lead holds.")
elif closest <= 4.0:
    print("=> R260 anchors a SALT BRIDGE. R260F loses it — verify/demote before trusting the design.")
else:
    print("=> borderline long-range contact. Mild concern, worth an eyeball.")