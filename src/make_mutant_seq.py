wt = ("TNPYARGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRL"
      "ASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNP"
      "SLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQ"
      "ALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEH")

# (PDB_position, wild-type, mutant)  -- your selected 5-mutant
mutations = [
    (233, "N", "K"),
    (77,  "T", "I"),
    (179, "A", "V"),
    (260, "R", "F"),
    (119, "Q", "D"),
]

OFFSET = 29  # string is 1-based PDB pos 29 -> index 0; so index = pos - 29
seq = list(wt)

for pos, wt_aa, mut_aa in mutations:
    idx = pos - OFFSET
    actual = seq[idx]
    if actual != wt_aa:
        raise SystemExit(f"MISMATCH at {pos}: expected {wt_aa}, found {actual} — STOP, do not use this sequence")
    seq[idx] = mut_aa
    print(f"{wt_aa}{pos}{mut_aa}: ok (index {idx})")

mutant = "".join(seq)
print(f"\nlength: {len(mutant)} (should be {len(wt)})")
print(f"WT  ends: ...{wt[-10:]}")
print(f"MUT ends: ...{mutant[-10:]}")
print(f"\nMUTANT SEQUENCE:\n{mutant}")