import pandas as pd
mpnn = pd.read_csv("data/ThermoMPNN_inference_6eqe_Repair.csv")

# negatives in PDB numbering: (pdb_position, wildtype, mutation)
negatives = [
    (257, "W", "A"),
    (106, "F", "A"),
    (101, "L", "A"),
    (68,  "V", "D"),
    (83,  "I", "R"),
    (145, "I", "D"),
]

for pdb_pos, wt, mut in negatives:
    row = mpnn[(mpnn["position"] == pdb_pos - 29) & (mpnn["wildtype"] == wt) & (mpnn["mutation"] == mut)]
    print(f"{wt}{pdb_pos}{mut}: ThermoMPNN ddG = {row['ddG_pred'].values}")