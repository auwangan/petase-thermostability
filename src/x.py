import pandas as pd
mpnn = pd.read_csv("data/ThermoMPNN_inference_6eqe_Repair.csv")

# controls in PDB numbering: (pdb_position, wildtype, mutation)
controls = [
    (121, "S", "E"),
    (186, "D", "H"),
    (224, "R", "Q"),
    (280, "R", "A"),
    (159, "W", "H"),
    (238, "S", "F"),
]

# offset sanity check using N233K (known ThermoMPNN ddG = -1.198)
chk = mpnn[(mpnn["position"] == 233 - 29) & (mpnn["wildtype"] == "N") & (mpnn["mutation"] == "K")]
print("N233K offset check:", chk["ddG_pred"].values)

# pull the six controls
for pdb_pos, wt, mut in controls:
    row = mpnn[(mpnn["position"] == pdb_pos - 29) & (mpnn["wildtype"] == wt) & (mpnn["mutation"] == mut)]
    print(f"{wt}{pdb_pos}{mut}: ThermoMPNN ddG = {row['ddG_pred'].values}")