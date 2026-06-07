#Load both CSVs: data/ThermoMPNN_inference_6eqe_Repair.csv and data/residue_features.csv
#Fix the numbering offset — add 29 to each ThermoMPNN position so it matches PDB numbering
#Join the two tables on PDB position, so each mutation row also carries that residue's RSA, location, distance-to-triad, and near_active_site flag
#Apply filters — keep only mutations that are:
#stabilizing: ThermoMPNN ddG_pred > 0.5
#not near the active site (near_active_site flag = "NO")
#not on the functional blacklist: exclude positions 185 and 257
#not catalytic triad: exclude 160, 206, 237 (the active-site filter likely catches these already, but be explicit)
#Rank survivors by ddG_pred, most stabilizing first
#Save to data/phase1_candidates.csv and print the top ~20
import pandas as pd
thermo_df = pd.read_csv("data/ThermoMPNN_inference_6eqe_Repair.csv")
matched_thermo_df = thermo_df.copy()
matched_thermo_df["position"] = matched_thermo_df["position"] + 29
residue_df = pd.read_csv("data/residue_features.csv")
merged_df = pd.merge(matched_thermo_df, residue_df, on="position")
candidates = merged_df[
    (merged_df["ddG_pred"] > 0.5) &
    (merged_df["near_active_site"] == "NO") &
    (~merged_df["position"].isin([160, 206, 237, 185, 257]))
]
print(candidates.shape)
candidates_sorted = candidates.sort_values("ddG_pred", ascending=False)
print(candidates_sorted[["position","wildtype","mutation","ddG_pred","location","rsa"]].head(20))
top_candidates = candidates_sorted.head(30)
top_candidates.to_csv("data/phase1_candidates.csv", index=False)
