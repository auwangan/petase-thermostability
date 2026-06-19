import numpy as np
from Bio.PDB import PDBParser

chain = PDBParser(QUIET=True).get_structure("x", "data/6eqe_Repair.pdb")[0]["A"]

# groove center = centroid of the binding-cleft residues
CLEFT = [87, 159, 160, 161, 185, 206, 208, 237, 238, 241, 280]
pocket = np.mean([chain[r]["CA"].coord for r in CLEFT], axis=0)

def orientation(resnum, label):
    res = chain[resnum]
    ca, cb = res["CA"].coord, res["CB"].coord
    v_side   = cb - ca                  # which way the side chain projects
    v_pocket = pocket - ca              # toward the groove center
    cos = np.dot(v_side, v_pocket) / (np.linalg.norm(v_side) * np.linalg.norm(v_pocket))
    angle = np.degrees(np.arccos(np.clip(cos, -1, 1)))
    d_ca = np.linalg.norm(ca - pocket)
    d_cb = np.linalg.norm(cb - pocket)
    facing = "TOWARD groove (activity risk)" if angle < 90 else "AWAY from groove (ok)"
    print(f"{label}: angle {angle:3.0f}°  | CA->groove {d_ca:.1f} A, CB->groove {d_cb:.1f} A  -> {facing}")

orientation(119, "Q119D")
orientation(127, "Q127L")