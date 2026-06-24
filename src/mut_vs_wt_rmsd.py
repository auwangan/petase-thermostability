from Bio.PDB import PDBParser, Superimposer

p = PDBParser(QUIET=True)
wt  = p.get_structure("wt",  "data/alphafold/ispetase_wt_val_924d7_2/ispetase_wt_val_924d7_2_relaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb")[0]
mut = p.get_structure("mut", "data/alphafold/ispetase_5mut_cfbb9_0/ispetase_5mut_cfbb9_0_relaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb")[0]

def ca(chain): return {r.id[1]: r["CA"] for r in chain if "CA" in r}
w, m = ca(next(wt.get_chains())), ca(next(mut.get_chains()))
keys = sorted(set(w) & set(m))
sup = Superimposer()
sup.set_atoms([w[k] for k in keys], [m[k] for k in keys])
print(f"matched: {len(keys)} | mutant-vs-WT Cα RMSD: {sup.rms:.3f} Å")