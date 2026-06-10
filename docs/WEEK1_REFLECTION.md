Week 1 Reflection — Phase 1 First Pass
Period: Days 8-14 (June 4-10, 2026)
Status: Phase 1 first pass complete. Pipeline built, validated, produced a real result.

What I built this week
A full predict -> filter -> validate pipeline for finding thermostabilizing mutations in IsPETase:

ThermoMPNN mutation generation — site-saturation scan of all 265 residues x 20 amino acids = 5300 predicted mutations, each with a predicted ddG. (Already validated on the known S238F mutation: both tools called it stabilizing.)
Structural filter (compute_residue_features.py) — classifies each residue by RSA (core vs surface, cutoff 0.20), distance to the catalytic triad (Ser160/Asp206/His237), and Kyte-Doolittle hydrophobicity. This caught the W257 blind spot: ThermoMPNN's top prediction was a substrate-binding tryptophan that my triad-distance filter missed (15.6 A from the triad), so I added an explicit functional-residue blacklist.
Candidate filter (filter_candidates.py) — merges ThermoMPNN predictions with the structural features (applying the +29 PDB/ThermoMPNN numbering offset), excludes active-site and blacklist positions, ranks by ddG, outputs top 30.
FoldX validation (run_foldx_batch.py) — physics-based double-check of the top 30 via BuildModel.


The main result
FoldX rejected all 30 candidates — in both v1 and v2. ThermoMPNN said every candidate was stabilizing; FoldX said every one was destabilizing. The two tools disagreed on the sign for all 30. No consensus stabilizer in ThermoMPNN's top ranks.
The pattern explained why: ThermoMPNN over-predicts stabilization for chemically implausible core mutations — mutations to proline, charged residues into the buried core, and polar/charged residues into the core generally. These are exactly the mutation types FAST-PETase (the gold-standard experimental result) avoided; all five of its mutations target flexible surface regions, none are core prolines or buried charges.
So I built v2 of the filter with chemical-sanity rules: drop mutations to proline, and drop hydrophilic residues (D/E/K/R/N/Q/H/S/T/Y) going into the core (the "don't bury a charge / oil-likes-oil" principle). v2 worked — the best FoldX score improved from +2.5 (v1) to +1.2 (v2, the position-52 surface mutations, near-neutral). But still nothing crossed into stabilizing.
Conclusion: ThermoMPNN's highest-ranked predictions for IsPETase do not survive physics-based validation. A single-tool ML pipeline would have advanced 30 false positives. The two-tool consensus design caught this.

What I got wrong / had to fix

Ranking by ThermoMPNN ddG magnitude was the wrong selection criterion. It systematically surfaced buried-position mutations (78, 81, 109, 149, 52...) that FoldX rejects. The positions themselves are the problem, not just the substitutions — every time I blocked one bad substitution type, ThermoMPNN proposed another at the same buried positions (proline -> charge -> polar -> glycine).
The triad-distance filter alone missed W257, a substrate-binding residue. Needed an explicit functional blacklist.
Chemical-sanity rules have diminishing returns — I could chase glycine-into-core, then alanine, forever. The honest move was to let FoldX adjudicate rather than encode every edge case.

What I actually learned (not just did)

The two-tool consensus principle, concretely. I don't automatically trust FoldX over ThermoMPNN — neither is ground truth, both are predictions with different failure modes. When two independent methods disagree completely, the prediction is unreliable, full stop. ThermoMPNN's picks lost here because three independent lines (FoldX physics, basic chemistry, and FAST-PETase experiment) all pointed the same way. If chemistry+literature+ThermoMPNN had agreed and FoldX disagreed, I'd distrust FoldX instead.
Core vs surface chemistry: the core is dry/oily and wants hydrophobic residues; the surface is wet and tolerates charges. Burying a charge or polar residue destabilizes. This is just oil-and-water applied inside a protein, and it explained every one of FoldX's rejections.
ThermoMPNN's scale is compressed — a real stabilizer (S238F) scored only +0.054, so the absolute ddG number is not trustworthy; ranking matters more than magnitude, and even the ranking is biased toward buried positions for this protein.

Coding honesty check
I wrote three real analysis scripts this week (compute_residue_features.py, filter_candidates.py/v2, run_foldx_batch.py) — with hints, not copied code. Bugs I debugged myself: initialize-list-outside-loop (multiple times), .extend vs .append, raw SASA vs normalized RSA, atom-atom distance, pandas boolean-filter precedence with ~, the FoldX append-to-file quirk, the iloc row-skip mismatch. Python feels meaningfully more fluent than Day 7 — I can now write a multi-step pandas + Biopython script and reason about why a bug happens, not just that it happened. Still lean on hints for new patterns (subprocess, complex boolean indexing), which is fine.

Phase 1 scope — finalized decision
The "no consensus winner" result changes the plan. For the next pass:

Stop ranking by ThermoMPNN ddG alone. Rank by consensus (require both tools to agree on sign) or run FoldX earlier on a broader, chemically-clean pool and rank by FoldX.
The good candidates are probably the modest-ThermoMPNN, surface mutations I'm currently filtering past — the FAST-PETase pattern. My threshold selected for the wrong signal.
Phase 2 direction: use MD-identified flexible regions as the mutation-site selector (FAST-PETase's mutations all hit thermolabile regions), rather than blind ddG ranking. This is a more principled site-selection strategy than ThermoMPNN magnitude.

Open items going into Week 2

 Friend skill-check / role assignment — deferred (he is overloaded with summer-school work; no pressure, revisit when he is free).
 Verify the flagged unverified claims in literature/lu-2022-fast-petase.md (Tm number, per-residue burial, 159->5 narrowing, epistasis) against the real Lu 2022 paper / free Supplementary Methods.
 Resend Leveson-Gower with an honest hook (the stability-predictor tool was not his).
 Mentor: Furst agreed to 10 min of feedback (reply 7 Jun) — next contact when the ~2-page methodology writeup is ready (late July).
 Decide the next-pass ranking strategy (consensus vs FoldX-first) and re-run.

Honest overall
Week 1 overshot the technical plan — I completed Phase 1 first-pass work that was scheduled for weeks later. The headline isn't a stabilizing mutation (I have none yet); it's a methodological result: ThermoMPNN's top-ranked IsPETase predictions are systematically rejected by FoldX, consistent with known false-positive modes and with the FAST-PETase experimental precedent. That is a genuine, defensible finding and a good thing to show Furst.