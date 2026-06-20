# Week 2 Reflection

*(Draft — edit in my own voice. Written at the close of Week 2, after Phase 1 came together.)*

## What I set out to do
Diagnose why Week 1's pipeline rejected all 30 of its top candidates, rebuild it, and pick a direction. Reading and a decision document were stretch goals, not expectations.

## What actually happened
A lot more than that. The Week 1 failure turned out to be a single inverted sign in my ThermoMPNN filter — I'd been selecting the *most destabilizing* mutations as if they were the best. Fixing it didn't just recover the week; it became a full rebuild: a two-method consensus pipeline (ThermoMPNN + FoldX), noise-hardening with 3-run averaging, function checks against the substrate-binding cleft, and multi-mutant design with epistasis testing. I finished the week with a validated 5-mutation lead (N233K / T77I / A179V / R260F / Q119D, −6.57 kcal/mol), structural validation of the load-bearing mutations, and a written Phase 1 summary.

## What I learned — technical
- **Conventions kill, silently.** One negative-vs-positive sign flip invalidated an entire candidate set without throwing a single error. Anchoring the pipeline to known experimental stabilizers as positive controls is what made the bug visible. I won't trust a predictor's output again without checking its sign convention against a case I already know the answer to.
- **Two methods beat one — because they fail differently.** Consensus rejected the physically implausible mutations ThermoMPNN over-scored, and the apparent ML-vs-physics disagreement was the exact clue that led me to the bug. Agreement between independent methods is worth more than confidence from one.
- **Test combinations, don't sum them.** Epistasis is real: the 6-mutant lost 2.4 kcal/mol to antagonism I would never have seen by adding up single mutations.
- **The biggest number isn't the best design.** The 7-mutant scored −8.28 but leaned on things FoldX can't see (a free cysteine, a groove-intruding residue). The 5-mutant is "worse" on paper and better in reality, because I can defend every mutation in it.
- **Stability and activity are different questions.** A mutation can stabilize the fold and still wreck the active site. The orientation check — does the side chain point into the groove or away from it? — is the entire reason Q119D stayed and Q127L was retired, despite both being near the cleft.

## What worked — process
- Asking "what does this do, and why?" before running each script. I actually understand my own pipeline now, which is the only reason I could write the Phase 1 summary from memory.
- Treating the failure as a diagnostic problem with named hypotheses, not just "it's broken."
- Committing and logging as I went, so nothing was lost across breaks.

## What was hard / what I'd change
- I let the day-numbering and the calendar drift apart and lost track of where I was more than once. A one-line "current state" note at the start of each session fixes this — which is exactly what the handoff doc is now for.
- A few times I ran a script before fully understanding its output, then had to backtrack. Reading the intent first is slower per step but faster overall.

## On pacing
I took breaks this week, some of them because the motivation wasn't there. I'm writing that down honestly instead of pretending it was a clean sprint — it wasn't, and the work still got done. The commit history shows it: the project sat untouched for days and was exactly where I left it when I came back. That's the thing worth remembering for next time the motivation dips — the work waits, it doesn't expire, and picking the thread back up after a gap is a skill, not a failure.

## Where I am now
Phase 1 is done and banked on GitHub. The lead is a validated 5-mutant. Next is outreach (first email drafted) and Phase 2 planning — MD, substrate docking, AlphaFold.