# Selecting the plate: diversity over score

## The argument

Once designs pass a deterministic epitope gate, every one of them satisfies the
geometric criterion. If no available score discriminates further — which is what
`scoring-is-not-location.md` establishes — then ranking the plate by score is
selection on a variable with no demonstrated relevance to the outcome.

Under "no score predicts success", the objective that *is* defensible is the
number of **independent shots on goal**. If designs sharing a binding mode
succeed or fail together, then what determines whether the plate returns
anything is how many distinct modes it contains, not how they are ordered.

This is a decision-theoretic argument, not an empirical one, and should be
presented as such. It becomes wrong the moment someone demonstrates a calibrated
score.

## Method

1. **Assemble the gate survivors** with their sequences, CDR1/CDR2/CDR3, lineage
   label, and epitope footprint.

2. **Composite distance.** Weight the features that plausibly determine whether
   two designs fail for the same reason:

   ```
   0.40 · CDR3 Levenshtein (normalised)
   0.15 · CDR3 length difference
   0.20 · CDR1+CDR2 Levenshtein
   0.25 · epitope footprint Jaccard distance
   ```

   CDR3 dominates because it usually dominates binding; footprint is included
   because two designs contacting different epitope residues are different
   experiments even with similar sequences. Adjust and state the weights.

3. **Quality floor before diversity.** Diversity without a floor becomes "the
   worst designs are the most unusual". Apply a floor within the survivors —
   `atom_purity ≥` the survivor median worked well. Report what it excludes.

   This floor matters more than it looks: it is what stops the plate filling
   with large, sticky, low-purity interfaces, the exact failure mode the
   ipSAE↔purity anticorrelation warns about.

4. **Cluster into N clusters** (average linkage or k-medoids on the composite
   distance) and take one representative per cluster. Break ties on **geometry**
   — recall, then atom_purity, then contact count — explicitly *not* on the score
   you have declined to select on.

5. **Cap lineage concentration** at ~30% of the plate, backfilling from the
   next-most-distant unused clusters. Even if the cap never binds, checking it
   is what tells you whether it needed to.

## What the comparison looked like

Diversity-first versus score-ranked, same 548 gate survivors, 96 seats:

| | diversity | score-ranked |
|---|---|---|
| effective lineages | **64.0** | 12.3 |
| largest lineage share | **6.2%** | 27.1% |
| distinct epitope footprints | 15 | 11 |
| min pairwise distance | **0.2535** | 0.1258 |
| mean nearest-neighbour distance | **0.3396** | 0.2368 |
| purity | **0.840** | 0.724 (p < 1e-6) |
| atom_purity | **0.962** | 0.913 |
| recall | 0.853 | 0.843 (n.s.) |
| median binding_confidence | 0.0022 | 0.6595 |

Two things to notice.

**The geometry did not get worse.** The diversity plate was equal or better on
every occupancy metric. The score-ranked plate led only on raw contact count —
at lower purity, which is the sticky-interface artifact, not an advantage. Only
35 of its 96 designs passed the atom_purity floor at all.

**The score-ranked plate was paying twice for near-duplicates.** Its minimum
pairwise distance was half the diversity plate's; the nearest-neighbour gap was
43% smaller.

The median score dropped 306×. That is the intended trade, and it costs nothing
*given* the null results — a caveat that belongs in the same sentence.

## Expected independent modes

Score this on a **neutral partition**: cluster all gate survivors once and count
how many of those clusters each plate occupies. Grading the diversity plate on
its own clustering would be circular.

With modes occupied 71 vs 19, and a per-mode success probability p:

| p | diversity | score-ranked | E[modes hit] |
|---|---|---|---|
| 0.10 | 0.99944 | 0.86492 | 7.1 vs 1.9 |
| 0.15 | 0.999990 | 0.954401 | 10.65 vs 2.85 |

**Report E[modes hit], not P(≥1).** Both plates probably return *something*;
the difference is how many independent starting points survive to lead
optimisation. P(≥1) makes the gap look more dramatic than the decision warrants.

The within-cluster-correlated / between-cluster-independent assumption is an
assumption. Label it. It is still the right decision procedure in the absence of
a calibrated score, and the conclusion is not sensitive to the exact p.

## Final QC — run it on the batch, after all edits

Diversity selection surfaces designs that score ranking had buried, so the
liability rate goes **up**: 19.8% framework anomalies versus 4.2% on the
score-ranked plate. This is the one genuine cost, and it is fixable.

For each design in the final plate check:
- framework canonical (FR4 motif intact, conserved positions present)
- exactly the canonical number of cysteines
- no N-glycosylation sequons
- CDR3 unchanged by any repair
- sequence length unchanged

Then **make the repairs and re-run the whole QC**, because a repair can create a
new defect — see `failure-modes.md`. Before accepting a repair, verify the edited
positions are outside the paratope by recomputing binder-side contacts from the
cached structure. Framework reversions normally are; showing it costs one loop.

Output both a table (with original sequence, repaired sequence, and the explicit
edit list) and a FASTA. The edit list is what lets someone else audit the repair
without re-deriving it.
