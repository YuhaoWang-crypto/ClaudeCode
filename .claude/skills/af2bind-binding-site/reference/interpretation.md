# Reading AF2BIND output honestly

## p(bind) is a ranking, not a probability

The head is a logistic regression fitted on one training distribution of
structures. Its output looks like a probability and is not one in any calibrated
sense. The absolute value moves with protein length, structure quality, whether
side chains were masked, and which AlphaFold parameter release produced the
features. A p(bind) of 0.6 on one target and 0.6 on another are not comparable
claims.

What survives all of that is the **ordering within one run**. Use the rank.

Practical consequences:

- Do not set a fixed p(bind) cut-off and carry it between targets.
- Do not report "residue X has a 72% chance of being a binding site".
- Do compare a residue against the rest of its own profile.

## What to actually report

**The top ~15 residues.** That is the unit the upstream notebook prints and the
one that behaves. Below roughly the top 25 the list starts collecting
second-shell and surface residues.

**The pockets, not the residues.** `report.json` clusters the top-scoring
residues into spatially distinct sites with single-linkage on side-chain
centroids, ranks them by summed score, and gives each a PyMOL selection. A
ranked residue list hides the fact that positions 3, 7 and 12 may be a second
site on the other face of the protein. Note this clustering is added by this
pipeline; it is not part of the published method, and `link_cutoff` (default
10 Å) is a tunable, not a physical constant.

**The shape of the profile.** A real site shows a handful of residues clearly
separated from the bulk. If the maximum sits close to the body of the
distribution, the honest answer is "no confident site found", not the top-ranked
residue. `report.json` reports `max`, `mean`, and counts above 0.5 and 0.9
precisely so this is checkable.

## Measured behaviour

Numbers below were produced by this pipeline, `mask_sidechains=True`, AlphaFold
parameters 2021-07-14, seed 0 unless stated. Ground truth is every residue with a
heavy atom within 5 Å of the crystallographic ligand, waters and crystallisation
additives excluded, ligands assigned to their nearest protein chain.

| | 6W70 chain A | 1M17 chain A |
|---|---|---|
| what it is | ABLE, a designed apixaban binder (the paper's own example) | EGFR kinase domain with erlotinib |
| residues scored | 126 | 312 |
| ligand | GG2 (apixaban) | AQ4 (erlotinib) |
| contact residues (positives) | 25 of 126 (19.8%) | 20 of 312 (6.4%) |
| ROC-AUC | 0.978 | 0.934 |
| average precision | 0.917 | 0.547 |
| precision@10 | 1.00 | 0.70 |
| precision@15 | 0.93 | 0.67 |
| enrichment over random @10 | 5.0× | 10.9× |
| max p(bind) | 0.980 | 0.916 |
| mean p(bind) | 0.198 | 0.087 |
| pockets found | 1 | 1 |

Three things in that table are worth internalising.

**The absolute scale moves with the target.** Max p(bind) is 0.98 on the small
designed binder and 0.92 on the kinase; mean p(bind) more than halves. Same
method, same settings, different scale. This is why a fixed threshold does not
transfer.

**Average precision falls while ROC-AUC barely does.** Going from 20% positives
to 6% positives costs 0.37 AP and only 0.04 AUC. AP tracks the base rate, so it
is the honest number on a large protein and the one to quote. A high AUC on a
sparse target can coexist with a top-15 list that is a third wrong.

**The top of the list is chemistry, not coincidence.** On EGFR the top-ranked
residues are L694, V702, K721, T766, M769, L820, T830, D831 and A719, which is
the canonical ATP pocket that erlotinib occupies, recovered without the ligand
ever being shown to the model. The bait-activation panel for those residues is
dominated by the aromatic baits F, W and Y, consistent with the quinazoline
ligand that actually binds there.

Figures: `figures/af2bind/af2bind_6w70_A.png`, `figures/af2bind/af2bind_1m17_A.png`.

Ensembling all ten folds instead of seed 0 changed these targets very little:
ROC-AUC moved by under 0.004 on both, and top-15 precision on EGFR fell slightly
(0.67 to 0.53) while top-10 was unchanged. Do not expect the ensemble to rescue
a weak result.

Read these as a smoke test of the implementation, not as a benchmark. Two targets
is not an evaluation, and both are holo structures, which is the easy case. The
paper's own benchmark numbers are the ones to cite.

## Ensembling across the 10 folds

The released archive contains ten independently trained folds. `--seeds
0,1,...,9` averages their logits before the sigmoid.

Ensembling reduces fold-to-fold variance in the ranking. It does **not** make
p(bind) calibrated, and it is not the configuration the paper reports, so use
seed 0 when the goal is to reproduce published behaviour and the ensemble when
the goal is a more stable ranking for a decision you are about to act on. Report
which you used.

## Failure modes, in the order they bite

**Disordered regions on predicted models.** The strongest known false-positive
source. Disordered AlphaFold tails are exposed and unpacked, which is the local
environment the score rewards. Always `--min-plddt 70` on an AlphaFold model, and
split multi-domain proteins by domain.

**Interfaces that are not small-molecule sites.** A buried protein-protein
interface, a crystal-contact surface, or a nucleic-acid binding groove can score
well. The method was trained on small-molecule pockets; it does not know what it
is looking at. Check whether a high-scoring patch is a genuine cavity.

**Very small or very large targets.** A short peptide has no pocket to find. A
large multi-domain protein dilutes the signal and strains the quadratic memory.
Domain-wise scoring is both cheaper and better.

**Apo relaxation.** Scores on an apo structure typically run lower than on the
matched holo structure. That is expected, not a bug, and it is exactly the
setting the method was built for. Do not compare apo and holo p(bind) values as
if they were on one scale.

## When the answer is "I don't know"

Say so. Cases that warrant it:

- The profile is flat: no residue stands clearly above the bulk.
- The top residues do not cluster: scattered singletons rather than a site.
- The input is an untrimmed AlphaFold model with long disordered stretches.
- A holo control on a related protein failed to recover its known site.

A ranked list always exists. That is not the same as a site existing.

## Downstream use

The natural next step is docking: take the top pocket's centroid from
`report.json` as the box centre and its residue extent as the box size. The
bait-activation profile in `bait_activations.csv` hints at pocket chemistry and
can help choose a library to screen, but it is a reported correlation, not a
calibrated ligand-property predictor. Treat it as a prior, then let docking or
experiment decide.
