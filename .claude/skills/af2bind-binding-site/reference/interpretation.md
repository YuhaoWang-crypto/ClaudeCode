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

## The calibration panel

`python -m af2bind_pipeline.benchmark` runs ten drug targets spanning target
classes, each with its ligand named explicitly. Reproduce with
`AF2BIND_GPU=L40S`; raw output is in `figures/af2bind/panel.json`.

| target | class | L | positives | ROC-AUC | AP | P@10 | enrich@10 | max p | mean p |
|---|---|---|---|---|---|---|---|---|---|
| 6W70 (GG2, apixaban) | designed binder | 126 | 25 (19.8%) | 0.977 | 0.916 | 1.00 | 5.0× | 0.980 | 0.197 |
| 3LN1 (CEL, celecoxib) | large glycoprotein | 552 | 24 (4.3%) | 0.959 | 0.544 | 0.70 | 16.1× | 0.963 | 0.119 |
| 3ERT (OHT, tamoxifen) | nuclear receptor | 247 | 22 (8.9%) | 0.954 | 0.728 | 0.90 | 10.1× | 0.956 | 0.108 |
| 2RH1 (CAU, carazolol) | GPCR | 442 | 20 (4.5%) | 0.938 | 0.628 | 0.80 | 17.7× | 0.961 | 0.217 |
| 1IEP (STI, imatinib) | kinase | 274 | 28 (10.2%) | 0.935 | 0.611 | 0.70 | 6.9× | 0.987 | 0.145 |
| 1M17 (AQ4, erlotinib) | kinase | 312 | 20 (6.4%) | 0.934 | 0.542 | 0.70 | 10.9× | 0.917 | 0.088 |
| 4EIY (ZMA) | GPCR | 390 | 16 (4.1%) | 0.918 | 0.519 | 0.50 | 12.2× | 0.929 | 0.160 |
| 1HXW (RIT, ritonavir) | obligate-dimer site | 99 | 15 (15.2%) | 0.866 | 0.504 | 0.40 | 2.6× | 0.970 | 0.260 |
| 1STP (BTN, biotin) | small soluble | 121 | 17 (14.1%) | 0.858 | 0.643 | 0.80 | 5.7× | 0.719 | 0.139 |
| 3HS4 (AZM, acetazolamide) | metalloenzyme | 257 | 32 (12.4%) | 0.779 | 0.424 | 0.70 | 5.6× | 0.968 | 0.090 |

Every target clears random by a wide margin, and enrichment at the top of the
list runs 2.6× to 17.7×. Four things in the spread are worth knowing before you
read a number on a new target.

**Deep, enclosed pockets are the best case.** The nuclear-receptor LBD and the
GPCR orthosteric sites — buried cavities lined on all sides — score at or near
the top, and the paper's own designed binder is a purpose-built version of the
same thing. That is where to trust a high score most.

**Shallow and solvent-exposed sites are the worst case.** Carbonic anhydrase is
the weakest target here (AUC 0.78) despite a textbook druggable active site: the
site is a wide cone whose specificity comes from zinc coordination, and AF2BIND
scores amino-acid company, not metal chemistry. Biotin in streptavidin is small
and largely backbone- and water-coordinated, and it is the only target whose max
p(bind) does not reach 0.9. Expect the method to underrate sites whose binding
is driven by a metal, a cofactor or ordered water rather than by side-chain
packing.

**Scoring one chain of an obligate oligomer surfaces the wrong pocket.** HIV-1
protease has the lowest precision@10 (0.40) and the lowest enrichment (2.6×), and
the reason is instructive rather than a metric artefact. The active site is
found: D25, I47, G49, I50 and I84 are all high, which is the textbook
pharmacophore. But I3, L5 and L97 score just as high, and those form the
interdigitated β-sheet that dimerises the enzyme — a hydrophobic surface that is
buried in the biological unit and exposed in the single chain AF2BIND was given.
On any obligate oligomer, score the assembly's chains knowing that the
oligomerisation interface will compete with the real site, and cross-check a top
pocket against the biological assembly before acting on it.

**A high AUC on a big protein still means a messy top-15.** 3LN1 has the second
best AUC in the panel (0.959) and an average precision of 0.54, because 24
positives among 552 residues is a 4.3% base rate. AUC is nearly base-rate
independent and AP is not, which is exactly why both are reported. On a large
target, read AP and precision@N; the AUC will flatter you.

Two caveats on the panel itself. Every target is a holo structure with the
ligand stripped, which is the easy setting — apo and predicted inputs are harder.
And ten targets chosen by hand is a calibration aid, not an evaluation; the
paper's benchmark is what to cite.

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
