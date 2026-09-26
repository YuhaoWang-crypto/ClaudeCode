# Does a structural ensemble buy anything over 2D chemistry?

A CDK2 benchmark that is built to be able to answer **no**.

The pitch for structure-based features in virtual screening is that a docked pose
tells you *where* a ligand touches the pocket, and that this is information a
molecular fingerprint cannot contain. The pitch is rarely tested against the
baseline that matters: a plain ECFP4 model on the same compounds, under a split
that does not leak analogue series between train and test.

This repository runs that test. Every design decision below exists to keep the
negative result reachable.

---

## The question, stated so it can fail

> Given the same compounds and the same labels, does a per-residue interaction
> fingerprint derived from ensemble docking predict CDK2 potency better than
> ECFP4 + physicochemical descriptors — under a scaffold split, against a
> permutation null, with error bars?

Three separate things are being asked, and they have different answers:

| # | Question | Where it is answered |
|---|---|---|
| 1 | Do the conformational slices actually differ in a way that matters? | `s03` cross-docking control |
| 2 | Does docking alone rank compounds? | `s07` zero-parameter baselines |
| 3 | Do 3D features add anything **over** 2D? | `s07` arm comparison |

Question 1 is a prerequisite: if any receptor reproduced any ligand's pose, the
ensemble would be redundant and question 3 would be uninteresting.

---

## Pipeline

```
s01  build the dataset      ChEMBL CDK2 (CHEMBL301) -> 2,200 compounds, 935 scaffolds
s02  build the ensemble     6 PDB slices, superposed on the ligand contact shell
s03  cross-docking control  every native ligand into every slice  <- the premise test
s04  pick + embed ligands   whole scaffold groups, label-blind, ETKDGv3 + MMFF94s
s05  ensemble docking       ligand x slice, checkpointed, smina
s06  interaction features   per-residue contact channels from the poses
s07  model + honesty checks scaffold split, permutation null, repeated CV
s08  Boltz-2 arms           affinity head vs pIC50; predicted conformers as slices
```

Reproduce with `bash setup.sh` then the stages in order. Every stage is
deterministic given the seeds recorded in its header.

---

## What the data actually is

ChEMBL CDK2, human, single-protein target, `pchembl_value` present, activity
type in {IC50, Ki, Kd}, no censored (`>`/`<`) records.

* **2,200 compounds** after aggregating replicate measurements by median.
* Replicates are real: 396 compounds carry more than one measurement, and their
  median absolute deviation is **0.30 log units**. That is the noise floor. Any
  model claiming RMSE far below ~0.3 on this data is fitting the assay, not the
  chemistry, and the benchmark reports it next to the model numbers rather than
  quietly beating it.
* 935 Bemis-Murcko scaffolds; 36 of them have 10+ members, 662 are singletons.

Aggregating by median rather than mean is deliberate — ChEMBL potency
distributions have outliers from single bad assays, and the mean chases them.

---

## The ensemble, and why these six structures

| slice | ligand | pocket RMSD to ref | global RMSD | what it contributes |
|---|---|---|---|---|
| 1HCK | ATP | 0.00 (reference) | 0.00 | the cofactor-bound ground state |
| 2VTA | LZ1 | 0.44 | 0.43 | fragment-bound, tight pocket |
| 1KE5 | LS1 | 0.54 | 0.93 | oxindole inhibitor |
| 1AQ1 | STU | 0.86 | 0.99 | staurosporine, pocket expanded |
| 3PXF | 2AN | 1.48 | 2.49 | **ATP site empty** — ANS bound allosterically |
| 1FIN | ATP | 2.01 | 4.73 | cyclin-A-bound active conformation |

Three decisions in that table are worth stating, because each was a measurement
rather than a default:

**The alignment frame is the ligand contact shell, not all CA atoms.** Aligning
globally put 1FIN 3.9 Å from the reference and walked its ATP site out of the
shared box. The culprit is real biology — cyclin A remodels the activation loop
by 8–15 Å at residues 145–165 — but those residues have no business defining the
frame for a docking box. Restricting the fit to the 30 residues within 6 Å of
the reference ATP drops 1FIN's *pocket* RMSD to 2.01 Å while leaving its global
RMSD at 4.73 Å. Both numbers are reported; the second is the honest statement of
how different these structures are.

**3PXF contributes a receptor but no box geometry.** Its detected ligand sits
11.5 Å from the hinge: both ANS copies are in the allosteric pocket under the
C-helix and the ATP site is empty. It stays in the ensemble deliberately, as the
negative control for "does an inappropriate conformer degrade the ensemble?",
but it is excluded from sizing the shared box.

**One shared box for every slice.** Per-slice boxes would make docking scores
incomparable across the ensemble, which is precisely the confound the benchmark
is testing for.

---

## Control: do the slices carry distinct information?

Every ATP-site co-crystal ligand docked into every slice, symmetry-corrected
RMSD to the crystallographic pose (`results/crossdock_control.json`):

```
             1HCK     1AQ1     2VTA     1FIN     1KE5     3PXF
   ATP      0.60*    2.38     2.87     3.57     3.72     4.79
   STU      9.38     0.27*    5.02     5.53     9.54    11.81
   LZ1       n/a     9.93     2.93*     n/a     3.14     1.25
   LS1       9.41     2.86     9.52     0.81     0.58*    6.64
                                            (* = self-dock)
```

**Self-docking: 3 of 4 pass under 2 Å, median 0.59 Å.** The setup can reproduce
crystallography. The failure is LZ1 at 2.93 Å — a 9-heavy-atom fragment with too
little shape complementarity for docking to pin down. That is a real limitation
of fragment docking, reported rather than hidden.

**Cross-docking: median 4.9 Å, only 11% under 2 Å.** A ligand docked into the
wrong conformer usually lands in the wrong place. The slices are *not*
interchangeable, so the ensemble premise holds and question 3 is worth asking.

The score matrix says the same thing in energy units — staurosporine scores
−13.3 kcal/mol in its own 1AQ1 and −6.8 to −9.7 everywhere else.

Two entries read `n/a`: LZ1's poses in 1HCK and 1FIN came back with a topology
Open Babel perceived differently from the reference, so no symmetry mapping
exists and no RMSD is claimed. Reporting `n/a` is the honest option; an
index-order RMSD would have produced a number that means nothing.

---

## Results

**Read the scale caveat first.** 60 compounds across 26 scaffolds went through the
docking ensemble; 65 through Boltz-2; 37 through both. Spearman's standard error
is ≈ 0.13 at n = 60 and ≈ 0.17 at n = 37. Differences smaller than about 0.3 are
not resolvable here. Nothing below establishes a ranking between methods — it
establishes what the pipeline measures, and which claims survive their own
controls.

### Question 2: does docking alone rank compounds?

Zero-parameter baselines, raw score as the ranker, no model fitted:

| slice | Spearman ρ | EF5% | what it is |
|---|---|---|---|
| 1AQ1 | **+0.351** | 2.31 | staurosporine-expanded pocket |
| 1FIN | +0.331 | 1.54 | cyclin-A-bound active form |
| 1HCK | +0.225 | 1.54 | ATP ground state (reference) |
| 1KE5 | +0.060 | 1.54 | oxindole |
| 2VTA | +0.054 | 0.77 | fragment-bound |
| 3PXF | **−0.115** | 1.54 | ATP site empty — the negative control |
| best-of-ensemble | +0.191 | 2.31 | |

Two things fall out that were not put in by hand. **3PXF ranks worst and
negative** — the slice whose ATP site is empty behaves like the negative control
it was included as, which is a check on the whole setup. And
**best-of-ensemble (+0.191) is worse than its best single slice (+0.351)**:
taking the most favourable score across an ensemble is a common heuristic, and
here it is dominated by whichever receptor is most permissive rather than most
relevant.

### Question 3: do 3D features add anything over 2D?

Scaffold-grouped 5-fold CV, 5 partitions, mean ± sd:

| arm | Spearman ρ | R² | RMSE | clears its null? |
|---|---|---|---|---|
| 2d (ECFP4 + descriptors) | +0.263 ± 0.047 | −0.019 | 1.252 | yes (null 95th pct +0.124) |
| 2d + ifp_ens | +0.275 ± 0.096 | +0.007 | 1.236 | — |
| ifp_ens | +0.098 ± 0.122 | −0.095 | 1.297 | **no** (null 95th pct +0.216) |
| ifp_ens + score | +0.059 ± 0.104 | −0.125 | 1.314 | — |
| ifp_1 | −0.022 ± 0.069 | −0.221 | 1.370 | no |
| score_ens | −0.049 ± 0.101 | −0.528 | 1.532 | no |
| score_1 | −0.299 ± 0.062 | −0.472 | 1.505 | no |

**The ensemble interaction fingerprint does not clear its own permutation null.**
At this scale it is indistinguishable from shuffled labels. Adding it to the 2D
baseline moves ρ from +0.263 to +0.275 — well inside the ±0.096 spread, so no
incremental information is demonstrated.

**Fitting a model on docking scores is worse than not fitting one.** `score_1`
reaches ρ = −0.299 where the same 1HCK score used raw gives +0.225. With 48
training compounds and one feature, the learner fits fold-specific structure and
inverts the sign. This is a result about small-N model fitting, not about
docking.

Note also that RMSE (1.24–1.53) sits above the 0.30-log assay noise floor by a
wide margin, so none of these arms is near the ceiling the data allows.

### What the random split would have told you instead

| arm | scaffold split | random split | inflation |
|---|---|---|---|
| 2d | +0.263 | +0.723 | 2.7× |
| ifp_ens | +0.098 | +0.437 | 4.5× |
| 2d + ifp_ens | +0.275 | +0.608 | 2.2× |
| score_ens | −0.049 | +0.331 | sign flip |

The arm that gains most from leakage is the one that fails its null under a
scaffold split. A random-split report of this same pipeline would have read as a
clear success for ensemble docking.

### Boltz-2 as a first-pass filter

Zero-shot, no training, so no split is needed — which makes this a fairer
comparison than the fitted arms above, not a worse one.

| score | ρ vs pIC50 | p | n |
|---|---|---|---|
| `optimization_score` | **+0.553** | 1.8e-06 | 65 |
| `binding_confidence` | +0.343 | 0.005 | 65 |
| `iptm` | +0.338 | 0.006 | 65 |
| `structure_confidence` | +0.124 | 0.33 | 65 |
| `complex_plddt` | +0.030 | 0.81 | 65 |

On the 37 compounds that went through **both** paths:

| | ρ | p |
|---|---|---|
| Boltz `iptm` | **+0.594** | 1.1e-04 |
| Boltz `binding_confidence` | +0.541 | 5.4e-04 |
| docking, 1AQ1 | +0.384 | 0.019 |
| docking, 1FIN | +0.379 | 0.021 |
| docking, best-of-ensemble | +0.138 | 0.42 |
| docking, 3PXF | −0.253 | 0.13 |

Boltz-2's interface confidence ranks CDK2 potency above every docking arm on the
same compounds. At n = 37 the gap (+0.594 vs +0.384) is about 1.2 standard
errors, so this is suggestive, not established.

The internal control matters more than the headline: **`structure_confidence`
and `complex_plddt` carry no signal (ρ = +0.12, +0.03) while the interface
scores do.** If the correlation were an artifact of easy-to-fold or
easy-to-dock chemistry, overall confidence would track potency too. It does
not, which is what you would want if the signal really is about the interface.

Three caveats that are not optional:

* **These are confidence scores, not affinities.** The library-screen endpoint
  returns `iptm` / `binding_confidence` / `optimization_score` and ADME. There
  is no calibrated affinity or predicted pIC50 in it. "Ranks potency" is the
  claim; "predicts potency" is not.
* **The default SMARTS alert filter silently dropped 35 of 100 compounds.** A
  third of a literature CDK2 set never gets scored unless you disable it. The
  dropped set is not potency-biased (Mann-Whitney p = 0.179, tested rather than
  assumed), so the comparison stands — but anyone running a real library through
  this endpoint should know the denominator changes underneath them.
* **`optimization_score` ranks best on 65 compounds and `iptm` best on the
  37-compound subset.** The ordering is not stable at this n; do not read a
  winner off it.

### Boltz-2 structures as ensemble slices

Predicted CDK2–ATP (`sab_pred_6eaMQZNUJ9s5FuoCxiL7`, structure_confidence 0.932,
iptm 0.952) against the 1HCK crystal, on the contact-shell residues:

* pocket CA RMSD **0.30 Å**
* global CA RMSD 0.51 Å

For scale, the crystal slices' pocket RMSD to the same reference run 0.44 (2VTA)
to 2.01 Å (1FIN). **The prediction is closer to 1HCK than any crystal form is** —
Boltz reproduces the ATP-bound pocket essentially exactly.

That is an accuracy result and a diversity problem at once. An ensemble is worth
something because its members differ; a predictor that lands on the canonical
structure adds a slice that is nearly redundant with the reference. Whether
co-folding a *different* ligand recovers genuine induced-fit spread — the
staurosporine-expanded pocket is the case that would show it — was **not
tested**. Those two predictions completed and their IDs are in
`results/boltz/job_ids.json`; the comparison is one download and one RMSD call,
and it is the first thing to run next.

### So, does the structural ensemble buy anything?

At validation scale, on this target: **not demonstrably.** The ensemble
interaction fingerprint fails its own permutation null, and fusing it with ECFP4
does not move the 2D baseline outside its error bars. The two structure-based
things that did carry signal were both zero-parameter: a single well-chosen
crystal conformer used as a raw score (ρ +0.35), and Boltz-2's interface
confidence (ρ +0.54–0.59).

What this does **not** say is that 3D features are useless — n = 60 with 26
scaffolds cannot support that conclusion, and the honest statement is that the
signal was not detected rather than that it is absent. What it does say is that
the benchmark is capable of returning a negative answer, which is the property
that makes a positive answer from it worth anything.

---

## Honesty machinery

These are the parts that exist specifically to stop the benchmark from
flattering itself.

**Scaffold-grouped CV is the headline; random CV is run alongside to show the
leak.** Analogue series split across train and test inflate every metric, and a
random split reports that inflation as skill. Both are computed; only the
scaffold split is quoted as a result.

**Compound selection is label-blind and keeps series intact.** A first attempt
sampled one compound per scaffold round-robin; with 935 scaffolds available the
quota filled on the first pass, so every selected compound had a unique
scaffold — which makes a scaffold split indistinguishable from a random split
and erases the within-series SAR that QSAR is supposed to learn. The benchmark
would have been measuring nothing. Whole scaffold groups are now taken, capped
at 25, and the selected set's potency distribution matches the source pool
(6.52 ± 1.13 vs 6.59 ± 1.14) because potency is never consulted during
selection.

**Zero-parameter baselines are reported before any model is fitted.** Docking
score as a ranker, per slice and ensemble-best. If a fitted model cannot beat
the raw score, that needs to be visible.

**A label-permutation null runs through the identical pipeline.** An arm is
interesting only if it clears its own null, not zero.

**The assay noise floor sits next to the model RMSE**, so "good" is judged
against what the data can support rather than against zero.

**Interaction channels are named for what they are.** The per-residue features
are `mind` / `ncon` / `apol` / `polr` — minimum distance, contact count, apolar
pairs, N/O-to-N/O pairs. `polr` is a distance proxy for a polar contact, **not**
a hydrogen bond: the receptors carry no explicit hydrogens and no per-residue
protonation states, so donor/acceptor roles and H-bond geometry are simply not
determined here. An aromatic-stacking channel was considered and dropped rather
than reported without the ring-normal angle that would make it meaningful.

---

## Things that were measured rather than assumed

Recorded because each one changed a decision, and because the wrong version of
each is a common way these pipelines go quietly wrong.

* **Shrinking the docking box does not speed docking up.** A 0.85× box was
  *slower* (45.5 s vs 39.7 s per dock): smina's step count is set by
  exhaustiveness, not box volume, and a tighter box rejects more moves.
  Exhaustiveness is the only real lever, and it trades directly against pose
  quality — which is what the interaction fingerprints are built from.
* **RMSD needs a symmetry-aware reference, and the reference must be the docking
  input.** Using the RCSB ideal-geometry SDF as a bond-order template returned
  `n/a` for every ligand, because Open Babel adds polar hydrogens writing PDBQT
  and RDKit silently keeps them when `sanitize=False`. Index-order RMSD is also
  wrong: LS1 reads 1.23 Å by atom order and 0.74 Å symmetry-corrected, because
  the sulfonamide-phenyl flip is an equivalence, not a displacement.
* **An early ATP "redocking failure" at 7.76 Å was an artifact** of that broken
  reference, not a docking failure. The corrected value is 0.60 Å.

---

## Scale, and what that means for the numbers

The full cross product is 784 compounds × 6 slices = 4,704 docking runs. On the
4 cores available here that is ~19 h, so **this run is validation-scale: 60
compounds × 6 slices = 360 docking runs**, with Boltz-2 on the first 100
compounds (65 surviving its filter, 37 overlapping the docked set).

The cross-docking and alignment controls do not depend on N and are fully valid.
The model comparison does depend on N, and at n = 60 / 37 it is underpowered by
design — see the scale caveat at the top of Results.

`s05` is checkpointed through a JSONL ledger and queues all slices of a ligand
together, so scaling up is `--max-ligands 784` and a re-run: everything already
computed is reused and the matrix is never left ragged. Poses are written as
plain appended SDF rather than gzip, because appending to gzip corrupts silently
across a restart — that bug cost a full re-dock and is documented in `s05`.

---

## Cost

| item | unit | this run |
|---|---|---|
| Boltz-2 library screen | $0.025 / compound | 100 submitted = $2.50 |
| Boltz-2 structure + binding | $0.02 / sample | 3 × 5 samples = $0.30 |
| Docking, dataset, features | free (local, smina + RDKit) | — |
| | | **$2.80** |

Job IDs and inputs are recorded in `results/boltz/job_ids.json`.

---

## Layout

```
data/        dataset, screening subset, 3D ligands (ligands/ gitignored)
structures/  PDB slices, superposed receptors, ensemble_manifest.json
pipeline/    s01..s08, one stage per file, each runnable alone
results/     controls, scores, fingerprints, model comparison
```
