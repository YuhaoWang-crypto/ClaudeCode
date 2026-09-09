# OpenKnot AI pseudoknot design benchmark — computational reproduction

Reproduction of the scoring and benchmark analysis of *De novo design of RNA
pseudoknots with deep learning* (Townley, Kladwang, Das et al., bioRxiv 2026,
[10.64898/2026.05.21.726960](https://doi.org/10.64898/2026.05.21.726960)), from
the released experimental data.

This is the analysis half of the paper. The SHAPE chemical mapping, the M2R-seq
compensatory mutagenesis and the cryo-EM are wet-lab work and are **not**
reproduced here; the released measurements are used as ground truth.

## What is in here

| Piece | Status |
|---|---|
| OpenKnot score (Eterna Classic Score + Crossed Pair Quality), dependency-free port | reproduces the released score for **every** design outside two documented puzzles |
| Per-method, per-round success rates (the paper's Fig. 1H/1J, 2C/2F) | recomputed scores give **identical** rates to the released ones, to 0.000 pp |
| Targets solved by AI vs by Eterna participants | matches the paper's headline **19/20 and 19/20** for Rounds 3 and 4 |
| RibonanzaNet (RNet) inference: SHAPE reactivity and secondary structure | official checkpoints, reproduced without a Kaggle account; the paper's RNet F1 filter agrees on **100%** of designs |
| gRNAde design, re-run on all 40 Round 3 and 4 targets | 320 designs on CPU; the gap to the paper's designs is a **search-budget** gap |
| Struct2SeQ design | blocked — weights are Kaggle-only, see *Not done yet* |

## Results

### Scoring

Run over all 36,761 released designs (`results/score_validation.csv`):

```
designs scored                              36761
exact vs own target structure               34415   93.618%
explained (exact, or within rounding ties)  36143   98.319%
unexplained                                   618    1.681%
```

Every one of the 618 unexplained designs belongs to puzzle **W09** (rounds 1–2)
or **P06** (round 3) — the two puzzles the v4.5.0 release notes single out as
having been released with two slightly different target structures and rescored
against the better of the two. Outside those, **100% of designs are accounted
for**: 93.6% reproduce bit-exactly, and the remainder differ only by residues
whose released reactivity, rounded to three decimals, lands exactly on a scoring
threshold (0.125, 0.25 or 0.5) where the original scoring saw the unrounded
value. Those residues are counted per design, and in every case the published
score sits inside the interval they could produce.

Targets solved by at least one design scoring above 90 (`figures/fig_ai_vs_human.png`):

| Round | AI (any method) | Eterna | Starting sequence |
|---|---|---|---|
| 1 | 11/16 | 17/17 | 7/12 |
| 2 | 17/17 | 17/17 | 7/12 |
| 3 | **19/20** | **19/20** | 12/18 |
| 4 | **19/20** | **19/20** | **4**/14 |

The Round 3 and Round 4 figures are the paper's central claim (AI reaching
parity with experienced human designers), and the Round 4 baseline of 4 successes
matches the paper's "the starting sequences were successful in only 4 of 20
targets" — the denominator differs because starting sequences for six of the
targets have no design passing the release's signal-to-noise gate.

## RibonanzaNet, reproduced

RNet is the model that carries the paper's argument: it stands in for the 3D
structure prediction that RNA design does not have, both to guide the design
methods and to filter their output. Two things had to be established.

**The checkpoints.** They are distributed through Kaggle, which needs an
account. Three HuggingFace mirrors avoid that, and they agree with each other:
`roos23/RibonanzaNet` carries the base model verbatim; `chaitjo/gRNAde` carries
both official checkpoints alongside gRNAde; and `multimolecule/ribonanzanet-ss`
carries the secondary-structure fine-tune under a different naming scheme.
`scripts/convert_rnet_weights.py` derives the rename rule between the two
layouts, requires it to map all 576 base tensors bit-identically, applies it to
the SS checkpoint, and then checks the result against the official file — all
578 tensors match. So the model here is the paper's model, not an approximation
of it.

**The inference.** Run on 399 released designs, stratified across the four
rounds (`results/rnet_validation.csv`):

| Quantity | Result |
|---|---|
| secondary structure vs the released `RNet_structure` | mean base-pair F1 **0.999**, identical pair sets on **96.9%** of designs |
| RNet F1 against the target vs the released `RNet_F1` | Spearman **0.982**, mean absolute difference **0.002** |
| the paper's filter decision at RNet F1 ≥ 0.8 | **100%** agreement (306 kept by both, 93 dropped by both) |

RNet must be run on the design sequence alone. The released columns are
design-length, and predicting on the padded construct and slicing agrees
distinctly worse.

### Scoring a design without an experiment

The main text says RNet "gave simulated scores largely reproducing experimental
scores, especially for poorly performing designs", which is what licenses
filtering designs before spending an experiment on them. Computing the OpenKnot
score from RNet-predicted reactivity instead of measured reactivity, on the same
399 designs:

- Spearman **0.56** against the experimental score (Pearson 0.49; 0.59 on the
  designs that pass the release's signal-to-noise gate).
- As a prospective filter at the cutoff of 90: precision **0.84**, recall
  **0.83**, against a base rate of 0.73.
- The error is **not** uniform, and not in the direction the sentence suggests
  to a first reading. On designs that measured below 70, the simulated score is
  **+17.4 points too generous** on average; between 80 and 90 the mean signed
  error is +1.3. RNet is least reliable exactly on the designs that failed.

The supplementary figure this refers to (Fig. S4) is not in the PDF used here,
so this tests the main-text sentence rather than the figure. Two caveats also
cut in the paper's favour: the experiment measures the design inside its
flanking pads while RNet here sees the design alone
(`scripts/rnet_padded_check.py` tests whether that explains the gap), and this
is a 399-design sample rather than the full release.

Either way the practical consequence stands: an in-silico score can rank and
pre-filter designs, but it cannot stand in for the measurement — which is what
the paper spent 50,000 experiments to establish.

## Designing with gRNAde

gRNAde runs here on CPU against the paper's own target set — the Round 3 and
Round 4 metadata and structures ship inside the gRNAde repository, under
`projects/openknot_benchmark`. 8 samples per target (4 each at temperature 0.1
and 0.5), in 2D mode, which conditions on the target secondary structure alone
and is the paper's `gRNAde-no3d` variant. 320 designs across 40 targets, about
40 minutes on four cores.

Those designs can only be scored in silico, so comparing them against the
paper's *measured* scores would confound two different things. The comparison in
`scripts/compare_designs.py` separates them by scoring the paper's own submitted
gRNAde design the same in-silico way (`figures/fig_grnade_designs.png`):

|  | mean simulated score | above 90 |
|---|---|---|
| best of 8 samples per target, here | 89.4 | 21/40 targets |
| the paper's submitted design, simulated | 93.6 | 35/40 targets |
| the paper's submitted design, **measured** | — | **27/40 targets** |

Two things fall out of that table.

**The gap to the paper's designs is a search-budget gap.** Same model, same
targets, same scorer: 4.2 points of simulated score separate 8 samples from a
search of up to a million. Our best beat theirs on 5 of 40 targets, which is
about what a tiny sample should manage against a large one.

**The simulation is optimistic, and by a measurable amount.** For the very same
molecules, the in-silico score clears 90 on 35 targets while the experiment
clears it on 27. That is the same bias the RNet section quantifies, seen from the
design side: a design pipeline scored only by RNet will believe it has succeeded
about a third more often than the experiment agrees.

## Getting the OpenKnot score right

The score is `0.5 × ECS + 0.5 × CPQ`, and four details decide whether a
re-implementation matches. All four were established by testing against the
released scores, not assumed:

1. **The ECS paired threshold is 0.5, not 0.25.** The dataset README says target
   paired residues "should have reactivity less than 0.25"; 0.25 is the CPQ
   threshold. The reference code uses 0.5 for ECS, and 0.5 is what reproduces the
   released column: on the first 6,000 designs, 0.5 matches 89.1% exactly and
   0.25 matches 5.5%.
2. **Missing reactivity counts as a miss in ECS.** The reference skips NaN
   positions in its loop but leaves them in the denominator.
3. **CPQ is computed after singlet filtering; ECS is not.** Helices of length 1
   are dropped before crossed pairs are enumerated for CPQ. This is the
   `filter_singlets` flag, which defaults to False in the reference function but
   was evidently on for this release. It matters on exactly the puzzles whose
   targets contain lone pairs — turning it off leaves W03, W04, W08, P01, Q18 and
   others systematically wrong while the simpler puzzles still look correct, so
   spot-checking one puzzle will not catch it.
4. **The base-pair list must be sorted by opening index before helices are
   grouped.** arnie's `get_helices` walks the list and only ever stacks a pair
   onto the previous one, so an unsorted list shatters every helix into singlets
   and singlet filtering then deletes the entire structure.

Point 3 is corroborated independently: gRNAde ships its own adaptation of the
same pipeline (`src/openknot_score.py`), and it too filters singlet helices
before enumerating crossed pairs while scoring the Eterna half on the
unfiltered structure. On 200 released designs that implementation reproduces
the published score exactly 129 times against this one's 189.

## Layout

```
openknot/score.py             OpenKnot score, no dependencies beyond the stdlib
openknot/data.py              loaders for the released CSVs
openknot/rnet.py              RibonanzaNet: reactivity and secondary structure
openknot/grnade.py            gRNAde: sample designs for a target
scripts/download_data.sh      fetch the release (Git LFS media endpoint, ~200 MB)
scripts/setup_rnet.sh         fetch RNet code and checkpoints, verify the conversion
scripts/setup_grnade.sh       fetch gRNAde and its checkpoint
scripts/validate_score.py     recompute every score, account for every difference
scripts/success_rates.py      per-method and AI-vs-human success rates
scripts/rnet_predict.py       run RNet over released designs
scripts/validate_rnet.py      check RNet against the release, test the Fig. S4 claim
scripts/rnet_padded_check.py  does the flanking pad explain the simulated-score gap?
scripts/grnade_design.py      design for the Round 3/4 targets, score with RNet
scripts/compare_designs.py    our designs vs the paper's, scored the same way
scripts/figures.py            the figures
tests/                        unit tests, plus regressions against the release
```

## Running it

Analysis of the release only:

```bash
pip install pandas numpy matplotlib pytest
bash openknot_repro/scripts/download_data.sh
python openknot_repro/scripts/validate_score.py     # ~30 s
python openknot_repro/scripts/success_rates.py
python openknot_repro/scripts/figures.py
```

Running the models as well (CPU is enough; no GPU is used anywhere here):

```bash
pip install torch arnie==0.1.9 scipy einops     # RNet
pip install torch_geometric torch_scatter torch_cluster \
            biotite biopython python-dotenv wandb cpdb-protein h5py   # gRNAde
bash openknot_repro/scripts/setup_rnet.sh
bash openknot_repro/scripts/setup_grnade.sh
python openknot_repro/scripts/rnet_predict.py --per-round 100   # ~20 min
python openknot_repro/scripts/validate_rnet.py
python openknot_repro/scripts/grnade_design.py --rounds 3 4 --samples 4  # ~30 min
python openknot_repro/scripts/compare_designs.py
python -m pytest openknot_repro/tests -q
```

`data/`, `weights/` and `third_party/` are not committed; the setup scripts
recreate them. Set `OPENKNOT_DATA` to move the data elsewhere. Long runs append
to their output CSV and skip what is already there, so they can be interrupted
and restarted.

## Data sources

* Designs, SHAPE profiles, scores, M2R and M2 data:
  [eternagame/OpenKnotAIDesignData](https://github.com/eternagame/OpenKnotAIDesignData)
  v4.5.2, archived at [Zenodo 10.5281/zenodo.20101887](https://doi.org/10.5281/zenodo.20101887).
  The Zenodo archive is a GitHub source snapshot, so its CSVs are Git LFS
  pointers — the download script fetches the real files.
* Reference scoring implementation:
  [eternagame/OpenKnotScorePipeline](https://github.com/eternagame/OpenKnotScorePipeline)
  (`src/openknotscore/pipeline/scoring.py`) and [DasLab/arnie](https://github.com/DasLab/arnie)
  (`convert_dotbracket_to_bp_list`, `get_helices`, `post_process_struct`).
* Raw chemical mapping: RMDB accessions OK45LIB, OK6LIB, OK7ALIB, OK7BLIB.

## What this does and does not show

It shows that the paper's experimental scoring and its benchmark conclusions
follow from the released data, and it gives a scorer that can be pointed at new
designs.

It now also generates designs, with the same model and on the same targets as
one of the paper's three AI methods.

What it cannot do is tell you whether a new design works. Scoring a design needs
its measured SHAPE profile, which needs the wet lab. RNet's predicted profile
stands in, but designing under RNet's guidance and then scoring with RNet is a
closed loop that flatters itself — and the numbers above put a size on how much:
the simulated score is systematically too generous on designs that fail. The
paper's value is precisely that its judge was an experiment.

## Not done yet

* **Struct2SeQ**, the third AI method and the strongest one in Round 4. Its code
  is open, but `Struct2SeQ.pt` and `Struct2SeQ_SHAPE.pt` are published only as a
  Kaggle dataset, and Kaggle's API needs an account — anonymous requests get a
  browser challenge. Two ways to unblock it: put `KAGGLE_USERNAME` and
  `KAGGLE_KEY` in the environment, or drop the two `.pt` files into `weights/`.
  Everything else it needs (RibonanzaNet.pt, RibonanzaNet-SS.pt) is already here.
* **NA-MPNN**, the Baker-lab method behind MPNN-fixbb and MPNN-RFdiff. Open
  source, not yet wired up.
* **A real sampling budget.** The paper searched up to a million designs per
  target; this samples a few dozen on four CPU cores. Modal is not reachable
  from this environment — its client speaks gRPC, which the network policy does
  not carry — so there is no GPU here.
* **M2R-seq stem-recovery analysis** (Fig. 3G) from `OK7a_M2R_data.v4.5.1.csv`,
  which is downloaded but not yet analysed.
