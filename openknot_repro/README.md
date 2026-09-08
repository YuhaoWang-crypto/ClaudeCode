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
| Design generation with gRNAde / NA-MPNN / Struct2SeQ, and RNet inference | not here — see *Not done yet* |

## Results

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

## Layout

```
openknot/score.py         OpenKnot score, no dependencies beyond the stdlib
openknot/data.py          loaders for the released CSVs
scripts/download_data.sh  fetch the release (Git LFS media endpoint, ~200 MB)
scripts/validate_score.py recompute every score and account for every difference
scripts/success_rates.py  per-method and AI-vs-human success rates
scripts/figures.py        the three figures
tests/test_score.py       unit tests plus a regression against released scores
```

## Running it

```bash
pip install pandas numpy matplotlib pytest
bash openknot_repro/scripts/download_data.sh
python openknot_repro/scripts/validate_score.py     # ~30 s
python openknot_repro/scripts/success_rates.py
python openknot_repro/scripts/figures.py
python -m pytest openknot_repro/tests -q
```

`data/` is not committed (200 MB); the download script recreates it. Set
`OPENKNOT_DATA` to put it elsewhere.

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

It does not show that any new design would work. Nothing here generates
sequences, and scoring a design needs its measured SHAPE profile — which needs
the wet lab. A predicted profile from RNet can stand in, but designing with RNet
and then scoring with RNet is a closed loop that will flatter itself; the paper's
value is precisely that its judge was an experiment.

## Not done yet

* RNet (RibonanzaNet) inference, so predicted SHAPE profiles and predicted
  secondary structures can be computed for new sequences rather than read from
  the release. Weights live on Kaggle and need an account.
* The design methods themselves — gRNAde, NA-MPNN, Struct2SeQ — all open source,
  all needing a GPU.
* M2R-seq stem-recovery analysis (Fig. 3G) from `OK7a_M2R_data.v4.5.1.csv`, which
  is downloaded but not yet analysed.
