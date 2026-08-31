# OpenKnot AI reproduction — does the in-silico screening score predict the wet lab?

A CPU-only, from-scratch reproduction of the **computational half** of

> Townley J, Kladwang W, Baker D, … Das R.
> *De novo design of RNA pseudoknots with deep learning.* **Science** 393:931–937 (2026).
> [10.1126/science.aeg6829](https://doi.org/10.1126/science.aeg6829) ·
> preprint [10.64898/2026.05.21.726960](https://doi.org/10.64898/2026.05.21.726960)

The paper's central computational claim is that RNA pseudoknot design can be driven by a
chemical-mapping foundation model (**RNet** = RibonanzaNet) *instead of* 3D structure
prediction. That claim rests on one link: **the in-silico OpenKnot score must predict the
experimental one.** Everything here tests that link against the authors' own released
experimental data.

## What was run

| Stage | What | Scale |
|---|---|---|
| Setup | gRNAde + RibonanzaNet + RibonanzaNet-SS, CPU only, no GPU | 4 cores, 15 GB RAM |
| Validation | re-implemented OpenKnot score vs. the published `target_openknot_score` | 17,710 designs |
| **Main test** | RibonanzaNet → predicted 2A3 SHAPE → in-silico OpenKnot score, compared to measured | **17,710 designs** (Rounds 1 + 3) |
| Design | gRNAde 2D-mode inverse folding for target P20 "Kissing Multiloops" | 1,024 designs |

Total compute: ~2.5 CPU-hours. No GPU was used or needed.

## Headline results

See `results/analysis_output.txt` for the full run and `figures/` for the plots.

**The in-silico score is a real but weak-to-moderate predictor.**

| | Round 1 (designers had **no** access to RNet) | Round 3 (designs pre-filtered on RNet) |
|---|---|---|
| n (SN_filter=1) | 8,571 | 8,946 |
| experimental success rate (OpenKnot > 90) | 19.9% | 39.0% |
| pooled Spearman (in-silico vs measured) | 0.60 | 0.40 |
| **within-target Spearman (median)** | **0.46** | **0.30** |
| AUROC, in-silico OpenKnot score | **0.749** | **0.671** |
| AUROC, published `RNet_F1` filter | 0.720 | 0.570 |
| AUROC, GC content alone (baseline) | 0.593 | 0.449 |

Round 3 looks worse than Round 1 only because of **selection bias**: Round 3 designs were
already screened on RNet before submission (7,437 / 8,946 have `RNet_F1` ≥ 0.8), which
compresses the range the score has left to discriminate over. On Rosetta designs — the one
method that did *not* use RNet — the `RNet_F1` ≥ 0.8 filter lifts experimental success from
**10.5% → 28.1%**.

**Screening enrichment** (pick top-k per target by in-silico score, then look at what the
wet lab measured):

| | Round 1 | Round 3 |
|---|---|---|
| base rate | 19.9% | 39.0% |
| top-1 | 47.1% (×2.4) | 65.0% (×1.7) |
| top-10 | 43.5% (×2.2) | 55.0% (×1.4) |
| targets with ≥1 success in top-10 | 88% (random-10: 65% ± 10%) | 90% (random-10: 81% ± 6%) |

**It is not just GC content.** GC is a genuine confound — GC-rich designs both look
unreactive to the model *and* actually fold better (measured Spearman(OpenKnot, GC) = 0.39
on target P20). But holding GC roughly constant inside narrow strata, the in-silico score
still separates successes from failures (AUROC 0.53–0.82 across strata, mostly 0.65–0.84).
The model is contributing beyond the trivial baseline.

**Per-design SHAPE profile accuracy** is modest: mean Pearson *r* = 0.50 between predicted
and measured 2A3 reactivity (median 0.51; only 10% of designs exceed r = 0.7). The
aggregate OpenKnot score is a far better signal than the profile it is computed from.

## Correctness checkpoints

Two things were verified rather than assumed:

1. **The scorer.** Feeding *experimental* SHAPE into the re-implemented OpenKnot score
   reproduces the published `target_openknot_score` with mean |error| = 0.48 points on a
   0–100 scale (80% of designs within 1 point; Pearson r = 0.995 Round 1, 0.872 Round 3 —
   the lower r reflects Round 3's compressed score range, not larger errors). A residual
   scorer-detail difference remains for ~20% of designs; it does not affect any conclusion
   here, because every comparison uses the *published* experimental score as ground truth.
2. **The output channel.** RibonanzaNet emits two channels. Channel 0 correlates with
   measured 2A3 at mean r = 0.59 vs 0.43 for channel 1, confirming gRNAde's `[:, :, 0]`
   indexing selects 2A3.

## Dependency shims (`shim/`)

`torch_cluster` and `torch_scatter` have no prebuilt wheels for torch 2.13 CPU. Both are
replaced by pure-PyTorch equivalents, each proven safe:

- `scatter_add` — matches a reference implementation bit-for-bit on 200 random 1-D and 2-D cases.
- `knn_graph` — in gRNAde's 2D design mode the featurizer output is **bit-identical** even
  when `knn_graph` is sabotaged to return a random graph, because all 3-D edges are masked
  out. The shim is provably inert on this path.

## Caveats

- **Training-set overlap is not fully excluded.** RibonanzaNet was trained on ~1M Eterna
  chemical-mapping measurements that include earlier OpenKnot rounds. Round 1 designs here
  were probed later and are the cleaner test, but sequence-level overlap was not audited.
- **This reproduces the screening step, not the paper.** The design methods' own success
  rates depend on generation *and* filtering; only the filtering link is tested here.
- The gRNAde P20 run is 1,024 designs against the paper's ~1,000,000, and the top design
  was never synthesised — its 95.8 is an in-silico number with no experimental status.
- 3 of 36,761 designs carry a malformed target dot-bracket in the released data
  (`W14`, 5 `[` vs 4 `]`) and are skipped.

## Reproducing

```sh
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install numpy pandas scipy pyyaml einops tqdm matplotlib torch_geometric \
            torchmetrics python-dotenv ml_collections biopython arnie wandb \
            biotite cpdb-protein draw_rna MDAnalysis

git clone https://github.com/chaitjo/geometric-rna-design.git gRNAde
# checkpoints (RibonanzaNet, RibonanzaNet-SS, gRNAde) — public, no auth
hf download chaitjo/gRNAde --local-dir gRNAde/checkpoints/

# experimental data (Git LFS)
git clone https://github.com/eternagame/OpenKnotAIDesignData.git data-okb
curl -L -o data-okb/Data/OpenKnotBench_data.v4.5.1.csv \
  https://media.githubusercontent.com/media/eternagame/OpenKnotAIDesignData/main/Data/OpenKnotBench_data.v4.5.1.csv

python scripts/extract.py                        # SHAPE + metadata -> okb_designs.pkl
python scripts/run_rnet.py --rounds 1 3          # ~50 min on 4 CPU cores
python scripts/analyze.py
python scripts/figures.py
python scripts/design_grnade.py --puzzle P20 --total 1024   # gRNAde inverse folding
```

Paths are absolute in the scripts (`/home/user/work/...`) — adjust for your layout.

## Sources

- Code archive (MIT): [Zenodo 20649966](https://zenodo.org/records/20649966)
- Experimental data: [eternagame/OpenKnotAIDesignData](https://github.com/eternagame/OpenKnotAIDesignData), RMDB `OK45LIB_2A3_0000` / `OK7ALIB_2A3_0000`
- gRNAde: [chaitjo/geometric-rna-design](https://github.com/chaitjo/geometric-rna-design)
- RibonanzaNet: [Shujun-He/RibonanzaNet](https://github.com/Shujun-He/RibonanzaNet)
