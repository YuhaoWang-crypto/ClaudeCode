# Technologies referenced by the Virtual Cell Challenge paper

Source: Roohani, Hua, Tung, … Goodarzi & Burke, *"Virtual Cell Challenge: Toward a
Turing test for the virtual cell"*, **Cell 188:3370–3374** (26 June 2025),
[doi:10.1016/j.cell.2025.06.008](https://doi.org/10.1016/j.cell.2025.06.008).

Each row records what the technology provides, and — because the brief was to
actually obtain and use these — whether it was reachable from this machine and
what was done with it. "Used" means it is wired into the model or the benchmark
in this directory; "surveyed" means read and drawn on for design, not executed.

## 1. First-generation whole-cell simulators (refs 1–5)

| Technology | What it is | Status here |
|---|---|---|
| **E-CELL** (Tomita 1999) | Whole-cell simulation environment; hand-built ODE/stochastic models of a minimal cell | Surveyed. Mechanistic, parameter-hungry; the paper cites it as the approach that *stalled* for want of data |
| **Virtual Cell / VCell** (Loew & Schaff 2001; Slepchenko 2003) | Compartmental reaction-diffusion simulator for cell biology | Surveyed. Same lineage; the name "virtual cell" is inherited from here, the method is not |
| **Whole-cell simulation as a grand challenge** (Tomita 2001) | Position paper | Surveyed |
| **Lessons from dynamical modelling** (Raue 2013) | Identifiability/sloppiness analysis in systems biology | Surveyed. Directly relevant: it is why this model is deliberately low-dimensional and its parameters are cross-validated rather than fit |

The paper's framing is that these failed on data scarcity, and that the modern
substitute is statistical learning over large perturbation atlases. That framing
is adopted here.

## 2. The AI-virtual-cell roadmap (ref 6)

**Bunne et al., "How to build the virtual cell with AI: priorities and
opportunities", Cell 187:7045 (2024).** Defines the universal cell embedding /
virtual cell concept, the decomposition into cell-state representation +
perturbation response, and the evaluation gap this challenge exists to close.
Surveyed; its decomposition (context representation, perturbation
representation, response decoder) is the skeleton of `model.py`.

## 3. Perturbation data resources

| Resource | Content | Status here |
|---|---|---|
| **Replogle 2022 genome-scale Perturb-seq** (ref 11) | CRISPRi Perturb-seq: K562 genome-wide (~9,900 targets), K562 essential, RPE1 essential | **Used.** Downloaded pseudobulk from [Figshare+ 20029387](https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387). K562-essential (2,057 targets) and RPE1 (2,393 targets) are two of the four benchmark contexts |
| **Nadig 2025 perturbation atlases** (ref 16) | Transcriptome-wide DE across perturbation atlases; contributes matched essential-gene CRISPRi in **Jurkat** and **HepG2** | **Used.** Downloaded single-cell h5ad from GEO **GSE264667** (5.6 GB + 9.4 GB), pseudobulked here by `prep_nadig.py`. These are the other two benchmark contexts |
| **Arc Virtual Cell Atlas — scBaseCount** (ref 14) | AI-agent-curated, uniformly reprocessed scRNA-seq; ~230M+ cells | Surveyed. Observational only (no perturbations), so it cannot supply the transfer signal this task needs; would serve a pretrained cell-state encoder, which the 2025 results suggest is not where the gains are |
| **Arc Virtual Cell Atlas — Tahoe-100M** (ref 15) | ~100M cells, *chemical* perturbations across ~50 cell lines | Surveyed. Wrong perturbation modality for a gene-knockdown task, but it is the only public resource with enough *contexts* to fit a real context encoder |
| **scPerturb** (ref 7) | Harmonised single-cell perturbation datasets across labs and technologies | Surveyed. Harmonisation is exactly the problem; not used because the four lines adopted here already share one library design, one lab and one pipeline, which removes a batch confound rather than adding one |
| **McFaline-Figueroa 2024** (ref 17) | Multiplex single-cell chemical genomics, kinase dependence | Surveyed; chemical modality |
| **Jiang 2025** (ref 18) | Scalable single-cell perturbation screens → pathway signatures | Surveyed |
| **VCC 2025 H1 hESC dataset** | ~300k cells, 300 CRISPRi targets, 10x Flex | Not obtained: distributed through the challenge portal behind registration. Note that VCC-2 explicitly permits its reuse, so a real entry would add it as a fifth context |

## 4. Experimental technology

| Technology | Role |
|---|---|
| **CRISPRi** (Gilbert 2014, ref 13) | The perturbation modality: dCas9-KRAB knocks *down* rather than out. Measured here at 82–91% median residual knockdown across the four lines, and per-gene efficiency transfers between lines (Spearman 0.38–0.54) — which the model exploits |
| **10x Genomics Flex** | Probe-based fixed-cell chemistry used for the challenge dataset; cell fixation reduces batch effects and improves residual knockdown. The four public lines used here predate it and are 3'-capture, at ~45–120 cells and ~11–14k UMI per perturbation versus the challenge's ~1,000 cells and >50k UMI |
| **RNAi caveats** (Kaelin 2012, ref 10) | Cited as the standards argument: off-target and incomplete-knockdown artefacts are why reproducibility, not just scale, is the bottleneck |

## 5. Evaluation

| Technology | What it is | Status here |
|---|---|---|
| **cell-eval** (Arc) | The official metric implementation | **Used.** `cell_eval` 0.8.2 pulled from PyPI and read; `metrics.py` is a port of `discrimination_score`, `mae`/`mae_delta`, `pearson_delta`, `overlap_at_k`, `de_direction_match`, `de_spearman_lfc`, and the `_score.py` baseline-normalised aggregation |
| **PerturBench** (Wu 2024, ref 12) | Origin of the perturbation discrimination score | Surveyed via the cell-eval implementation, which is the operative definition |
| **NeurIPS-Kaggle 2023** (Szałata 2024, ref 9) | Prior benchmark, chemical perturbation in immune cells | Surveyed |
| **Cancer Immunotherapy Grand Challenge** (ref 8) | Prior benchmark, phenotype-level not expression-level | Surveyed |

The three headline metrics named in the commentary map onto the port as:
differential expression score → `overlap_at_k`; perturbation discrimination
score → `discrimination_score` (L1, target gene excluded); MAE → `mae`.

## 6. What the 2025 results changed about the design

From Arc's [2025 wrap-up](https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up):

- **1st, BioMap `xTrimoSCPerturb`** — deep learning *plus* classical statistics,
  protein embeddings, public perturbation data. Arc's stated conclusion:
  *"purely AI-based approaches did not consistently outperform statistical
  baselines."*
- **2nd, XLearning Lab** — pseudobulk representation, a plain fully-connected
  network, ESM-2 embeddings, public Perturb-seq only.
- **3rd, `TransPert`** — summary-level statistics, similarity-aware aggregation,
  global linear scaling for metric optimisation.
- **Generalist prize, Altos `go-with-the-flow`** — flow matching over
  heterogeneous single-cell responses.
- Most submissions were **worse than baseline on MAE**, so entrants optimised
  the discrimination and DE scores instead.

Three of the four winning designs are summary-level statistical transfer with
context conditioning, which is what `model.py` implements: similarity-aware
aggregation (3rd place's idea), pseudobulk representation (2nd place's), and an
explicit global scale (3rd place's). The concurrent
[*Virtual Cells Need Context, Not Just Scale*](https://pmc.ncbi.nlm.nih.gov/articles/PMC12919078/)
(2026) argues the same from the data side: DEG recovery improves with the number
of *contexts* a perturbation is observed in, not with cell count, and it
recommends making context explicit rather than inferring it from expression.
This benchmark is built to test exactly that claim, with four contexts.

## 7. Deliberately not used

- **scGPT / Geneformer / scFoundation and other single-cell foundation models.**
  No GPU on this machine, and the 2025 evidence is that they do not beat
  statistical baselines on this task. A pretrained gene embedding would be a
  cheap addition (2nd place used ESM-2) and is the most defensible next step.
- **GEARS / CPA / scLDM / STATE.** Surveyed as design references. GEARS' idea —
  reach unseen perturbations through a gene-similarity graph — appears here as
  the perturbation-neighbour smoothing that makes the double-blind regime
  possible at all.
