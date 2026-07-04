# lipidlib — AI-guided lipid & targeting-ligand library screening for LNP delivery

An end-to-end, reproducible pipeline for using machine learning to **design and
screen molecular libraries for LNP (lipid nanoparticle) mRNA/RNA delivery**, built
around the two AI-lipid papers (Witten 2024 "LiON"; Su 2026 spatial-conformation)
and public predicted-interaction databases.

Everything here is scripted, committed, and figured. This file is the map; each
step links to a focused doc under [`docs/`](docs/).

---

## The one idea that organises everything: two separate problems

The project deliberately splits what is often conflated (full rationale in
[`docs/PLAN.md`](docs/PLAN.md)):

| | **Track A — ionizable-lipid design** | **Track B — active-targeting ligand** |
|---|---|---|
| optimize | the lipid *inside* the LNP | a ligand *conjugated to the LNP surface* (via PEG) |
| governs | encapsulation, endosomal escape, **passive** organ tropism (protein corona) | **active** receptor binding → receptor-mediated endocytosis |
| model | structure→property QSAR (graph NN) | protein–ligand binding / retrieval / docking |
| basis | **LiON** (Witten 2024) & Su 2026 | **DrugCLIP**, humanPPI, ChEMBL, Boltz |

Both attached papers are Track A. The "small molecule binds a CD ectodomain →
conjugate to the liposome" idea is Track B. The ionizable lipid is **not** the
recognition element and is **not** extracted from the particle.

---

## Results at a glance

| # | What | Result | Doc / figure |
|---|---|---|---|
| A1 | Reproduce LiON on Modal | lite (T4, 1-fold) RMSE 0.90 in 6 min; full 5-fold ensemble RMSE 0.78–0.86 | [PLAN](docs/PLAN.md) |
| A2 | Enumerate + screen combinatorial library (323 lipids) | model **recovers literature SAR**: unsaturated tails win, 3-tail optimal | [COMBO_LIBRARY](docs/COMBO_LIBRARY.md) · `combo_v1_screen.png` |
| A3 | 5-fold ensemble re-rank + confidence | reranks vs lite (ρ=0.78); confidence-aware leads = ester-linoleyl H7/H10 | `combo_v1_ensemble.png` |
| A4 | Organ re-screen (lung/liver/spleen) | SAR shifts: spleen→2-tail/small-heads; unsat disfavoured in vivo | [ORGAN_SCREEN](docs/ORGAN_SCREEN.md) · `combo_v1_organ.png` |
| A5 | Expanded library (538) + spleen deep-screen | degradable S–S/thioether tails **help spleen**; leads shortlisted | [COMBO_V2_SPLEEN](docs/COMBO_V2_SPLEEN.md) · `combo_v2_spleen.png` |
| B1 | GLP1R pilot: mine 3 sources | ChEMBL 1,422 · DrugCLIP 173 · humanPPI 259 | [PROBLEM_B_GLP1R](docs/PROBLEM_B_GLP1R.md) |
| B2 | Cross-compare ChEMBL × DrugCLIP | **no consensus** (fragments ≠ drug-sized actives), triangulated | [GLP1R_crosscompare](docs/GLP1R_crosscompare.md) · `glp1r_crosscompare.png` |
| B3 | humanPPI protein-layer enrichment | GLP1R interactome **cell-surface enriched** (p<1e-9); leads DLK1, PAM | `glp1r_humanppi_enrichment.png` |
| B4 | Generic pipeline + validation | any-target one-command; enrichment machinery detects real similarity at **1674×** | [PIPELINE](docs/PIPELINE.md) · `pipeline_validation.png` |
| M1 | Mechanistic delivery kinetics (Phase 4) | uptake→escape→translation ODE; ~1.8% escape bottleneck, ~13 h peak; links LiON potency → expression dynamics | [DELIVERY_KINETICS](docs/DELIVERY_KINETICS.md) · `delivery_kinetics.png` |

---

## Track A — ionizable-lipid design (the two papers)

Flow: **featurize → reproduce LiON → enumerate library → screen → rank → per-organ / expand**

1. **Featurize** any SMILES: `lipidlib/featurize.py` (Morgan/MACCS/RDKit-2D) + `lipidlib/cli.py`.
2. **Reproduce LiON** (Witten 2024): `scripts/fetch_lnp_ml.sh` clones the MIT repo
   (13,331 datapoints); `modal_app/lion.py` runs its Chemprop D-MPNN on Modal GPU
   (Python 3.8 in a `uv` venv; `cv_num` made configurable). Two tiers:
   `train_lite` (T4, 1-fold, 6 min) · `train` / `train_resilient` (A10G, 5-fold,
   preemption-safe).
3. **Enumerate** a combinatorial library by aza-Michael chemistry:
   `lipidlib/lion_library.py::enumerate_michael_lipids` (amine head N-H + acrylate/
   acrylamide tail → ester/amide-linked lipid). Drivers: `analysis/enumerate_library.py`
   (v1, 323) and `analysis/enumerate_library_v2.py` (v2, 538, + degradable tails).
4. **Screen** on Modal: `modal_app/lion.py::screen(folds=…, context=…)` at KK/HeLa
   or an organ context; builds the exact LiON screen input via `lion_library.build_library`.
5. **Rank / analyse**: outputs under `data/libraries/combo_v*/` (`*_ranked*.csv`,
   `*_shortlist_top50.csv`, `*_organ_compare.csv`).

**Headline**: the model independently recovered known SAR (unsaturated tails,
3-tail optimum), the 5-fold ensemble gives confidence-aware leads, organ context
shifts the SAR, and degradable disulfide/thioether tails help spleen.

## Track B — active-targeting ligand design (GLP1R pilot)

Flow: **mine (ChEMBL + DrugCLIP + humanPPI) → cross-compare → enrichment → (Boltz validate)**

1. **Mine** by UniProt (APIs reverse-engineered, generic):
   `scripts/fetch_glp1r_ligands.py` (ChEMBL), `scripts/fetch_drugclip.py`
   (drug-the-whole-genome), `scripts/fetch_humanppi.py` (predicted PPIs). See
   [`docs/RESOURCES.md`](docs/RESOURCES.md).
2. **Cross-compare** ChEMBL × DrugCLIP: `analysis/crosscompare_glp1r.py` +
   `analysis/fragment_match_glp1r.py`. Finding: DrugCLIP is a *fragment* library, so
   its hits don't match drug-sized measured actives (no consensus, triangulated).
3. **Protein-layer enrichment**: `analysis/humanppi_enrichment_glp1r.py` — GLP1R's
   interactome is significantly cell-surface enriched; leads DLK1/PAM (β-cell
   membrane proteins). See [`docs/PROBLEM_B_GLP1R.md`](docs/PROBLEM_B_GLP1R.md).
4. **Generic + validated**: `lipidlib/targetpipe.py` + `analysis/run_target.py`
   run all of the above for **any** target in one command (demoed on ASGR1, PSMA);
   `analysis/validate_pipeline.py` proves the enrichment fires (1674× positive
   control). See [`docs/PIPELINE.md`](docs/PIPELINE.md).
5. **Next**: Boltz co-folding (GLP-1 peptide ↔ GLP1R ECD) — hosted MCP, no local GPU.

---

## Repository map

```
lipidlib/            core library
  featurize.py         SMILES -> Morgan/MACCS/RDKit-2D
  lion_library.py      LiON screen-input builder + aza-Michael enumerator
  targetpipe.py        Track-B engine: mine + enrichment maths (any UniProt)
  cli.py               featurizer CLI
modal_app/lion.py    Modal GPU app: train / train_lite / train_resilient / screen
scripts/             data fetchers (LNP_ML, PubChem, ChEMBL, DrugCLIP, humanPPI)
analysis/            enumerations, cross-comparisons, enrichments, run_target, validate
data/                reference_lipids, libraries/combo_v*, targets/<NAME>
results/figures/     all 11 figures
results/reports/     per-target auto reports (run_target)
docs/                PLAN, RESOURCES, PIPELINE, PROBLEM_B_GLP1R, GLP1R_crosscompare,
                     ORGAN_SCREEN, COMBO_LIBRARY, COMBO_V2_SPLEEN
tests/               schema regression tests
external/LNP_ML/     upstream LiON repo (gitignored; fetch via scripts/)
```

## Quickstart

```bash
pip install -r requirements.txt && pip install -e .
python -m pytest tests/ -q

# Track A: enumerate + (on a Modal box) train + screen
python analysis/enumerate_library.py
modal run modal_app/lion.py::train_lite                 # ~6 min, T4
modal run modal_app/lion.py::screen --folds 1 --name combo_v1 \
    --smiles-json "$(python -c 'import pandas,json;print(json.dumps(pandas.read_csv("data/libraries/combo_v1/combo_v1.csv").smiles.tolist()))')"

# Track B: everything for one target, one command
python analysis/run_target.py --uniprot P43220 --name GLP1R
```

## Compute
- **This/any CPU box**: featurization, enumeration, mining, all enrichment analyses, figures.
- **Modal** (token via `MODAL_TOKEN_ID/SECRET`): LiON training + screening (see
  two-tier table in [PLAN](docs/PLAN.md)); Chemprop 1.7 runs in a `uv`-built py3.8
  venv, training is preemption-resilient.
- **Hosted MCP** (live in-session): Boltz (structure + binding), ChEMBL.

## Status & next steps
Track A is a closed loop (enumerate → screen → rank → organ/expand). Track B has
mining + cross-comparison + enrichment + a validated generic pipeline. Open threads:
Boltz structural validation for Track B; batch `run_target` over the receptor table;
larger / more-degradable libraries per organ. Caveats (in-vivo data sparsity,
lite-vs-ensemble, model-conditioning) are noted in each doc — treat predictions as
ranked hypotheses, favour low-ensemble-variance picks.
