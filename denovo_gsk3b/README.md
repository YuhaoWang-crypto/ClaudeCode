# GSK3B de novo design campaign

A graph-genetic-algorithm pipeline that generates novel small molecules against
GSK3-beta (GSK3B, UniProt P49841, ChEMBL CHEMBL262), scores them on a
multi-objective function, filters and clusters them into a diverse shortlist,
proposes retrosynthetic routes, and builds a PDF report.

**Everything it produces is computational. No compound here has been synthesised
or assayed.** Read the Limitations section of the generated report before acting
on any ranking — in particular, the activity score is a surrogate model that the
generator optimised against directly.

## Results

The campaign output lives in `results/GSK3B_generative_design/`:

| File | What it is |
|---|---|
| `GSK3B_de_novo_design_report.pdf` | The report (8 pages, figures + references) |
| `top10_designs.csv` / `.json` | Selected designs with properties, novelty, alerts |
| `cascade_survivors.json` | All 248 molecules passing the full filter cascade |
| `ga_all_scored.json` | Every unique molecule generated and scored (8,175) |
| `ga_history.json` | Per-generation convergence statistics |
| `cascade_log.json` | Attrition at each filter stage |
| `retro_results.json` | AiZynthFinder route search results |
| `figures/` | Figures used in the report |

## Modules

| Module | Role |
|---|---|
| `scoring.py` | TDC oracles (GSK3B / QED / SA) + the conjunctive objective and property gate |
| `graph_ga.py` | Crossover and mutation operators on molecular graphs |
| `run_campaign.py` | The GA loop, seeding, and run artefacts |
| `filters.py` | Novelty, PAINS, ring sanity, property calculation, filter cascade |
| `select.py` | Butina clustering and chemotype-aware representative selection |
| `run_retro.py` | AiZynthFinder retrosynthesis (separate venv — see below) |
| `figures.py` | Report figures |
| `report.py` | PDF assembly |

## Reproducing

Two virtualenvs are needed because AiZynthFinder's pins conflict with the
scikit-learn version required to unpickle the TDC oracles.

```bash
# main env
uv venv .venv --python 3.11 && . .venv/bin/activate
uv pip install rdkit PyTDC "setuptools<81" "scikit-learn==1.2.2" \
               matplotlib reportlab pypdf pillow
```

`scikit-learn==1.2.2` is **required**: the TDC GSK3B oracle ships a forest
pickled before `missing_go_to_left` was added to sklearn's tree node dtype
(sklearn 1.3), so newer versions raise `ValueError` on load.

```bash
# 1. fetch reference actives from ChEMBL (pChEMBL >= 7 against CHEMBL262)
#    -> gsk3b_actives.json
# 2. run the GA
python -m denovo_gsk3b.run_campaign \
    --actives gsk3b_actives.json \
    --out results/GSK3B_generative_design \
    --pop-size 200 --generations 40 --seed 42
# 3. cascade + selection + figures + report (see the campaign driver)
```

Retrosynthesis, in its own env:

```bash
uv venv .retro --python 3.11 && . .retro/bin/activate
uv pip install "aizynthfinder[all]==4.4.1"
mkdir -p aizynth_models && cd aizynth_models && download_public_data .  # ~770 MB
python denovo_gsk3b/run_retro.py \
    --designs results/GSK3B_generative_design/top10_designs.json \
    --config aizynth_models/config.yml \
    --out results/GSK3B_generative_design/retro_results.json
```

The GA is seeded (`--seed 42`) and runs on CPU in about 4 minutes for
200 x 40. The TDC oracles download ~36 MB of pickles to `./oracle/` on first
use; that directory is gitignored.

## Design notes

**The objective is conjunctive on purpose.** Score is the geometric mean of
activity, QED and a normalised SA term, so a zero in any term zeroes the score —
the GA cannot buy predicted potency with unmakeability. A hard property gate
(MW, heavy atoms, cLogP, charge, fragment count, ring size, elements) backs this
up, because a GA optimising a surrogate otherwise drifts into high-scoring
nonsense.

**Selection is chemotype-aware, and this matters.** Taking the ten
highest-scoring survivors yields ten decorations of a single core. Bemis-Murcko
scaffold uniqueness does not fix it — swapping a pendant ring changes the Murcko
scaffold while leaving the recognition motif intact. Selection therefore clusters
survivors (Butina, Tanimoto 0.4) and takes the best representative per cluster,
which costs ~0.03 mean score and buys an actual portfolio.

**Structural alerts are annotated, not filtered.** Several validated GSK3B
chemotypes are maleimides, so filtering Michael acceptors would discard real
chemistry; instead Brenk alerts and Michael-acceptor SMARTS are recorded per
design and surfaced in the report. In this run 6 of 10 designs tripped an alert
and 5 were maleimides — a finding, not a footnote.
