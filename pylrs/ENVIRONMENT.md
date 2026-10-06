# Running wendao/PylRS-libY in 2026

The published model is an AutoGluon **0.4.0** pickle from 2022-06-29 that needs
**Python 3.8**. The upstream Colab notebook gets there via `condacolab`, which is
slow and fragile. `uv` does it in about two minutes and does not touch the system
Python.

✅ Verified working on linux-x86_64, 2026-10-06: the pipeline reproduces the
published Colab output **digit for digit** (see `run_prediction.py`).

## Recipe

```bash
# 1. the repo (~800 MB: the AutoGluon model files dominate)
git clone https://github.com/wendao/PylRS-libY.git

# 2. a Python 3.8 environment -- uv downloads the interpreter itself
uv venv --python 3.8 venv38
uv pip install --python venv38/bin/python "autogluon.tabular[all]==0.4.0"

# 3. reproduce
cd PylRS-libY/4.prediction
../../venv38/bin/python -W ignore -c "
from autogluon.tabular import TabularPredictor
p = TabularPredictor.load('AutogluonModels/ag-20220629_024330/')
print(p.problem_type, p.eval_metric, p.get_model_best())"
```

`autogluon.tabular[all]` is enough — the full `autogluon` meta-package also pulls
the text/vision stacks, which this model does not use. The upstream notebook's
`spacy==3.2.3` / `tokenizers==0.10.1` pins exist only to make the full
meta-package resolve, and are unnecessary here.

## Expected noise

The load emits a wall of
`UserWarning: Trying to unpickle estimator ... from version 0.24.2 when using
version 1.0.2`. These are benign here — the bit-exact reproduction is the proof.
`-W ignore` silences them.

## What needs more than this environment

The inference path runs on the shipped pre-computed features. **Scoring a new
ncAA does not**, because the 18 substrate-aware features come from Rosetta:

| Step | Tool | Obtainable? |
|---|---|---|
| ncAA conformers | RDKit | ✅ free |
| ncAA-AMP ligand params, docking, `cartesian_ddg` | **Rosetta** | ⚠️ free for academics, licence required; not installable unattended |
| RESP charges (optional, improves params) | ORCA + Multiwfn | ⚠️ free registration |
| variant structures | ColabFold / AlphaFold2 | ✅ free |
| ESM-1v ×5, ESM-msa-1b | fair-esm | ✅ free, needs GPU for comfort |
| 256-d eUniRep | jax-unirep, **fine-tuned on a PylRS MSA** | ✅ free; the shipped weights are PylRS-specific |

So the honest status is: **inference reproduces today; a new substrate needs a
Rosetta licence and a cluster.** Budget that before planning around this model.
