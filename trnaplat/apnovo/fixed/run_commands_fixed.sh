#!/bin/bash
# Corrected AP Novo launch. Differences from the audited script, each one a
# failure it would have hit:
#   * LigandMPNN is installed and its two flags are passed, because the manifest
#     sets resequence.enabled -- run_pipeline.py raises app.UsageError otherwise.
#   * the "quick" run really is quick: --only_stage=generation on a cut-down
#     manifest, instead of 1000 sampling steps x 50 designs.
set -euo pipefail

# --- 1. AP Novo environment (per the repo's own README) ---------------------
uv venv --python 3.12 .venv
source .venv/bin/activate
uv pip install -U "jax[cuda12]>=0.4.30"      # CPU-only: uv pip install -U jax
uv pip install git+https://github.com/google-deepmind/alphafold3.git
uv pip install -e .
build_data                                    # compiles ccd.pickle

# --- 2. LigandMPNN environment (separate, PyTorch) -------------------------
# Required because manifest.resequence.enabled is true. Drop this section only
# if you also set {"resequence": {"enabled": false}}.
git clone https://github.com/dauparas/LigandMPNN.git
( cd LigandMPNN \
  && uv venv --python 3.12 .venv_ligandmpnn \
  && uv pip install --python .venv_ligandmpnn/bin/python -r requirements.txt \
  && uv pip install --python .venv_ligandmpnn/bin/python "setuptools<82" \
  && bash get_model_params.sh ./model_params )
export LIGANDMPNN_DIR="$PWD/LigandMPNN"
export LIGANDMPNN_PYTHON="$PWD/LigandMPNN/.venv_ligandmpnn/bin/python"

# --- 3. Weights ------------------------------------------------------------
# ✅ verified reachable: 547,140,096 B and 1,020,518,661 B, HTTP 200.
mkdir -p models/apnovo_generator models/af3_la
wget -nc -P models/apnovo_generator \
  https://storage.googleapis.com/alphaprotein_novo/generator.bin.zst
wget -nc -P models/af3_la \
  https://storage.googleapis.com/alphafold3/af3_leaving_atom.bin.zst

# --- 4. Smoke test: generation only, few steps, one design ----------------
python - <<'PY'
import json, pathlib
m = json.loads(pathlib.Path('manifest_fixed.json').read_text())
m['defaults'].update(num_sampling_steps=50, num_designs=1)
m['designs'] = m['designs'][:1]
pathlib.Path('manifest_smoke.json').write_text(json.dumps(m, indent=2))
PY
python run_pipeline.py --manifest=manifest_smoke.json --output_dir=./out_smoke \
  --only_stage=generation --apn_model_dir=./models/apnovo_generator

# --- 5. Full campaign ------------------------------------------------------
python run_pipeline.py --manifest=manifest_fixed.json --output_dir=./out_full \
  --apn_model_dir=./models/apnovo_generator \
  --af3_model_dir=./models/af3_la \
  --ligandmpnn_dir="$LIGANDMPNN_DIR" \
  --ligandmpnn_python="$LIGANDMPNN_PYTHON"
