#!/usr/bin/env bash
# RiboNN demo: install a CPU env, fetch pretrained weights (Zenodo), and predict
# translation efficiency (TE) for example human transcripts.
#
# Uses `uv` for a fast isolated venv. Falls back to python -m venv + pip if no uv.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${WORK:-$HERE/work}
mkdir -p "$WORK" && cd "$WORK"

# 1. Clone
[ -d RiboNN ] || git clone --depth 1 https://github.com/Sanofi-Public/RiboNN

# 2. CPU env with the versions RiboNN pins (torch 1.13.1, pytorch-lightning 1.8.5)
VENV="$WORK/venv-ribonn"
if command -v uv >/dev/null; then
    uv venv --python 3.10 "$VENV"
    uv pip install --python "$VENV/bin/python" --index-strategy unsafe-best-match \
        --extra-index-url https://download.pytorch.org/whl/cpu \
        "torch==1.13.1+cpu" pytorch-lightning==1.8.5 lightning-utilities==0.10.1 \
        torchmetrics==1.3.1 mlflow==2.18.0 numpy==1.22.4 pandas==2.2.3 \
        scikit-learn==1.0.2 "setuptools<82"
else
    python3.10 -m venv "$VENV"
    "$VENV/bin/pip" install --extra-index-url https://download.pytorch.org/whl/cpu \
        "torch==1.13.1+cpu" pytorch-lightning==1.8.5 lightning-utilities==0.10.1 \
        torchmetrics==1.3.1 mlflow==2.18.0 numpy==1.22.4 pandas==2.2.3 \
        scikit-learn==1.0.2 "setuptools<82"
fi
PY="$VENV/bin/python"

# 3. Predict. The Makefile target auto-downloads weights.zip from Zenodo 17258709
#    (90 human + 90 mouse nested-CV models) if models/ is absent.
cd RiboNN
if [ ! -d models/human ]; then
    mkdir -p tmp
    wget -t 0 -O tmp/weights.zip "https://zenodo.org/records/17258709/files/weights.zip?download=1"
    python3 -c "import zipfile; zipfile.ZipFile('tmp/weights.zip').extractall('models')"
    rm -rf tmp
fi
export MLFLOW_TRACKING_URI=sqlite:///mlruns.db
"$PY" -m src.main --predict human     # use 'mouse' for mouse models

echo "Done. Predictions: $WORK/RiboNN/results/human/prediction_output.txt"
