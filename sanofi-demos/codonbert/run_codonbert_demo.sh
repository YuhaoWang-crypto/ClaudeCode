#!/usr/bin/env bash
# CodonBERT demo: install a CPU env, fetch the pretrained model, and benchmark
# frozen CodonBERT embeddings on the mRFP expression dataset.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${WORK:-$HERE/work}
mkdir -p "$WORK" && cd "$WORK"

# 1. Clone
[ -d CodonBERT ] || git clone --depth 1 https://github.com/Sanofi-Public/CodonBERT

# 2. CPU env. poetry's lockfile is old; we install equivalent pinned versions.
VENV="$WORK/venv-codonbert"
if command -v uv >/dev/null; then
    uv venv --python 3.10 "$VENV"
    uv pip install --python "$VENV/bin/python" \
        --extra-index-url https://download.pytorch.org/whl/cpu \
        "torch==1.13.1+cpu" "transformers==4.29.2" "tokenizers==0.13.3" \
        "numpy<2" "pandas<2" scikit-learn scipy biopython
else
    python3.10 -m venv "$VENV"
    "$VENV/bin/pip" install --extra-index-url https://download.pytorch.org/whl/cpu \
        "torch==1.13.1+cpu" "transformers==4.29.2" "tokenizers==0.13.3" \
        "numpy<2" "pandas<2" scikit-learn scipy biopython
fi
PY="$VENV/bin/python"

# 3. Download + extract the pretrained CodonBERT model
if [ ! -f cb_model/codonbert/pytorch_model.bin ]; then
    mkdir -p cb_model
    curl -sSL -o cb_model/CodonBERT.zip https://cdn.prod.accelerator.sanofi/llm/CodonBERT.zip
    "$PY" -c "import zipfile; zipfile.ZipFile('cb_model/CodonBERT.zip').extractall('cb_model')"
    rm -f cb_model/CodonBERT.zip
fi

# 4. Run the frozen-feature benchmark
export CODONBERT_SRC="$WORK/CodonBERT/benchmarks"
export CODONBERT_MODEL="$WORK/cb_model/codonbert"
"$PY" -I "$HERE/codonbert_mrfp_demo.py"
