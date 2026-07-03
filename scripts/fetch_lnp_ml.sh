#!/usr/bin/env bash
# Clone the LiON / LNP_ML repo (Witten et al. 2024, MIT-licensed) next to this
# project so we can reuse data/all_data.csv and its Chemprop D-MPNN pipeline.
#
# NOTE: LNP_ML pins Chemprop 1.7.0 / Python 3.8. Install it in its OWN venv,
# not the repo's main env, to avoid dependency clashes with modern RDKit.
set -euo pipefail

DEST="${1:-external/LNP_ML}"
mkdir -p "$(dirname "$DEST")"

if [ -d "$DEST/.git" ]; then
  echo "already cloned at $DEST; pulling latest"
  git -C "$DEST" pull --ff-only
else
  git clone --depth 1 https://github.com/jswitten/LNP_ML.git "$DEST"
fi

echo
echo "Key files:"
echo "  $DEST/data/all_data.csv         # merged >9k-datapoint training set"
echo "  $DEST/main_script.py            # split / train / predict / screen"
echo
echo "Next (in a separate py3.8 venv):"
echo "  pip install chemprop==1.7.0"
echo "  python main_script.py --help"
