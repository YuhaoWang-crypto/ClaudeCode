#!/usr/bin/env bash
# Download the Sade-Feldman et al. (Cell 2018) melanoma scRNA-seq data
# from GEO series GSE120575. Files are large and are NOT committed
# (see ../.gitignore); run this once before analysis/run_scrna_analysis.py.
set -euo pipefail

DEST="$(cd "$(dirname "$0")/.." && pwd)/data"
mkdir -p "$DEST"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl"

echo "Downloading per-cell clinical metadata (small)..."
curl -fsSL -o "$DEST/GSE120575_patient_ID_single_cells.txt.gz" \
  "$BASE/GSE120575_patient_ID_single_cells.txt.gz"

echo "Downloading single-cell TPM matrix (~120 MB gz)..."
curl -fsSL -o "$DEST/GSE120575_TPM.txt.gz" \
  "$BASE/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"

echo "Done. Files in $DEST:"
ls -lh "$DEST"
