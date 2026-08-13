#!/bin/bash
# Re-fetch the raw Perturb-seq single-cell files and subsample each to a size
# that fits alongside everything else, deleting the raw file before the next
# download starts.  Peak disk ~11 GB instead of the 34 GB the four raws total.
#
# Stack consumes cells; this project's own pipeline is pseudobulk, so these are
# needed only for the Stack comparison.
set -u
cd /home/user/ClaudeCode
export OMP_NUM_THREADS=2
PY=/home/user/venv/bin/python
DEST=/home/user/vcc_data/singlecell
mkdir -p "$DEST" /home/user/vcc_data/raw

GEO=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl
FIG=https://ndownloader.figshare.com/files

fetch() {   # $1 url  $2 raw filename  $3 output npz/h5ad basename
  local url="$1" raw="/home/user/vcc_data/raw/$2" out="$DEST/$3.h5ad"
  if [ -s "$out" ]; then echo "== have $3"; return; fi
  echo "== downloading $3"
  curl -sSL --retry 4 --retry-delay 5 -o "$raw" "$url" || { echo "  DOWNLOAD FAILED $3"; return; }
  ls -la "$raw"
  echo "== subsampling $3"
  $PY -m virtualcell.prep_singlecell "$raw" "$out" --per-pert 20 --per-control 2000 \
    2>&1 | grep -vE "^    [0-9]" || echo "  SUBSAMPLE FAILED $3"
  rm -f "$raw"
  df -h / | tail -1
}

fetch "$GEO/GSE264667_hepg2_raw_singlecell_01.h5ad"  hepg2.h5ad   hepg2
fetch "$GEO/GSE264667_jurkat_raw_singlecell_01.h5ad" jurkat.h5ad  jurkat
fetch "$FIG/35775606"                                rpe1.h5ad    rpe1
fetch "$FIG/35773219"                                k562.h5ad    k562

echo "== done"; ls -la "$DEST"; df -h / | tail -1
