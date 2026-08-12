#!/usr/bin/env bash
# Rebuild every input this project uses, from primary sources.
#
# Downloads ~16 GB but never holds more than ~10 GB at once: each GEO
# single-cell file is pseudobulked to ~37 MB and deleted before the next one
# starts.  Budget ~30 GB of free disk and 30-60 minutes on a decent connection.
#
#   ./fetch_data.sh [destination]        # default /home/user/vcc_data
set -euo pipefail

DEST="${1:-/home/user/vcc_data}"
PY="${PYTHON:-python3}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

mkdir -p "$DEST/replogle" "$DEST/nadig" "$DEST/pseudobulk"

echo "==> Replogle et al. Cell 2022 pseudobulk (Figshare+ 20029387)"
# Pseudobulk, not single-cell: the K562 genome-wide single-cell matrix alone is
# 66 GB and nothing here needs it.
fetch_figshare() {   # $1 = ndownloader file id, $2 = output name
  local out="$DEST/replogle/$2"
  [ -s "$out" ] && { echo "    have $2"; return; }
  echo "    $2"
  curl -sSL --retry 4 --retry-delay 3 -o "$out" \
    "https://ndownloader.figshare.com/files/$1"
}
fetch_figshare 35773070 K562_essential_raw_bulk.h5ad     #  80 MB
fetch_figshare 35775581 rpe1_raw_bulk.h5ad               #  95 MB
fetch_figshare 35774443 K562_gwps_raw_bulk.h5ad          # 375 MB, optional
fetch_figshare 35780870 K562_essential_norm_bulk.h5ad    #  80 MB, optional
fetch_figshare 35775512 rpe1_norm_bulk.h5ad              #  95 MB, optional

echo "==> Nadig et al. Nat Genet 2025 single-cell (GEO GSE264667)"
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl
for pair in "hepg2:5.6GB" "jurkat:9.4GB"; do
  name="${pair%%:*}"; size="${pair##*:}"
  target="$DEST/pseudobulk/$name.npz"
  [ -s "$target" ] && { echo "    have $name.npz"; continue; }

  raw="$DEST/nadig/GSE264667_${name}_raw_singlecell_01.h5ad"
  echo "    downloading $name ($size) ..."
  curl -sSL --retry 4 --retry-delay 5 -o "$raw" \
    "$GEO/GSE264667_${name}_raw_singlecell_01.h5ad"

  echo "    pseudobulking $name ..."
  ( cd "$REPO" && "$PY" -m virtualcell.prep_nadig "$raw" "$target" )

  # Free the raw file immediately; both raws together will not fit alongside
  # the rest on a typical session allowance.
  rm -f "$raw"
done

echo
echo "==> done"
du -sh "$DEST"/* 2>/dev/null || true
echo
echo "Verify with:  cd $REPO && $PY -m virtualcell.data"
echo "Expected:     4 cell lines, 6642 shared genes, 2053 shared perturbations"
