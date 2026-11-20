#!/bin/bash
# Pull the two single-cell files the DE-recall measurement needs, subsample each,
# and delete the raw file only once the subsample is readable.
#
# Target is Jurkat, which the 2026 wrap-up de-anonymised as one of the three
# validation cell lines, so the measurement happens in a context the challenge
# actually scored.  Source is K562, the canonical reference screen -- 277 of the
# 300 test targets are in Replogle's K562 arm.
#
# The target keeps few perturbations with ALL their cells, because the metric's
# four DE members are a per-group statistical test and a thin group is worth
# less than no group.  Nadig's Jurkat holds 262,956 cells over ~2,000 targets,
# so 400 cells each was never on offer; --max-perts takes the largest instead.
# The source only ever becomes a pseudobulk mean, so it keeps many
# perturbations thinly.
#
# Two failure modes this guards against, both of which bit the first attempt:
# the subsample's exit status was being read through a pipe (so `tail`'s zero
# hid an out-of-memory kill), and the raw file was deleted on the failure path
# (so a retry had to re-download 9.4 GB).
set -u
cd /home/user/ClaudeCode
export OMP_NUM_THREADS=2
DEST=/home/user/vcc_data/singlecell
RAW=/home/user/vcc_data/raw
mkdir -p "$DEST" "$RAW"

GEO=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl
FIG=https://ndownloader.figshare.com/files

readable() {   # the only check that means anything: can anndata open it?
  python -c "
import anndata, sys
a = anndata.read_h5ad(sys.argv[1])
assert a.n_obs > 0 and a.n_vars > 0
print(f'  verified {a.n_obs:,} cells x {a.n_vars:,} genes')
" "$1"
}

fetch() {   # $1 url  $2 raw name  $3 out base  $4 per-pert  $5 per-control  $6 max-perts
  local url="$1" raw="$RAW/$2" out="$DEST/$3.h5ad"
  if [ -s "$out" ] && readable "$out" >/dev/null 2>&1; then
    echo "== have $3"; readable "$out"; return 0
  fi
  if [ ! -s "$raw" ]; then
    echo "== downloading $3"
    curl -sSL --retry 4 --retry-delay 5 -o "$raw" "$url" \
      || { echo "  DOWNLOAD FAILED $3"; return 1; }
  else
    echo "== reusing raw $2"
  fi
  ls -la "$raw"
  echo "== subsampling $3: per-pert $4, control $5, max-perts ${6:-all}"
  local extra=()
  [ -n "${6:-}" ] && extra=(--max-perts "$6")
  if ! python -m virtualcell.prep_singlecell "$raw" "$out" \
        --per-pert "$4" --per-control "$5" "${extra[@]}"; then
    echo "  SUBSAMPLE FAILED $3 -- keeping $raw for a retry"; return 1
  fi
  if ! readable "$out"; then
    echo "  OUTPUT UNREADABLE $3 -- keeping $raw for a retry"; rm -f "$out"; return 1
  fi
  rm -f "$raw"
  df -h / | tail -1
}

fetch "$GEO/GSE264667_jurkat_raw_singlecell_01.h5ad" jurkat.h5ad jurkat_de 400 8000 150 \
  || { echo "ABORT: target missing"; exit 1; }
fetch "$FIG/35773219" k562.h5ad k562_de 60 4000 "" \
  || { echo "ABORT: source missing"; exit 1; }

echo "== done"; ls -la "$DEST"; df -h / | tail -1
