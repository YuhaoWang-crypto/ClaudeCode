#!/usr/bin/env bash
# Official Virtual Cell Challenge 2026 pipeline, end to end.
#
# Prerequisite, run in your own terminal so the token never touches a chat log:
#
#     vcc login --token-stdin
#     vcc datasets download controls -o vcc_2026_controls.zip
#     unzip -o -j vcc_2026_controls.zip -d /home/user/vcc_official
#
# This script does not submit.  Submission is an outward-facing action under
# your account and is left to you; the command is printed at the end.
set -euo pipefail

PY=${PY:-/home/user/venv/bin/python}
VCC=${VCC:-/home/user/venv/bin/vcc}
OFFICIAL=${OFFICIAL:-/home/user/vcc_official}
OUT=${OUT:-/home/user/vcc_submission}
cd "$(dirname "$0")"
mkdir -p "$OUT" results

echo "=== 1. what the bundle contains, and what the atlas can reach in it"
$PY -m virtualcell.vcc2026 --identify --describe 2>&1 | tee results/vcc2026_describe.txt

echo
echo "=== 2. single-source transfer benchmarked on lines whose truth we have"
if [ ! -f results/gwps_single_source.json ]; then
  $PY -m virtualcell.gwps --benchmark 2>&1 | tee results/gwps_bench.log
else
  echo "    results/gwps_single_source.json exists, skipping"
fi
$PY -m virtualcell.report --vcc2026

echo
echo "=== 3. cell-level dry run: the submission pipeline, scored on real cells"
$PY -m virtualcell.vcc_eval --target Jurkat --source gwps --match-panel \
    --n-knockdowns 40 2>&1 | tee results/vcc2026_celllevel.txt

echo
echo "=== 4. build the submission (360,000 cells, ~2.1e9 stored entries)"
$PY -m virtualcell.vcc2026 --build --out "$OUT/prediction.h5ad" \
    2>&1 | tee results/vcc2026_build.txt

echo
echo "=== 5. validate with the official validator"
$VCC prep "$OUT/prediction.h5ad" \
    -g "$OFFICIAL/gene_names.csv" --perts "$OFFICIAL/pert_counts.csv" \
    -o "$OUT/prediction.vcc" -f 2>&1 | tee results/vcc2026_prep.txt

echo
echo "submit it yourself when you are ready:"
echo "    vcc submit $OUT/prediction.vcc -m 'ContextTransfer, genome-wide K562 source'"
