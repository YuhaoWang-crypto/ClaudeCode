#!/usr/bin/env bash
# Launch an mdscreen Modal run and stream progress + failure signatures.
#
# Usage:
#   scripts/run_and_watch.sh protein --pdb 1ubq --ns 1
#   scripts/run_and_watch.sh complex --pdb 3poz.pdb --ligand lig.sdf --ns 2
#   scripts/run_and_watch.sh screen  --config examples/screen_config.yaml
#
# Env:
#   MDSCREEN_GPU=A100   # override GPU (default A10G, set in modal_app.py)
set -euo pipefail

ENTRY="${1:?usage: run_and_watch.sh <protein|ligand|complex|screen> [args...]}"
shift
LOG="$(mktemp -t mdscreen.XXXXXX.log)"
echo "==> Launching: modal run modal_app.py::${ENTRY} $*"
echo "==> Log: ${LOG}"

# Run in background; the grep below is the live event stream.
modal run "modal_app.py::${ENTRY}" "$@" > "${LOG}" 2>&1 &
RUN_PID=$!

# Stream the meaningful lines; covers progress AND every failure signature so a
# crash is never silent. Exits when the app completes or the run process dies.
tail -f "${LOG}" --pid="${RUN_PID}" | grep -E --line-buffered \
  "run_id|minimisation|equilibration|production|Analysing|Analysis summary|ΔG|MM/GBSA|SCREEN RANKING|Pull results|Error|Traceback|Exception|UnicodeEncode|Failed|App completed" || true

wait "${RUN_PID}"
echo "==> Run finished. Full log: ${LOG}"
echo "==> Pull results with the 'modal volume get mdscreen-outputs <run_id> ./results' line above."
