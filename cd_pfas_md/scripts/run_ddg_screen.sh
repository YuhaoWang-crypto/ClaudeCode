#!/usr/bin/env bash
# Batch FEP/TI ΔΔG screen of the CD modification library.
# Requires that parameterize.py has already produced the guest parameters
# (run_calibration.sh does this, or run step 1 below on its own).
#
# Usage:  scripts/run_ddg_screen.sh [GUEST]   (default = screen_guest in modifications.yaml)
set -euo pipefail
cd "$(dirname "$0")/.."
export CDPFAS_LOGLEVEL="${CDPFAS_LOGLEVEL:-INFO}"

GUEST_ARG=()
if [[ $# -ge 1 ]]; then GUEST_ARG=(--guest "$1"); fi

echo "[1/2] ensuring guest parameters exist ..."
python -m src.parameterize --config config/system.yaml --out work/params

echo "[2/2] running ΔΔG screen ..."
python -m src.fep_ti_ddg \
    --config config/system.yaml \
    --modifications config/modifications.yaml \
    --params work/params \
    --work work/ddg \
    "${GUEST_ARG[@]}"

echo "done -> work/ddg/ddg_screen_*.json"
