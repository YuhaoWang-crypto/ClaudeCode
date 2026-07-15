#!/usr/bin/env bash
# End-to-end APR calibration on the dye·β-CD reference system.
#   1. parameterize host + guests
#   2. build APR windows for the dye
#   3. run all windows (GPU)
#   4. analyze + compare to the measured dye Ka
#
# Usage:  scripts/run_calibration.sh [GUEST]     (default GUEST=dye)
set -euo pipefail
cd "$(dirname "$0")/.."

GUEST="${1:-dye}"
export CDPFAS_LOGLEVEL="${CDPFAS_LOGLEVEL:-INFO}"

echo "[1/4] parameterizing host + guests ..."
python -m src.parameterize --config config/system.yaml --out work/params

echo "[2/4] building APR windows for guest=${GUEST} ..."
python -m src.build_apr --config config/system.yaml --params work/params \
    --work work --guest "${GUEST}"

echo "[3/4] running MD for all windows (this needs a GPU) ..."
python -m src.run_apr --config config/system.yaml \
    --work "work/apr_${GUEST}" --window all

echo "[4/4] analyzing + calibrating vs experiment ..."
python -m src.analyze_apr --config config/system.yaml \
    --work "work/apr_${GUEST}" --guest "${GUEST}"

echo "done -> work/apr_${GUEST}/binding_result.json"
