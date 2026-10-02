#!/usr/bin/env bash
# Install the BootLoops toolkit (github.com/BootLoops-ai/bootloops) beside this
# repository (../bootloops, the sibling layout its INSTALL.md assumes), pinned
# to the commit grn_pipeline was validated against, and run the batteries of
# the packages this project uses.
#
#   bash scripts/setup_bootloops.sh            # install + targeted selftests
#   BOOTLOOPS_DIR=/opt/bootloops bash scripts/setup_bootloops.sh
#   FULL=1 bash scripts/setup_bootloops.sh     # also run all 49 batteries (~10 min)
set -euo pipefail

PIN=66b680ce742e654cfe86da4f072a69061fe182b1
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${BOOTLOOPS_DIR:-$(dirname "$REPO_ROOT")/bootloops}"

if [ ! -d "$DEST/.git" ]; then
  git clone https://github.com/BootLoops-ai/bootloops "$DEST"
fi
git -C "$DEST" fetch --quiet origin "$PIN" 2>/dev/null || git -C "$DEST" fetch --quiet origin
git -C "$DEST" checkout --quiet "$PIN"
echo "BootLoops at $DEST ($(git -C "$DEST" rev-parse --short HEAD))"

python3 -m pip install --quiet mpmath sympy numpy scipy python-flint pytest pyyaml networkx

# packages grn_pipeline drives: baller (ball arithmetic, tripwire, render),
# emitall (claims integrity); clinch/posq are the next candidates (see REPORT M23)
cd "$DEST"
if [ "${FULL:-0}" = "1" ]; then
  python3 run_selftests.py --par 8
else
  python3 run_selftests.py baller emitall
fi
echo
echo "done. grn_pipeline finds it via \$BOOTLOOPS_DIR or ../bootloops."
