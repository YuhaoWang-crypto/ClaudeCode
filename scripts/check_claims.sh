#!/usr/bin/env bash
# Claims integrity (BootLoops emitall): re-emit REPORT.md's headline numbers
# from the certified M23 receipt and fail on any drift. Also checks that every
# quoted string still appears verbatim in REPORT.md.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BL="${BOOTLOOPS_DIR:-$(dirname "$REPO_ROOT")/bootloops}"
cd "$REPO_ROOT"
[ -f receipts/m23_certified_folds.json ] || python3 -m grn_pipeline.m23_certified_folds
python3 "$BL/tools/emitall/run.py" receipts/claims.json
python3 - <<'PY'
import json, sys
report = open("REPORT.md", encoding="utf-8").read().replace("\u2212", "-")
missing = [c["quoted"] for c in json.load(open("receipts/claims.json"))["claims"]
           if c["quoted"] not in report]
if missing:
    sys.exit(f"quotes no longer present in REPORT.md: {missing}")
print("all quoted numbers present verbatim in REPORT.md")
PY
