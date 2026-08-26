#!/usr/bin/env bash
# Verify this skill's miner reproduces the published pipeline.
#
# Fetches the 200-structure demo set from the paper's reference repository
# (10 experimentally characterized radical halogenases + 190 negatives), runs
# `mcmine benchmark`, and asserts the published numbers:
#
#     203 2-His sites | 192 facial-triad sites | 10 hits | recall 1.0 | specificity 1.0
#
# Usage:  bash validate_against_paper.sh [workdir]
# Needs:  git, python3 with biopython + pandas + numpy (~50 MB download).

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${1:-${TMPDIR:-/tmp}/mcmine-validation}"
REPO="https://github.com/yannikipouros/hal-discovery.git"
PY="${PYTHON:-python3}"

mkdir -p "$WORK"
if [ ! -d "$WORK/hal-discovery" ]; then
  echo "Cloning reference repository (demo dataset ~47 MB)…"
  git clone --depth 1 "$REPO" "$WORK/hal-discovery"
fi

DEMO="$WORK/hal-discovery/demo/demo_dataset"
[ -d "$DEMO" ] || { echo "Demo dataset not found at $DEMO" >&2; exit 1; }

"$PY" "$HERE/mcmine.py" benchmark \
  --struct-dir "$DEMO" \
  --motif "$HERE/motifs/fe_akg_radical_halogenase.json" \
  --positives "$HERE/benchmarks/radical_halogenase_positives.txt" \
  --outdir "$WORK/out" > "$WORK/benchmark.log"

"$PY" "$HERE/mcmine.py" mine \
  --struct-dir "$DEMO" \
  --motif "$HERE/motifs/fe_akg_radical_halogenase.json" \
  --outdir "$WORK/out" > "$WORK/mine.log"

"$PY" - "$WORK/out" <<'PY'
import glob, json, sys

outdir = sys.argv[1]
summaries = [json.load(open(p)) for p in sorted(glob.glob(f"{outdir}/summary_*.json"))]
mine = next(s for s in summaries if "hits" in s and "total_sites" in s)
bench = next(s for s in summaries if "baseline" in s)

expected = {
    "total_sites": 203,
    "sites_2His_1AspGlu": 192,
    "hits": 10,
    "recall": 1.0,
    "specificity": 1.0,
}
actual = {
    "total_sites": mine["total_sites"],
    "sites_2His_1AspGlu": mine["contrast_sites"],
    "hits": mine["hits"],
    "recall": bench["baseline"]["recall"],
    "specificity": bench["baseline"]["specificity"],
}

ok = True
for key, want in expected.items():
    got = actual[key]
    flag = "OK  " if got == want else "FAIL"
    ok &= got == want
    print(f"{flag} {key}: expected {want}, got {got}")

if bench["baseline"]["missed"]:
    ok = False
    print("FAIL missed positives:", bench["baseline"]["missed"])
if bench["baseline"]["false_positives"]:
    ok = False
    print("FAIL false positives:", bench["baseline"]["false_positives"])

print("\nVALIDATION PASSED" if ok else "\nVALIDATION FAILED")
sys.exit(0 if ok else 1)
PY
