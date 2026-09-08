#!/usr/bin/env bash
# Fetch the OpenKnot AI design benchmark release (eternagame/OpenKnotAIDesignData).
#
# The CSVs are stored with Git LFS, so a plain clone or the Zenodo source archive
# yields pointer stubs, not data. The media endpoint below serves the real files.
#
# Total download ~200 MB into openknot_repro/data/.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${OPENKNOT_DATA:-$ROOT/data}"
BASE="https://media.githubusercontent.com/media/eternagame/OpenKnotAIDesignData/main"

FILES=(
  "Data/OpenKnotBench_data.v4.5.1.csv"      # SHAPE profiles + scores, rounds 1-4
  "Data/OK7a_M2R_data.v4.5.1.csv"           # mutate-map-rescue, round 3 top designs
  "Data/OK7a_M2_data.v4.5.2.csv"            # single-mutant libraries, round 3
  "Targets/Rounds1and2_targets.csv"
  "Targets/Round3_targets.csv"
  "Targets/Round4_targets.csv"
)

mkdir -p "$DEST/Data" "$DEST/Targets"
for file in "${FILES[@]}"; do
  out="$DEST/$file"
  if [[ -s "$out" ]] && ! head -c 64 "$out" | grep -q "git-lfs.github.com"; then
    echo "have    $file"
    continue
  fi
  echo "fetch   $file"
  curl -fsSL -o "$out" "$BASE/$file"
done
echo "data in $DEST"
