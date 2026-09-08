#!/usr/bin/env bash
# Fetch gRNAde's source and the checkpoint used for the OpenKnot benchmark.
#
# The repository carries the paper's own benchmark project under
# projects/openknot_benchmark: the Round 3 and Round 4 target metadata and the
# 3D structures that were handed to designers.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THIRD_PARTY="$ROOT/third_party"
WEIGHTS="$ROOT/weights/grnade"

GRNADE_COMMIT="2453b18778a1981ac25f314790e5395f1e9ab218"
DEST="$THIRD_PARTY/geometric-rna-design"

mkdir -p "$THIRD_PARTY" "$WEIGHTS"
if [[ -d "$DEST/.git" ]]; then
  echo "have    geometric-rna-design"
else
  echo "clone   geometric-rna-design"
  GIT_LFS_SKIP_SMUDGE=1 git clone -q https://github.com/chaitjo/geometric-rna-design "$DEST"
  git -C "$DEST" checkout -q "$GRNADE_COMMIT"
fi

fetch_hf() {
  local file="$1" out="$2"
  if [[ -s "$out" ]]; then
    echo "have    $(basename "$out")"
    return
  fi
  echo "fetch   $(basename "$out")"
  local auth=()
  [[ -n "${HF_TOKEN:-}" ]] && auth=(-H "Authorization: Bearer $HF_TOKEN")
  curl -fsSL "${auth[@]}" -o "$out" "https://huggingface.co/chaitjo/gRNAde/resolve/main/$file"
}

fetch_hf "gRNAde_drop3d@0.75_maxlen@500.h5" "$WEIGHTS/gRNAde_drop3d@0.75_maxlen@500.h5"
echo "gRNAde ready in $DEST"
