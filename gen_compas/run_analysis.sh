#!/usr/bin/env bash
# Everything downstream of a finished Gen-COMPAS run.
set -u
cd "$(dirname "$0")"

TAGS="${*:-main rep1}"

for tag in $TAGS; do
  if [ -f "results/gencompas_${tag}_store.npz" ]; then
    echo "=== analysis: $tag ==="
    python3 src/analysis.py --tag "$tag" || echo "  analysis failed for $tag"
  fi
done

MAIN=$(echo "$TAGS" | awk '{print $1}')
if [ -f "results/gencompas_${MAIN}_store.npz" ]; then
  echo "=== committor validation: $MAIN ==="
  python3 src/validate_committor.py --tag "$MAIN" \
    || echo "  validation failed"
  echo "=== generator comparison: $MAIN ==="
  python3 src/compare_generators.py "$MAIN" || echo "  comparison failed"
  echo "=== figures: $MAIN ==="
  python3 src/figures.py "$MAIN" || echo "  figures failed"
fi

echo "=== summary ==="
python3 src/summarize.py | tee results/summary.txt
