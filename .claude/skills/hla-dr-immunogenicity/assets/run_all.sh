#!/usr/bin/env bash
# Full pipeline, in order. Run from the project root (the parent of scripts/).
#
#   bash run_all.sh            # assessment only (M0-M9) + report
#   bash run_all.sh --calibrate  # also re-run the M10-M13 calibration layer
#
# The calibration layer is a property of the panel and thresholds, not of the
# test article: skip it when only the sequence changed.
set -euo pipefail

if [ ! -f data/human_sprot.fasta ]; then
  echo "fetching human proteome (once, ~20400 entries)"
  curl -L 'https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta&compressed=true' \
    | gunzip > data/human_sprot.fasta
fi

run() { echo "=== $1"; python3 "scripts/$1"; }

run m0_fetch_sequences.py
run m1_sequence_qc.py
run m2_panel_design.py
run m3_binding_prediction.py        # ~40 min, resumable
run m4_tolerance_filter.py
run m5_risk_scoring.py
run m6_benchmark_calibration.py     # system suitability: batch is reportable or not
run m7_bcell_layer.py
run m8_exposure_context.py
run m9_deimmunization_scan.py       # optional

if [ "${1:-}" = "--calibrate" ]; then
  run m10_benchmark_fetch.py        # ~10 min
  run m11_threshold_calibration.py  # ~1-2 h, resumable
  run m12_tolerance_weight.py
  run m13_promiscuity_vs_bestrank.py  # ~2 h, resumable
fi

run make_figures.py
run make_report.py
run make_deck.py
run check_deck.py                   # must report no layout issues
echo "report.html and report.pptx written"
