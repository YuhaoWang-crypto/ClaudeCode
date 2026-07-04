#!/usr/bin/env bash
# Finish the two pending guide-QC steps (Rule Set 3 on-target + off-target
# specificity) for the kinome library. Runs on a Latch task/Pod or any machine
# with an hg38 index. Inputs are already committed under data/.
#
# Prereqisites (install once):
#   pip install rs3                       # Broad Rule Set 3 on-target model
#   conda install -c bioconda guidescan   # GuideScan2 off-target/specificity
#   # hg38 reference FASTA (e.g. from UCSC/GENCODE), or a prebuilt GuideScan2 index
#
# Usage:
#   bash scripts/run_scoring.sh /path/to/hg38.fa
set -euo pipefail
HG38="${1:?pass path to hg38.fa (or a prebuilt GuideScan2 index dir)}"
cd "$(dirname "$0")/.."
mkdir -p data/scoring

echo "==> [1/3] Rule Set 3 on-target scoring (rs3) on 2,054 contexts"
python3 scripts/score_guides.py ontarget \
    --context data/rs3_context.tsv \
    --out data/scoring/ontarget_rs3.tsv

echo "==> [2/3] Off-target enumeration + specificity (GuideScan2)"
# Build an index once (skip if you already have one):
if [ ! -e hg38_guidescan.index ]; then
  guidescan index --index hg38_guidescan.index "$HG38"
fi
# Enumerate off-targets for our spacers and emit a spacer-keyed specificity table.
# The exact GuideScan2 sub-commands vary by version; the contract for the next
# step is only that the output TSV has columns:  spacer  specificity  n_offtargets
guidescan enumerate \
    --index hg38_guidescan.index \
    --mismatches 3 --pam NGG \
    --kmers data/all_spacers.txt \
    --format tsv \
    --output data/scoring/offtarget_guidescan.tsv
# (Alternative tool: CRISPOR CLI also emits per-guide CFD specificity; or
#  Cas-OFFinder for raw off-target sites if you compute CFD separately.)

echo "==> [3/3] Merge + apply CFD specificity filter (>= 0.2)"
python3 scripts/score_guides.py merge \
    --context data/rs3_context.tsv \
    --ontarget data/scoring/ontarget_rs3.tsv \
    --offtarget data/scoring/offtarget_guidescan.tsv \
    --out data/guides_scored.tsv

echo "Done -> data/guides_scored.tsv"
echo "Next: within each gene keep the top 6 by rs3_score among PASS guides;"
echo "flag any gene left with < 6 for manual review."
