#!/usr/bin/env bash
# End-to-end ClairS quick demo: fetch data -> call somatic variants -> benchmark
# against the SEQC2 truth set.  Assumes install_clairs.sh has already run.
#
# Usage: run_demo.sh <ont|ilmn|pacbio_hifi> [WORKDIR] [THREADS]
#   env: CLAIRS_PREFIX (default $HOME/clairs) -- where install_clairs.sh put things
set -euo pipefail

PLATFORM_KEY="${1:?usage: run_demo.sh <ont|ilmn|pacbio_hifi> [WORKDIR] [THREADS]}"
: "${CLAIRS_PREFIX:=$HOME/clairs}"
WORK="${2:-$HOME/clairs_demo/$PLATFORM_KEY}"
THREADS="${3:-4}"
PREFIX="$CLAIRS_PREFIX"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[[ -x "$PREFIX/bin/micromamba" ]] || {
  echo "No ClairS install at $PREFIX -- run install_clairs.sh, or set CLAIRS_PREFIX." >&2
  exit 1
}

# demo dataset -> (reference filename, --platform value for run_clairs)
case "$PLATFORM_KEY" in
  ont)         REF=GRCh38_no_alt_chr17.fa; CLAIRS_PLATFORM=ont_r10_guppy ;;
  pacbio_hifi) REF=GRCh38_no_alt_chr17.fa; CLAIRS_PLATFORM=hifi_revio ;;
  ilmn)        REF=GRCh38_chr17.fa;        CLAIRS_PLATFORM=ilmn ;;
  *) echo "unknown platform '$PLATFORM_KEY'" >&2; exit 2 ;;
esac

IN="$WORK/input"; OUT="$WORK/output"
mkdir -p "$IN" "$OUT"
bash "$HERE/fetch_demo_data.sh" "$PLATFORM_KEY" "$IN"

export MAMBA_ROOT_PREFIX="$PREFIX/mamba"
MM="$PREFIX/bin/micromamba"
RUN=("$MM" run -n clairs)

TRUTH=SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz
HCBED=SEQC2_High-Confidence_Regions_v1.2_chr17.bed
REGION=chr17:80000000-80100000

echo "==> Calling somatic variants (platform=$CLAIRS_PLATFORM, threads=$THREADS)"
"${RUN[@]}" "$PREFIX/ClairS/run_clairs" \
  --tumor_bam_fn  "$IN/HCC1395_tumor_chr17_demo.bam" \
  --normal_bam_fn "$IN/HCC1395BL_normal_chr17_demo.bam" \
  --ref_fn        "$IN/$REF" \
  --threads "$THREADS" \
  --platform "$CLAIRS_PLATFORM" \
  --output_dir "$OUT" \
  --region "$REGION" > "$WORK/run_clairs.stdout" 2>&1
grep -E '^\[INFO\] (STEP|Total time|Finish)' "$WORK/run_clairs.stdout" || true

# run_clairs exits 0 even when a sub-step dies, so check for the artefact itself.
[[ -s "$OUT/output.vcf.gz" ]] || {
  echo "No $OUT/output.vcf.gz -- calling failed; see $WORK/run_clairs.stdout" >&2
  tail -20 "$WORK/run_clairs.stdout" >&2
  exit 1
}

echo "==> Benchmarking against SEQC2 truth"
"${RUN[@]}" python3 "$PREFIX/ClairS/clairs.py" compare_vcf \
  --truth_vcf_fn "$IN/$TRUTH" \
  --input_vcf_fn "$OUT/output.vcf.gz" \
  --bed_fn "$IN/$HCBED" \
  --output_dir "$OUT/benchmark" \
  --input_filter_tag 'PASS' \
  --ctg_name chr17 --ctg_start 80000000 --ctg_end 80100000 \
  2>&1 | tee "$WORK/benchmark.txt"

echo
echo "Somatic VCF : $OUT/output.vcf.gz"
echo "Benchmark   : $WORK/benchmark.txt"
