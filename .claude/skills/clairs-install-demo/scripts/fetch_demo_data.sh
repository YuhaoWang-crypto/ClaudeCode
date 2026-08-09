#!/usr/bin/env bash
# Download an official ClairS quick-demo dataset (HCC1395 / HCC1395BL tumour-normal
# pair, chr17:80,000,000-80,100,000, with the SEQC2 truth VCF + high-confidence BED).
#
# Usage: fetch_demo_data.sh <ont|ilmn|pacbio_hifi> [DEST_DIR]
set -euo pipefail

PLATFORM="${1:?usage: fetch_demo_data.sh <ont|ilmn|pacbio_hifi> [DEST_DIR]}"
DEST="${2:-$HOME/clairs_demo_data/$PLATFORM}"
# https, not the http:// in upstream docs -- same files, works behind 443-only egress.
BASE="https://www.bio8.cs.hku.hk/clairs/quick_demo/$PLATFORM"

case "$PLATFORM" in
  ont|pacbio_hifi) REF=GRCh38_no_alt_chr17.fa ;;
  ilmn)            REF=GRCh38_chr17.fa ;;
  *) echo "unknown platform '$PLATFORM' (ont|ilmn|pacbio_hifi)" >&2; exit 2 ;;
esac

mkdir -p "$DEST"
FILES=(
  "$REF" "$REF.fai"
  HCC1395BL_normal_chr17_demo.bam HCC1395BL_normal_chr17_demo.bam.bai
  HCC1395_tumor_chr17_demo.bam    HCC1395_tumor_chr17_demo.bam.bai
  SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz
  SEQC2_high-confidence_sSNV_in_HC_regions_v1.2_chr17.vcf.gz.tbi
  SEQC2_High-Confidence_Regions_v1.2_chr17.bed
)
for f in "${FILES[@]}"; do
  if [[ -s "$DEST/$f" ]]; then echo "skip  $f"; continue; fi
  printf 'get   %s ... ' "$f"
  curl -fsSL --retry 4 --retry-delay 2 -o "$DEST/$f" "$BASE/$f" && echo ok
done

printf 'chr17\t80000000\t80100000\n' > "$DEST/quick_demo.bed"
echo
echo "Demo data in $DEST  (reference: $REF)"
du -sh "$DEST"
