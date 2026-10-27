#!/usr/bin/env bash
# =============================================================================
# 03_run_pipeline.sh
#   编辑前(WT_parent) vs 编辑后(EDIT_pool) 细菌全基因组重测序比对流程
#
#   QC -> 比对 -> 去重 -> 变异检测(SNV/InDel) -> 过滤 -> 两样本比较
#      -> 结构变异(SV) -> 功能注释 -> 覆盖度统计
#
#   每一步都把"数据此刻长什么样"的快照写到 work/snapshots/, 供报告引用。
# =============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$PWD"
REF="$ROOT/ref/DEMO_chr1.fasta"
THREADS=4
SAMPLES=(WT_parent EDIT_pool)

mkdir -p work/{qc,align,vcf,sv,cov,snapshots}
SNAP="$ROOT/work/snapshots"

log() { printf '\n\033[1;36m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*"; }

# -----------------------------------------------------------------------------
log "STEP 0 | 原始数据快照"
# -----------------------------------------------------------------------------
for S in "${SAMPLES[@]}"; do
  { zcat "rawdata/${S}_R1.fastq.gz" || true; } | head -8 > "$SNAP/00_raw_${S}_R1.txt"
  { zcat "rawdata/${S}_R2.fastq.gz" || true; } | head -8 > "$SNAP/00_raw_${S}_R2.txt"
  echo "$S  reads: $(( $(zcat "rawdata/${S}_R1.fastq.gz" | wc -l) / 4 )) pairs"
done

# -----------------------------------------------------------------------------
log "STEP 1 | fastp 质控与接头切除"
# -----------------------------------------------------------------------------
for S in "${SAMPLES[@]}"; do
  fastp \
    -i "rawdata/${S}_R1.fastq.gz" -I "rawdata/${S}_R2.fastq.gz" \
    -o "work/qc/${S}_R1.clean.fastq.gz" -O "work/qc/${S}_R2.clean.fastq.gz" \
    --detect_adapter_for_pe \
    --trim_poly_g --poly_g_min_len 10 \
    --cut_tail --cut_tail_mean_quality 20 \
    --qualified_quality_phred 15 --unqualified_percent_limit 40 \
    --length_required 50 --thread $THREADS \
    -j "work/qc/${S}.fastp.json" -h "work/qc/${S}.fastp.html" \
    2> "work/qc/${S}.fastp.log"
  echo "  $S  -> work/qc/${S}.fastp.json"
done
cp work/qc/EDIT_pool.fastp.json "$SNAP/01_fastp_EDIT_pool.json"

# -----------------------------------------------------------------------------
log "STEP 2 | bwa-mem 比对 -> 排序 -> 标记重复"
# -----------------------------------------------------------------------------
[ -f "${REF}.bwt" ] || bwa index "$REF" 2>/dev/null
[ -f "${REF}.fai" ] || samtools faidx "$REF"

for S in "${SAMPLES[@]}"; do
  RG="@RG\tID:${S}\tSM:${S}\tPL:ILLUMINA\tLB:lib_${S}\tPU:HXXKMDRXY"
  bwa mem -t $THREADS -R "$RG" "$REF" \
      "work/qc/${S}_R1.clean.fastq.gz" "work/qc/${S}_R2.clean.fastq.gz" \
      2> "work/align/${S}.bwa.log" \
    | samtools fixmate -m -@ 2 - - \
    | samtools sort -@ 2 -o "work/align/${S}.sorted.bam" -
  samtools markdup -@ 2 -s "work/align/${S}.sorted.bam" "work/align/${S}.bam" \
      2> "work/align/${S}.markdup.log"
  samtools index "work/align/${S}.bam"
  rm -f "work/align/${S}.sorted.bam"

  samtools flagstat "work/align/${S}.bam" > "work/align/${S}.flagstat.txt"
  samtools stats    "work/align/${S}.bam" > "work/align/${S}.stats.txt"
  echo "  $S  $(grep -m1 'mapped (' "work/align/${S}.flagstat.txt")"
done

# 比对后的数据快照: SAM 记录本身 + 编辑位点附近带 InDel 的 CIGAR
{ samtools view -h work/align/EDIT_pool.bam || true; } | grep -m3 "^@" > "$SNAP/02_bam_header.txt"
{ samtools view work/align/EDIT_pool.bam || true; } | head -3 | cut -c1-400 > "$SNAP/02_bam_records.txt"
samtools view work/align/EDIT_pool.bam DEMO_chr1:12600-12640 \
  | awk '$6 ~ /I/' | head -3 | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' \
  > "$SNAP/02_bam_cigar_insertion.txt" || true
samtools view work/align/EDIT_pool.bam DEMO_chr1:25880-25920 \
  | awk '$6 ~ /S/ || $2 >= 2048' | head -3 | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' \
  > "$SNAP/02_bam_split_reads.txt" || true

# -----------------------------------------------------------------------------
log "STEP 3 | 碱基堆叠 (pileup) 快照"
# -----------------------------------------------------------------------------
for RGN in "DEMO_chr1:12612-12620:on-target_pheA" \
           "DEMO_chr1:26310-26320:on-target_lacZ" \
           "DEMO_chr1:38224-38232:off-target_yjgF"; do
  R="${RGN%:*}"; TAG="${RGN##*:}"
  {
    echo "# region $R  ($TAG)"
    echo "# columns: CHROM POS REF  [WT: DP BASES QUALS]  [EDIT: DP BASES QUALS]"
    samtools mpileup -f "$REF" -r "$R" -q 20 -Q 20 \
        work/align/WT_parent.bam work/align/EDIT_pool.bam 2>/dev/null \
      | cut -f1-9
  } > "$SNAP/03_pileup_${TAG}.txt"
done

# -----------------------------------------------------------------------------
log "STEP 4 | bcftools 联合检测 SNV / 小 InDel"
# -----------------------------------------------------------------------------
# 细菌是单倍体, 但"编辑后"是一个混合群体(存在亚克隆), 所以这里用二倍体模型
# 做联合检测以保留等位基因频率信息, 再用 AD 自行计算 VAF。
bcftools mpileup \
    -f "$REF" -a AD,ADF,ADR,DP,SP -q 20 -Q 20 -d 1000 -L 1000 -Ou \
    work/align/WT_parent.bam work/align/EDIT_pool.bam 2> work/vcf/mpileup.log \
  | bcftools call -m --ploidy 2 -v -Ou 2> work/vcf/call.log \
  | bcftools norm -f "$REF" -m -any -Oz -o work/vcf/joint.raw.vcf.gz 2> work/vcf/norm.log
bcftools index -f work/vcf/joint.raw.vcf.gz
echo "  raw calls: $(bcftools view -H work/vcf/joint.raw.vcf.gz | wc -l)"
bcftools view work/vcf/joint.raw.vcf.gz | grep -v "^##" | cut -c1-240 > "$SNAP/04_vcf_raw.txt"

# 单倍体模型对照 (说明为什么亚克隆事件会被漏掉)
for S in "${SAMPLES[@]}"; do
  bcftools mpileup -f "$REF" -a AD,DP -q 20 -Q 20 -d 1000 -L 1000 -Ou \
      "work/align/${S}.bam" 2>/dev/null \
    | bcftools call -m --ploidy 1 -v -Ou 2>/dev/null \
    | bcftools norm -f "$REF" -m -any -Oz -o "work/vcf/${S}.haploid.vcf.gz" 2>/dev/null
  bcftools index -f "work/vcf/${S}.haploid.vcf.gz"
done
echo "  haploid-model calls  WT=$(bcftools view -H work/vcf/WT_parent.haploid.vcf.gz|wc -l)  EDIT=$(bcftools view -H work/vcf/EDIT_pool.haploid.vcf.gz|wc -l)"

# -----------------------------------------------------------------------------
log "STEP 5 | 变异过滤 (QUAL / DP / MQ / 链偏好)"
# -----------------------------------------------------------------------------
bcftools filter -Ou \
    -s LowQual  -e 'QUAL<30'                   work/vcf/joint.raw.vcf.gz \
  | bcftools filter -Ou -m + -s LowDP  -e 'INFO/DP<20' \
  | bcftools filter -Ou -m + -s LowMQ  -e 'INFO/MQ<40' \
  | bcftools filter -Ou -m + -s LowAltDepth -e 'MAX(FORMAT/AD[*:1]) < 5' \
  | bcftools filter -Oz -m + -s StrandBias -e 'INFO/DP4[0]+INFO/DP4[1] > 0 && INFO/DP4[2]+INFO/DP4[3] > 0 && (INFO/DP4[2]/(INFO/DP4[2]+INFO/DP4[3]+0.001) > 0.95 || INFO/DP4[3]/(INFO/DP4[2]+INFO/DP4[3]+0.001) > 0.95)' \
    -o work/vcf/joint.filtered.vcf.gz
bcftools index -f work/vcf/joint.filtered.vcf.gz

# 只保留 PASS 的位点进入下游分析; 被过滤掉的记录单独留档
bcftools view -f PASS -Oz -o work/vcf/joint.pass.vcf.gz work/vcf/joint.filtered.vcf.gz
bcftools index -f work/vcf/joint.pass.vcf.gz
bcftools view -H -i 'FILTER!="PASS"' work/vcf/joint.filtered.vcf.gz \
  | cut -c1-200 > work/vcf/filtered_out.txt
echo "  PASS: $(bcftools view -H work/vcf/joint.pass.vcf.gz | wc -l) / $(bcftools view -H work/vcf/joint.filtered.vcf.gz | wc -l)  (剔除 $(wc -l < work/vcf/filtered_out.txt) 条)"
bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/MQ\t%FILTER\n' \
    work/vcf/joint.filtered.vcf.gz > "$SNAP/05_vcf_filtered.txt"

# -----------------------------------------------------------------------------
log "STEP 6 | 两样本比较: 哪些变异是编辑带来的?"
# -----------------------------------------------------------------------------
# 6a. 单样本分别 call 后取交集/差集 (bcftools isec)
for S in "${SAMPLES[@]}"; do
  bcftools mpileup -f "$REF" -a AD,DP -q 20 -Q 20 -d 1000 -L 1000 -Ou \
      "work/align/${S}.bam" 2>/dev/null \
    | bcftools call -m --ploidy 2 -v -Ou 2>/dev/null \
    | bcftools norm -f "$REF" -m -any -Oz -o "work/vcf/${S}.vcf.gz" 2>/dev/null
  bcftools index -f "work/vcf/${S}.vcf.gz"
done
rm -rf work/vcf/isec
bcftools isec -p work/vcf/isec -Oz work/vcf/WT_parent.vcf.gz work/vcf/EDIT_pool.vcf.gz 2>/dev/null
echo "  WT-only   : $(bcftools view -H work/vcf/isec/0000.vcf.gz | wc -l)"
echo "  EDIT-only : $(bcftools view -H work/vcf/isec/0001.vcf.gz | wc -l)"
echo "  shared    : $(bcftools view -H work/vcf/isec/0002.vcf.gz | wc -l)"
{
  echo "# bcftools isec  (0000 = WT 独有, 0001 = EDIT 独有, 0002 = 两者共有)"
  for F in 0000 0001 0002; do
    echo "## $F"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "work/vcf/isec/${F}.vcf.gz"
  done
} > "$SNAP/06_isec.txt"

# -----------------------------------------------------------------------------
log "STEP 7 | delly 结构变异检测 (大片段缺失)"
# -----------------------------------------------------------------------------
delly call -g "$REF" -o work/sv/joint.bcf \
    work/align/WT_parent.bam work/align/EDIT_pool.bam > work/sv/delly.log 2>&1
bcftools index -f work/sv/joint.bcf
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%FILTER\t%INFO/PE\t%INFO/SR[\t%SAMPLE=%GT:RV=%RV:DV=%DV:RR=%RR:DR=%DR]\n' \
    work/sv/joint.bcf > "$SNAP/07_delly_sv.txt"
echo "  SV records: $(bcftools view -H work/sv/joint.bcf | wc -l)"
cat "$SNAP/07_delly_sv.txt"

# -----------------------------------------------------------------------------
log "STEP 8 | snpEff 功能注释"
# -----------------------------------------------------------------------------
# snpEff 需要一份带密码子表的配置; 从系统默认配置派生, 再追加本 demo 的基因组条目
mkdir -p snpeff/data/DEMO
cp -f ref/DEMO_chr1.gff3 snpeff/data/DEMO/genes.gff
# cds.fa / protein.fa 让 snpEff build 能做完整校验, 而不是跳过
python3 scripts/02b_make_cds_protein.py > work/vcf/cds_protein.log
if [ ! -f snpeff/snpEff.config ]; then
  SYS_CFG=""
  for c in /etc/snpeff/snpEff.config /usr/share/snpeff/snpEff.config "$HOME/.snpEff/snpEff.config"; do
    [ -f "$c" ] && SYS_CFG="$c" && break
  done
  [ -n "$SYS_CFG" ] || { echo "找不到系统 snpEff.config"; exit 1; }
  cp "$SYS_CFG" snpeff/snpEff.config
  sed -i 's|^data.dir.*|data.dir = ./data/|' snpeff/snpEff.config
  printf '\n# ---- demo genome ----\nDEMO.genome : %s\nDEMO.DEMO_chr1.codonTable : Bacterial_and_Plant_Plastid\n' \
    "Escherichia coli str. DEMO-WT (demo mini-chromosome)" >> snpeff/snpEff.config
fi

pushd snpeff > /dev/null
[ -f data/DEMO/snpEffectPredictor.bin ] || snpEff build -gff3 -v DEMO -c snpEff.config > ../work/vcf/snpeff_build.log 2>&1
snpEff ann -c snpEff.config -noLog \
    -csvStats ../work/vcf/snpeff_stats.csv \
    -stats ../work/vcf/snpeff_summary.html \
    DEMO ../work/vcf/joint.pass.vcf.gz 2> ../work/vcf/snpeff.log \
  | awk 'f||/^##fileformat/{f=1;print}' > ../work/vcf/joint.ann.vcf
popd > /dev/null
bgzip -f -c work/vcf/joint.ann.vcf > work/vcf/joint.ann.vcf.gz
bcftools index -f work/vcf/joint.ann.vcf.gz
bcftools query -f '%POS\t%REF\t%ALT\t%INFO/ANN\n' work/vcf/joint.ann.vcf.gz \
  | awk -F'\t' '{split($4,a,","); split(a[1],b,"|");
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,b[2],b[4],b[10],b[11]}' \
  > "$SNAP/08_snpeff_annotation.txt"
cat "$SNAP/08_snpeff_annotation.txt"

# -----------------------------------------------------------------------------
log "STEP 9 | 覆盖度统计"
# -----------------------------------------------------------------------------
for S in "${SAMPLES[@]}"; do
  samtools depth -a -q 0 -Q 20 "work/align/${S}.bam" > "work/cov/${S}.depth.txt"
  awk -v s="$S" '{t+=$3; n++; if($3>=10)c10++; if($3>=30)c30++}
       END{printf "  %-10s mean=%.1fx  >=10x=%.2f%%  >=30x=%.2f%%\n", s, t/n, 100*c10/n, 100*c30/n}' \
    "work/cov/${S}.depth.txt"
done

log "流程结束。中间数据快照: work/snapshots/"
