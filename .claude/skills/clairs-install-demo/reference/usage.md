# ClairS usage reference

Entry point: `run_clairs` (a Python script in the ClairS checkout). Sub-commands
of the pipeline live behind `python3 clairs.py <module>`; the two you call
directly are `compare_vcf` (benchmarking) and, rarely, `cal_metrics_in_af_range`.

## Required

```
-T, --tumor_bam_fn    tumour BAM, samtools-indexed
-N, --normal_bam_fn   normal BAM, samtools-indexed
-R, --ref_fn          reference FASTA, samtools-indexed
-o, --output_dir      output directory
-t, --threads         max threads
-p, --platform        see the platform table in SKILL.md
```

Output is `${output_dir}/output.vcf.gz` (+ `.tbi`), plus `run_clairs.log`,
`logs/` (per-step logs) and `tmp/` (intermediates — `--remove_intermediate_dir`
to drop them).

## Scoping the run

| Flag | Use |
|---|---|
| `-c/--ctg_name chr21,chr22` | whole contigs |
| `-r/--region chr20:1000000-2000000` | one region, 1-based start |
| `-b/--bed_fn regions.bed` | many regions — preferred for panels/exomes |
| `--include_all_ctgs` | otherwise only chr1–22 (+X/Y) are processed |
| `--chunk_size` | default 5 Mb; lower it when a contig is short and parallelism is poor |

## Calling modes

- **de novo** (default) — call anywhere in scope.
- **genotyping**, `-G/--genotyping_mode_vcf_fn sites.vcf` — report only at the
  given sites. Use to force-call a hotspot panel or to compare cohorts at fixed
  loci.
- **hybrid**, `-H/--hybrid_mode_vcf_fn sites.vcf` — de novo results *plus*
  genotyping at those sites, in one pass.

## Sensitivity / quality knobs

| Flag | Default | Effect |
|---|---|---|
| `--snv_min_af` | 0.05 | lower → more recall in the low-VAF tail, worse precision + slower |
| `--indel_min_af` | 0.1 (ONT) | same trade-off for Indels |
| `--min_coverage` | 4 | minimum depth to consider a site |
| `--snv_min_qual` / `--indel_min_qual` | model config | QUAL above which a call is `PASS` rather than `LowQual`; supersedes the deprecated `--qual` |
| `--print_ref_calls` / `--print_germline_calls` | off | emit `RefCall` / `Germline` records too |

Since v0.4.3 a model directory may ship `model_specific_settings.conf` carrying
`snv_min_qual=` / `indel_min_qual=`; those are read automatically and are why
two platforms can behave differently at identical CLI flags.

## Indels

`--enable_indel_calling`. Supported on ONT R10 and PacBio HiFi; **not** on
`ont_r9_guppy` or `ilmn`. Calling time rises substantially. Indel F1 ≈ 73 % at
50×/50× HCC1395/BL. Output goes to `${output_dir}/indel.vcf.gz` by default
(`--indel_output_prefix`).

## Phasing

ClairS phases germline hets and haplotags the tumour reads; the `H` INFO flag on
a call means it was seen on a single haplotype. Defaults (v0.3.1+):

- hets from the **normal** sample, applied to the **tumour** BAM
- LongPhase for both phasing and haplotagging (`--use_longphase_for_intermediate_haplotagging`)
- heterozygous Indels included (`--use_heterozygous_indel_for_intermediate_phasing`)

Switch the source with the four
`--use_heterozygous_snp_in_{normal,tumor}_sample_and_{normal,tumor}_bam_for_intermediate_phasing`
flags. If normal coverage is low, prefer `..._in_tumor_sample_...`. If the
tumour BAM is already haplotagged (WhatsHap or LongPhase),
`--haplotagged_tumor_bam_provided_so_skip_intermediate_phasing_and_haplotagging`
skips the whole stage. `--disable_phasing` exists but costs real accuracy.

`--whatshap_for_phasing` reverts to WhatsHap (slower).

## Skipping work you already have

- `--normal_vcf_fn normal.vcf.gz` — skip Clair3 germline calling on the normal.
- `--enable_clair3_germline_output` — use Clair3's *default* (not fast) settings
  and emit both germline VCFs; ~40 % slower.

## Verdict

`--enable_verdict` tags calls as germline / somatic / subclonal somatic using a
CNV profile plus a tumour-purity estimate. Suggested only for purity < 0.8.
Requires the extra install step in `reference/installation.md`. Docs:
`docs/verdict.md` in the repo.

## LongPhase-S post-filter

Separate from ClairS itself: LongPhase-S reconstructs somatic haplotypes and
re-flags calls inconsistent with them as `LowQual`. See
`docs/longphase-s_post-filter.md`. Worth it on ONT tumour-normal where FP rate
matters more than the last point of recall.

## Output VCF fields

| Field | Meaning |
|---|---|
| `FILTER` | `PASS` somatic · `LowQual` · `RefCall` · `Germline` |
| `AF` `DP` `AD` | tumour allele fraction, depth, (ref,alt) depths |
| `NAF` `NDP` `NAD` | same in the matched normal — near 0 for a true somatic call |
| `AU` `CU` `GU` `TU` | per-base read counts, tumour |
| `NAU` `NCU` `NGU` `NTU` | per-base read counts, normal |
| `FAU`…`FTU` / `RAU`…`RTU` | forward / reverse strand splits in tumour (strand-bias check) |
| `H` (INFO) | variant confined to one haplotype in the phased reads |
| `GQ` | genotype quality |

A quick sanity filter for a candidate list:

```bash
bcftools view -f PASS out/output.vcf.gz \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF\t%NAF\t%DP\t%NDP]\n' \
  | awk '$6 < 0.01 && $7 >= 10'
```

## Benchmarking

```bash
python3 "$CLAIRS_HOME/clairs.py" compare_vcf \
  --truth_vcf_fn truth.vcf.gz \
  --input_vcf_fn out/output.vcf.gz \
  --bed_fn high_confidence.bed \
  --output_dir out/benchmark \
  --input_filter_tag PASS \
  --ctg_name chr17 --ctg_start 80000000 --ctg_end 80100000
```

Prints precision / recall / F1 / TP / FP / FN and writes `tp.vcf`, `fp.vcf`,
`fn.vcf`, `fp_fn.vcf`. Restricting to the truth set's high-confidence BED is not
optional — outside it, "FP" is meaningless.

For SEQC2-comparable numbers use `som.py` from `jmcdani20/hap.py:v0.3.12` with
`-T target.bed -f high_confidence.bed -r ref.fa`.

## Related tools

| Tool | For |
|---|---|
| Clair3 | germline small variants, DNA-seq |
| Clair3-RNA | germline small variants, long-read RNA-seq |
| ClairS-TO | somatic calling from a **tumour-only** sample (no matched normal) |
| ClairS-series-model-training | training your own `ss` / `ssrs` models |
