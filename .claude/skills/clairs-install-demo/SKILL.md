---
name: clairs-install-demo
description: >-
  Install ClairS (HKU-BAL tumour-normal somatic small-variant caller) natively
  on Linux without Docker — micromamba from source — and run/benchmark the
  HCC1395 quick demo on ONT, PacBio HiFi or Illumina. Use when there is no
  Docker daemon (agent sandboxes, CI runners, HPC login nodes), when a
  ClairS/Clair3 install is failing (`pileup.pt` not found, missing Boost
  headers, 403 on model downloads), when reproducing the SEQC2/HCC1395 chr17
  demo and its expected precision/recall, when picking the right `--platform`
  model, or when reading ClairS VCF tags (AF/NAF, AU/CU/GU/TU, H). For running
  the official Docker image on Modal or another Linux host instead, see the
  `clairs-somatic` skill.
---

# ClairS without Docker — native install + validated demo

ClairS calls somatic SNVs (and, on some platforms, Indels) from a **paired
tumour + normal BAM**. It runs Clair3 on both BAMs to establish germline
background, phases and haplotags the tumour reads, then ensembles a pileup and
a full-alignment neural network over tumour/normal paired tensors, and applies
haplotype + post-calling filters.

- Repo: <https://github.com/HKU-BAL/ClairS> · paper: Nat. Methods (2026)
- Tumour-**only** samples → ClairS-TO. Germline-only → Clair3. Long-read RNA → Clair3-RNA.

## Decide the route first

| Situation | Route |
|---|---|
| Docker daemon available | `docker run hkubal/clairs:latest /opt/bin/run_clairs ...` — fastest, nothing to build |
| macOS, or you want it on Modal | the **`clairs-somatic`** skill — dispatches the official image to a Linux host |
| HPC, no root | Singularity: `singularity pull docker://hkubal/clairs:latest` |
| Linux sandbox/CI, **no Docker daemon** | **this skill** — `scripts/install_clairs.sh`, micromamba from source |

The rest of this skill is the last row: a native Linux install with no
container at all, and the demo that proves it matches the official image.

## Install from source (no Docker)

```bash
bash scripts/install_clairs.sh ~/clairs      # ~10 min, ~8 GB
source ~/clairs/env.sh                       # puts run_clairs on PATH
```

`~/clairs` is also the prefix `run_demo.sh` assumes; install anywhere else and
export `CLAIRS_PREFIX` to match.

The script pins micromamba, creates the conda env, clones ClairS, compiles the
two C++ helpers, and unpacks the 832 MB model bundle into
`$CONDA_PREFIX/bin/clairs_models/` (the only place `run_clairs` looks).

### Three failure modes worth knowing before you start

1. **`clair3` must be pinned to 1.x.** Plain `conda install clair3` now resolves
   to 2.x, which is PyTorch-based and wants `pileup.pt`; the Clair3 models
   inside ClairS's `clairs_models.tar.gz` are TensorFlow checkpoints
   (`pileup.index` + `pileup.data-*`). Result partway through a run:
   `FileNotFoundError: .../clair3_models/ont_r104_e81_sup_g5015/pileup.pt`.
   Use `clair3=1.2.0`. ClairS itself is PyTorch — both frameworks coexist in
   the one env. (Upstream's Dockerfile says just `clair3`; it was written when
   `clair3` still meant 1.x.)
2. **Boost headers** are needed for `debruijn_graph` (upstream apt-installs
   `libboost-graph-dev`); from conda use `libboost-headers` plus
   `-I$CONDA_PREFIX/include`.
3. **Model/demo downloads.** Upstream docs use `http://www.bio8.cs.hku.hk/...`.
   Rewrite to `https://` — identical files, and it survives egress policies
   that only allow 443.

Full detail, including the Verdict extras: `reference/installation.md`.

## Run

```bash
run_clairs \
  --tumor_bam_fn  tumor.bam \      # samtools-indexed
  --normal_bam_fn normal.bam \
  --ref_fn        ref.fa \         # samtools-indexed
  --threads 8 \
  --platform ont_r10_dorado_sup_5khz_ssrs \
  --output_dir out
# -> out/output.vcf.gz  (PASS = somatic call)
```

Restrict work with `--region chr17:80000000-80100000`, `--ctg_name chr21,chr22`,
or `--bed_fn regions.bed`. Add `--enable_indel_calling` for somatic Indels
(ONT R10 / HiFi only; substantially slower).

### Picking `--platform`

Choose by chemistry **and** basecaller — the wrong model silently costs recall.

| Data | `--platform` |
|---|---|
| ONT R10.4.1 5 kHz, Dorado SUP | `ont_r10_dorado_sup_5khz_ssrs` (default choice) |
| …same, cancer type absent from training a concern | `ont_r10_dorado_sup_5khz_ss` |
| ONT R10.4.1 4 kHz, Dorado SUP | `ont_r10_dorado_sup_4khz` |
| ONT R10.4.1 5 kHz / 4 kHz, Dorado HAC | `ont_r10_dorado_hac_5khz` / `_4khz` |
| ONT R10.4 Guppy5 SUP | `ont_r10_guppy` |
| ONT R9.4.1 Guppy5 SUP | `ont_r9_guppy` |
| PacBio Revio (SMRTbell 3.0) | `hifi_revio_ssrs` |
| PacBio Sequel II | `hifi_sequel2` (experimental; downsample to ~40×) |
| Illumina NovaSeq/HiSeqX | `ilmn_ssrs` |

`ssrs` = synthetic + real-sample training (better, use by default); `ss` =
synthetic only. Reference should be GRCh38_no_alt for long reads, GRCh38 for
Illumina; models were trained with chr20 held out.

Every other flag, plus genotyping/hybrid modes, phasing switches, Verdict and
tuning advice: `reference/usage.md`.

## Reading the output VCF

`FILTER` is `PASS` (somatic), `LowQual`, `RefCall`, or `Germline`. Per-record:
`AF`/`DP`/`AD` are tumour, `NAF`/`NDP`/`NAD` the matched normal — a real
somatic call has `NAF` at or near 0. `AU/CU/GU/TU` and `NAU/...` give per-base
counts in tumour/normal; `FAU/RAU/...` split them by strand (use them to spot
strand bias). The `H` INFO flag means the variant sits on a single haplotype in
the phased reads.

## Demo — HCC1395 / HCC1395BL, chr17:80.0–80.1 Mb

```bash
bash scripts/run_demo.sh ont ~/clairs_demo/ont 4    # ont | ilmn | pacbio_hifi
```

Fetches the tumour/normal BAM pair, the chr17 reference and the SEQC2 v1.2
truth VCF + high-confidence BED, calls, then benchmarks with
`clairs.py compare_vcf`. **Measured** on 4 CPU threads, 29 truth SNVs:

| Dataset | `--platform` | Precision | Recall | F1 | TP | FP | FN | wall |
|---|---|---|---|---|---|---|---|---|
| ONT R10.4.1 | `ont_r10_guppy` | 1.0 | 0.9655 | 0.9825 | 28 | 0 | 1 | 1m40s |
| Illumina NovaSeq | `ilmn` | 1.0 | 1.0 | 1.0 | 29 | 0 | 0 | 0m25s |
| PacBio Revio | `hifi_revio` | 1.0 | 0.9655 | 0.9825 | 28 | 0 | 1 | 1m38s |

All three reproduce upstream's published expected output exactly — that match
is the check that a from-source install is sound. Both long-read runs miss the
same site (`chr17:80,094,483`, TVAF≈0.096), the lowest-VAF truth variant in the
window; Illumina gets it at AF=0.064. Details: `reference/demo-results.md`;
committed artefacts under `clairs_demo/`.

## Benchmarking your own calls

```bash
python3 $CLAIRS_HOME/clairs.py compare_vcf \
  --truth_vcf_fn truth.vcf.gz --input_vcf_fn out/output.vcf.gz \
  --bed_fn high_confidence.bed --output_dir out/benchmark \
  --input_filter_tag PASS --ctg_name chr17
```

Writes `tp.vcf` / `fp.vcf` / `fn.vcf` next to the metrics — always read `fn.vcf`
before concluding a model is underperforming; most misses are low-VAF.
`som.py` (hap.py image) is the alternative when you need SEQC2-comparable
numbers.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `FileNotFoundError: .../pileup.pt` | clair3 2.x installed — pin `clair3=1.2.0` |
| `Cannot find clair3 main entry in ...` | env's `bin/run_clair3.sh` missing; `--clair3_path` points elsewhere |
| `Cannot find clair3 model path in ...` | model bundle not unpacked to `$CONDA_PREFIX/bin/clairs_models` |
| `boost/graph/adjacency_list.hpp: No such file` | install `libboost-headers`, compile with `-I$CONDA_PREFIX/include` |
| 403 fetching models/demo data | you used `http://`; switch to `https://` |
| `ModuleNotFoundError: sklearn` | only Verdict needs it; `conda install scikit-learn` |
| Empty VCF, no error | region has no candidates, or `--platform` mismatches the data |

## Scripts

| Script | Purpose |
|---|---|
| `scripts/install_clairs.sh [PREFIX]` | micromamba env + ClairS + compiled helpers + models |
| `scripts/fetch_demo_data.sh <ont\|ilmn\|pacbio_hifi> [DIR]` | official HCC1395 quick-demo data |
| `scripts/run_demo.sh <platform> [DIR] [THREADS]` | fetch → call → benchmark, end to end |

HKU-BAL also publishes agent skills for the Clair family at
<https://github.com/HKU-BAL/Clair-skills>.
