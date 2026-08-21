---
name: neoantigen-selection
description: Select the 20-34 patient-specific neoantigens for an individualized cancer vaccine and assemble them into one mRNA construct. Use when the task is to prioritize/rank/filter neoantigens from a patient's somatic variants, decide which mutations go into a personalized mRNA or peptide vaccine (the intismeran autogene / mRNA-4157 / autogene cevumeran class of product), predict HLA-I/HLA-II presentation of mutant peptides with NetMHCpan, compute agretopicity / foreignness / clonality / expression features, design a concatemeric minigene payload and check its junction epitopes, or benchmark a neoantigen ranking against validated epitopes. Also use for "which mutations should we put in the vaccine", "rank these neoepitopes", "build the vaccine construct", "how many neoantigens can we fit". Runs on open data only - cBioPortal, UniProt, IEDB NetMHCpan cloud, IEDB T-cell assays - no licensed binaries.
---

# Personalized neoantigen selection -> mRNA construct

The hard part of an individualized neoantigen therapy is not *generating*
candidates. A melanoma exome yields tens of thousands of mutant peptides in
minutes. The hard part is **spending a 34-slot payload**: choosing the
mutations most likely to be expressed by the tumor, presented by *this*
patient's HLA, seen by a T cell, and safely encodable in one mRNA.

This skill implements that selection layer end to end, with every threshold and
weight explicit and editable, and with a benchmark so the ranking can be
checked instead of trusted.

## What it does not do

It does not reproduce any company's proprietary score, training set,
thresholds or construct rules. The public description of intismeran autogene
(mRNA-4157/V940) fixes the *workflow* -- tumor/normal DNA + tumor RNA -> somatic
variants -> HLA typing -> presentation prediction -> <=34 neoantigens -> one
LNP-formulated mRNA -- and nothing more. Everything numeric here is a
literature-grounded stand-in, labelled as such in the output.

## The eight steps

| step | module | what it produces |
|---|---|---|
| 1-3 somatic variants, expression, clonality | `variants.py` | gated variant table + a gate waterfall |
| 4 mutant/WT peptide pairs | `peptides.py` | every window covering the mutation, plus its self counterpart |
| 5 HLA presentation | `presentation.py` | NetMHCpan-4.1 EL %rank per peptide x allele (IEDB cloud) |
| 6a features + score | `features.py`, `score.py` | 8 features in [0,1] -> one composite score |
| 6b constrained selection | `select.py` | <=34 slots under gene caps, allele spread, clonal preference |
| 7 construct | `construct.py` | junction-minimizing minigene order + codon-optimized CDS + QC |
| 8 report | `report.py` | markdown with `[computed]` / `[assumed]` / `[unverified]` labels |
| - benchmark | `benchmark.py` | AUC vs validated IEDB neoepitopes and matched decoys |
| - TESLA | `tesla.py` | AP / top-N recovery on 522 real T-cell-assayed pMHC pairs |

## Run the demo

```bash
python -m neoantigen_pipeline.selftest                       # offline, ~2 s, no network
python -m neoantigen_pipeline.run_demo --out demo_out --benchmark --tesla
```

Real TCGA-SKCM melanoma tumor (cBioPortal open API), real UniProt proteome,
real NetMHCpan calls, real IEDB T-cell-assay ground truth. First run downloads
and caches the proteome (~5 min); everything after that is cached.

## Run it on a patient

```python
from neoantigen_pipeline import fetch, pipeline, variants as V, report
from neoantigen_pipeline.config import PipelineConfig, PatientConfig

var = V.from_maf("patient.maf")                    # from your somatic caller
var = V.add_expression(var, tpm_by_entrez)         # tumor RNA-seq
var = V.add_clonality(var, purity=0.62)            # from ABSOLUTE / FACETS

cfg = PipelineConfig(patient=PatientConfig(
    "PT-014",
    hla_class1=["HLA-A*02:01","HLA-A*24:02","HLA-B*07:02",
                "HLA-B*44:02","HLA-C*05:01","HLA-C*07:02"],   # from OptiType/xHLA
    hla_class2=["HLA-DRB1*04:01"], tumor_purity=0.62))
cfg.selection.force_include_genes = ("BRAF","NRAS")           # mandate drivers
cfg.construct.linker = "GPGPG"                                # or "" for direct fusion

res = pipeline.run_pipeline(var, cfg)
pipeline.write_outputs(res, "PT-014_out")
open("PT-014_out/REPORT.md","w").write(report.build_report(res, cfg))
```

`res` keys: `gate_waterfall`, `peptides`, `peptides_skipped`, `predictions`,
`candidates`, `ranked`, `selected`, `coverage`, `selection_qc`, `construct`.

## The five things that actually decide the outcome

1. **HLA typing is an input, not a prediction.** Four-digit class-I typing from
   OptiType/xHLA/HLA-HD on the *normal* sample. A wrong allele silently
   invalidates the whole ranking. The demo has to assume one and says so.
2. **Expression is a gate, not a feature.** A beautifully-binding peptide from a
   gene at 0 TPM is not a neoantigen. Gate first, score second.
3. **Agretopicity needs the wild-type peptide.** Every mutant window is emitted
   with its positionally matched self counterpart so `WT rank / MUT rank` is
   exact rather than approximated.
4. **Junctions are new epitopes.** Fusing 34 minigenes creates 33 sequences that
   were never in the tumor. `construct.py` orders the minigenes to minimize
   predicted junction binding (greedy + 2-opt over a junction-cost matrix) and
   rescans the final junctions. Skipping this step is how a payload ends up
   with a strong non-tumor decoy epitope in it.
5. **Selection is constrained, not just sorted.** Gene caps, allele spread and
   clonal preference change which 34 you get. `select.py` records
   `why_selected` for every slot.

## What the benchmarks say about the score (read before quoting it)

Two benchmarks ship, and they agree.

**TESLA mirror** (`tesla.py`) — 522 peptide-HLA pairs, 35 experimentally
immunogenic, 6 patients, from the public DeepImmuno mirror. A label of 0 here
means *assayed and negative*, so these are real metrics, not lower bounds.
Random baseline AP = 0.067.

| score | AP | AUC | positives in a 34-slot budget |
|---|---|---|---|
| **NetMHCpan-4.1 EL %rank alone** | **0.207** | 0.791 | **31 / 35** |
| this package's composite | 0.149 | 0.729 | 24 / 35 |
| best published column in the mirror (`cnn_regress`) | 0.132 | 0.654 | 19 / 35 |
| DeepImmuno `immunogenic score` | 0.083 | 0.477 | 13 / 35 |

**IEDB-mined benchmark** (`benchmark.py`) — with decoys matched only on allele
and length, NetMHCpan alone scores AUC 0.966; that is a trap, not a result,
because IEDB epitopes were largely discovered *because* they bind. With binding
controlled for it falls to 0.599 and the composite to 0.576.

Both say the same thing: **a current presentation predictor is the strongest
single signal, and the extra peptide-intrinsic features do not beat it.** So the
default weights are presentation-dominant (0.45), `config.PRESENTATION_ONLY`
ships as a preset because it scored highest on TESLA, and
`config.LITERATURE_BALANCED` keeps the old set so the change is reproducible.

What neither benchmark can evaluate — and where the pipeline's real work
happens — is the expression / tumor-specificity / clonality gates and the
payload constraints. Full numbers, per-patient breakdown and caveats in
`reference/benchmark.md`.

## Backends

| layer | default | swap-in |
|---|---|---|
| MHC-I binding | NetMHCpan-4.1 EL via IEDB cloud REST | `mhcflurry` (auto-detected), local NetMHCpan |
| MHC-II binding | NetMHCIIpan via IEDB cloud REST | local NetMHCIIpan |
| immunogenicity cross-check | - | the `immunogenicity-multimodel` skill (DTU ImmunoGeNN pIRS) |
| codon optimization / QC | built into `construct.py` | the `codon-optimize-qc` skill for the full 20-enzyme panel + CAI |
| somatic calling upstream | not included | the `clairs-somatic` skill (tumor/normal BAM -> VCF) |
| delivery modelling downstream | not included | the `lnp-delivery-kinetics` skill (LNP uptake/escape/expression ODE) |

All prediction results are cached on disk per batch, so re-runs are free.

## Honesty rules (enforced in the report)

- `[computed]` -- produced in this run by a real predictor on real data.
- `[assumed]` -- a configuration stand-in (purity, an unavailable HLA type, a weight).
- `[unverified]` -- a prediction with no experimental confirmation. **Every
  immunogenicity claim this skill makes is `[unverified]`.**
- Variants that cannot be tiled (frameshift without transcript annotation,
  gene missing from the proteome, isoform mismatch) are *reported with the
  reason*, never silently dropped.
- Benchmark decoys are unlabelled, not verified negative, so reported AUCs are
  described as lower bounds.

## Reference

- `reference/scoring.md` -- every feature, its formula, its citation, and what
  would break it.
- `reference/workflow.md` -- the public workflow this mirrors, what is open vs
  proprietary, and where real patient data has to come from.
- `reference/benchmark.md` -- how the ground truth is mined and why the numbers
  are lower bounds.
