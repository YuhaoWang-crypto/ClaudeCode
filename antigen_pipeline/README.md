# cell-surface-antigen-discovery

Nominate antibody-accessible tumour surface antigens for ADC and CAR-T
development in lung adenocarcinoma, from public data only, with the
validation set locked before anything is ranked.

```bash
pip install numpy scipy pandas matplotlib requests openpyxl reportlab pillow cellxgene-census
python3 -m antigen_pipeline.run_all          # steps 1-9, reuses cached downloads
python3 -m antigen_pipeline.s6_score         # or any single step, off the files on disk
python3 -m antigen_pipeline.s2_census --refetch   # force a fresh Census pull
```

## The one idea

Targets are ranked by **tumour surface specificity x normal-tissue therapeutic
index**, never by tumour essentiality.

An antibody therapeutic kills cells that display the antigen. Whether the
tumour *needs* the gene is irrelevant to efficacy and actively misleading as a
filter, because an essential gene is essential in normal tissue too. The run
confirms it: every clinically validated antigen in the validation set is
non-essential in DepMap (TROP2 -0.10, c-MET -0.12, HER2 -0.34). Gating on
essentiality would drop exactly the targets that work and enrich for
housekeeping genes. DepMap is therefore carried as annotation and never enters
the score.

```
score = geometric_mean( tumour_quality , safety_coefficient , consensus_multiplier )

tumour_quality = 0.30 relative TME specificity
               + 0.20 expression intensity
               + 0.20 expression uniformity
               + 0.15 ectodomain accessibility
               + 0.15 antibody druggability
```

Safety is multiplicative, so a target that is loud in a vital organ is pushed
down however good its tumour profile looks. The geometric mean is a monotone
transform of the product: it cannot reorder anything, it only puts the score on
the same 0-1 scale as its three factors, which is the scale the absolute tier
thresholds are stated against.

## Steps

| Step | Module | What it does |
|---|---|---|
| 1 | `s1_surfaceome` | SURFY in-silico surfaceome -> 2,799 surface genes, each with its own topology; **locks the validation set before any ranking** |
| 2 | `s2_census` | CZ CELLxGENE Census, whole-cell LUAD atlases, fixed-seed sampling, four TME compartments, cross-atlas consensus |
| 3 | `s3_topology` | extracellular topology gate: drops ER/Golgi/endosome/cytoplasmic/secreted-only proteins and anything with no usable ectodomain |
| 4 | `s4_druggability` | Open Targets antibody tractability, independent localisation, known drugs; DepMap essentiality as annotation only |
| 5 | `s5_safety` | HPA consensus RNA + normal-tissue IHC, organ-weighted, **conservative minimum of the two arms** |
| 6 | `s6_score` | antigen display threshold, composite score, absolute tiers |
| 7 | `s7_validate` | recall of the pre-registered antigens, negative-control PASS/FAIL, rank stability under six alternative schemes |
| 8 | `s8_literature` | Europe PMC triage of the head of the ranking into four evidence classes, with an impossible-symbol calibration query |
| 9 | `s9_figures`, `s9_export`, `s9_report` | four figures (PNG+SVG), eight consistency gates, English and Chinese PDF reports |

## Design decisions that matter

**Single-nucleus data is excluded.** Nuclear preparations systematically
under-detect membrane-protein transcripts, which is precisely the signal this
pipeline ranks on.

**The safety coefficient is the minimum of the protein and RNA arms**, not the
mean. A target that looks clean in RNA but stains in heart muscle is treated as
dangerous. Genes with no HPA record get a flagged neutral default rather than a
free pass or a silent drop.

**No organ can zero a target outright.** Organ weights are capped below 1, so
the safety coefficient has a floor. A score of exactly zero would make the rest
of the evidence unreadable, and toxicity that is managed in the clinic is a
cost, not an absolute veto.

**An antigen must actually be displayed.** Candidates detected in under 10% of
epithelial/malignant cells are set aside before ranking. Without that rule a
gene no tumour cell expresses wins on a safety coefficient earned purely by
being silent everywhere.

**Gene symbols are pinned to title/abstract in the literature queries.**
Unrestricted full-text matching over-counts by one to two orders of magnitude,
and a paper that merely cites a sequence accession is not evidence about a
target.

**The report reads every number back from disk.** `s9_report` loads the result
files and provenance sidecars; it never recomputes a statistic or carries one
over from the run that produced it, so a report regenerated later must agree
with the results directory.

## Outputs

- `results/` — one CSV plus a provenance JSON per step, and `results/export/`
  with the bundle and its `manifest.json` (gate outcomes and per-file checksums)
- `figures/antigen/` — four figures as PNG and SVG
- `reports/antigen_discovery_report_en.pdf`, `..._zh.pdf`

## Data sources

All public, all fetched at run time into `data/raw/` (git-ignored):
SURFY surfaceome (Bausch-Fluck et al., PNAS 2018), CZ CELLxGENE Census LTS
release, Human Protein Atlas consensus RNA and normal-tissue IHC, Open Targets
Platform GraphQL, DepMap CRISPR gene effect, Europe PMC.
