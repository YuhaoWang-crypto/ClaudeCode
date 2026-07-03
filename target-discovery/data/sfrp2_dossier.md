# SFRP2 — target dossier (Geneformer #1 hit)

SFRP2 (secreted frizzled-related protein 2) was ranked **#1** by the Geneformer
in-silico deletion (KO pushes IPF fibroblast/epithelial cells toward normal).
The literature independently validates it as a **causal driver of the exact
aberrant-basaloid biology** we flagged as whitespace — a strong convergence
between the model's nomination and human functional data.

*Based on articles retrieved from PubMed.*

## Genetic / functional evidence

- **sFRP2 is a TGF-β1 fibroblast target that drives basal metaplasia in IPF.**
  In IPF/ILD, fibroblast-secreted sFRP2 is induced near AEC2 cells and activates
  a mature **KRT5+ basal-cell program** via Frizzled-5 → calcineurin/NFATc3 —
  i.e. it drives the KRT17+ → KRT5+ **aberrant basaloid** transition, our exact
  whitespace compartment. sFRP2 was *necessary* for KRT5 induction in
  AEC2-fibroblast organoids and precision-cut lung slices. Fibroblast-selective
  TGF-β1 inhibition (EGCG) reduced sFRP2 and IPF transcriptional changes.
  Cohen et al., *J Clin Invest* 2024
  ([DOI](https://doi.org/10.1172/JCI174598); preprint
  [DOI](https://doi.org/10.1101/2023.08.02.551383)); commentary Burgy &
  Königshoff ([DOI](https://doi.org/10.1172/JCI183970)).
- **SFRP2 marks fibroblast progenitors of myofibroblasts.** In systemic
  sclerosis skin/lung fibrosis, myofibroblasts arise from an SFRP2/DPP4+
  progenitor fibroblast population (scRNA-seq). Tabib et al., *Nat Commun* 2021
  ([DOI](https://doi.org/10.1038/s41467-021-24607-6)).
- **Context nuance (opposite direction in heart).** In cardiac fibrosis sFRP2
  can act as an *anti*-fibrotic Wnt antagonist (inactivating extracellular Wnt,
  preventing TGF-β myofibroblast transformation). Blyszczuk et al.,
  *Eur Heart J* 2017 ([DOI](https://doi.org/10.1093/eurheartj/ehw116)). →
  sFRP2 biology is tissue-context-dependent; a lung-directed biologic is the
  rational modality.

## Why it fits the pipeline
- **Model → biology convergence**: Geneformer nominated SFRP2 #1 from an unbiased
  candidate set; human functional studies show sFRP2 is *causally required* for
  the aberrant-basaloid transition we defined as the IPF whitespace.
- **Novelty**: not among the 55 in-clinic IPF Ph2/3 programs (crowded set).
- **Modality**: secreted protein; ChEMBL has **zero small-molecule ligands** →
  **biologic** (neutralizing antibody / nanobody / trap). See the Boltz nanobody
  design run (`data/sfrp2_nanobody_design.csv`).
- **Mechanistic hypothesis to test**: a lung-delivered anti-sFRP2 biologic blocks
  the fibroblast→AEC2 paracrine axis, preventing KRT5+ basaloid metaplasia.

## Target for biologic design (UniProt Q96HF1)
Mature sFRP2 (signal peptide trimmed) = cysteine-rich domain (Frizzled-like) +
netrin domain; the CRD is the Wnt/Fzd-interaction surface and the rational
epitope for a blocking binder.
