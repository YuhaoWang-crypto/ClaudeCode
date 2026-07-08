---
name: alphagenome-modal
description: Score genomic variants (regulatory / splice / promoter / eQTL effects) with AlphaGenome running on Modal GPU. Use when you need variant-effect predictions from DNA sequence — splice-site disruption, promoter/enhancer/5'UTR expression changes, chromatin accessibility (ATAC/DNase), CAGE, RNA-seq tracks — for a gene of interest, a VCF/variant list, or building a companion-diagnostic (CDx) enrichment panel. Deploys the genomicsxai AlphaGenome PyTorch port + gtca weights to Modal, fetches hg38 windows from UCSC, and returns per-track deltas with tissue mapping. NOT for coding missense/nonsense/frameshift effects — route those to the evo2-modal skill.
---

# AlphaGenome on Modal — variant-effect scoring

Validated, working deployment. Predicts how a DNA variant changes regulatory output (expression, splicing, chromatin) from sequence alone.

## When to use
- "What does this splice/promoter/5'UTR variant do to expression?"
- Score a list/VCF of regulatory variants for a CDx enrichment panel.
- Characterise a gene's regulatory/tissue landscape from sequence.
- **Do NOT** use for protein-coding LoF (missense/nonsense/frameshift) — AlphaGenome is structurally blind to those; use `evo2-modal` instead. (Empirically confirmed: a coding missense scored ~0 on every regulatory head.)

## Prerequisites
- Modal auth in env: `MODAL_TOKEN_ID`, `MODAL_TOKEN_SECRET` (Modal reads them automatically). `HF_TOKEN` for weight download.
- `pip install modal`. Outbound HTTPS may go through a proxy (`$HTTPS_PROXY`).
- Weights: `gtca/alphagenome_pytorch` on HF (ungated safetensors), architecture = `github.com/genomicsxai/alphagenome-pytorch` (NOT the identically-named PyPI `alphagenome-pytorch` — that gives a total key-mismatch). Input is one-hot `(B, S, 4)`, A/C/G/T = channels 0–3; window up to 2^20, validated at 131,072 bp; `organism_index=0` = human.

## How to run
The deployment is in `scripts/alphagenome_modal.py` (Modal app `alphagenome-pytorch`, Volume `alphagenome-weights`, A100-80GB, model loaded once per container).

```bash
modal run scripts/alphagenome_modal.py::download_weights          # CPU, one-time → commits weights to Volume
modal run scripts/alphagenome_modal.py::predict                   # GPU smoke test (forward pass, prints head shapes)
# batch a variant list (chrom,pos,ref,alt): edit the entrypoint or call batch_variant_effect(variants)
modal run scripts/alphagenome_modal.py::score
```

- `variant_effect_at(chrom, pos, ref, alt, window=131072)` — fetches the hg38 window from UCSC REST, builds ref/alt, returns per-head deltas.
- `batch_variant_effect(variants)` — loads the model ONCE and loops a table (VCF/list); use this for many variants (amortises model load + cold start).

## Output heads (human)
`rna_seq` (768 tracks), `cage` (640), `atac` (256), `dnase` (384), `procap` (128), `chip_tf` (1664, 128bp), `chip_histone` (1152, 128bp), `contact_maps` (64×64×28), `splice_sites` (5-way classifier), `splice_site_usage` (734), `splice_junctions`. Deltas = per-track mean|alt−ref| over the window; `@var` = |Δ| summed at the variant base. Splice classifier gives max|Δprob| (0–1).

## Tissue mapping
The port ships `data/track_metadata_human.parquet` + `TrackMetadataCatalog` → resolve top-track index to assay/biosample/ontology (e.g. INHBE top rna_seq/cage → liver, UBERON:0002107). `splice_sites` is class-indexed (Donor±/Acceptor±/None), not tissue-indexed.

## Interpreting scores
- Splice hits are highest-confidence: max|Δprob| near 1.0 on a canonical donor/acceptor + a matching expression shift = strong functional call.
- CAGE/RNA-seq deltas have large dynamic range — treat as **relative/directional** rankings.
- **Fidelity caveat:** community JAX→PyTorch port, not DeepMind's official serving stack. Confirm top CDx candidates against the **official AlphaGenome API** before decision-grade use.

## Cost
CPU download (free-ish) + ~0.8 s/forward on A100; a batch of ~20 variants is a few minutes of GPU (<$5). `references/example_cdx_variants.tsv` is a worked 18-variant INHBE/ACVR1C/GPR75 example.
