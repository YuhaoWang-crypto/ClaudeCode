# Methodology — the layered framework and scoring math

## Why layers, not a single similarity number

"Most similar to human" is the wrong target. Two organisms can share every
pathway gene yet behave differently because the *wiring*, *regulatory logic*, or
*dynamics* differ; and the best organism for a question is often the one that is
deliberately *unlike* human (a natural knockout, a regeneration champion). So
PSMOS scores five to six semi-independent layers and reports role-specific
winners.

## The six layers (PSMOS G/D/N/R/E/X)

```
PSMOS(s, P) = wG·G + wD·D + wN·N + wR·R + wE·E + wX·X
```

| Layer | Question | How PSMOS fills it |
|---|---|---|
| **G** | Are the core genes present and sequence-conserved? | ortholog search (gate) + **Evo2** log-likelihood constraint |
| **D** | Are domains / active sites / interfaces conserved? | curated (Evo2 partial); productionise with domain-restricted Evo2 |
| **N** | Is the network topology + feedback preserved? | curated (topology edit distance) |
| **R** | Is the regulatory grammar (promoter/enhancer/TF/splice) preserved? | **AlphaGenome** (human/mouse); curated elsewhere |
| **E** | Does it express in the right tissue/stage/cell? | curated; productionise with single-cell / AlphaGenome expression |
| **X** | Can you actually do the experiment (tools/throughput/cost)? | curated priors |

Weights change with the question. The package ships profiles: `default`,
`regeneration`, `drug-target`, `enhancer/cis-reg`, `development`,
`signal-dynamics`. Re-weighting reorders the ranking — that is the feature, not
a bug (`build_hippo_dashboard.py` shows this as tabs).

The Notch pilot uses a lighter 6-axis variant (completeness, fidelity,
simplicity, arch_match, tract, thru) with three research goals; Hippo uses the
full G/D/N/R/E/X. Both share the gate, Evo2, and Compara machinery.

## The hard gate (done honestly)

```
found_gate = gate_families ∩ orthologs_found(species)
if |found_gate| == |gate_families|:  PASS  (computed:uniprot)
elif |found_gate| == 0 and curated_completeness == 0:
                                     FAIL  (computed:uniprot, confirmed absent)  # yeast/plant
elif |found_gate| == 0:              curated fallback (no hit — likely query miss)
else:                                curated fallback (annotation gap), flagged
```

The middle branches encode *absence of evidence ≠ evidence of absence*. A true
negative control (yeast/plant on Notch) needs **both** an empirical zero-hit AND
a curated prior of absence. Everything else that misses is treated as an
annotation gap.

## Fidelity vs utility (the scatter)

Two composite axes make "no single best organism" visual:
- **Biological fidelity** = f(G, D, N, R, E) — how faithfully the pathway is
  reproduced.
- **Experimental utility** = X — how doable the experiment is.

Human sits top-left (max fidelity, min utility — the target, not a model); the
useful models spread along the trade-off. There is no top-right point.

## Evo2 log-likelihood — what it is and is NOT

- **Is**: mean per-token log-likelihood of a CDS under a genome model trained
  across the tree of life → "sequence constraint / naturalness". Higher (less
  negative) = more constrained/natural.
- **Is NOT**: "similarity to human". Evo2 is pan-species. So it is reported as an
  independent axis, and its Pearson correlation with the curated divergence-time
  fidelity proxy is stated (Notch r=+0.73 over 9 species). r≈0.7 is the useful
  regime: it corroborates direction but adds signal — e.g. Evo2 rates fly Notch
  coding sequence far more constrained than divergence time implies.
- Normalise across the live species (min-max on mean-LL) to a 0..1 constraint
  score for display; keep the raw meanLL in the tooltip for auditability.

## AlphaGenome R — human↔mouse regulatory concordance

Only human and mouse are in-scope for AlphaGenome. R(gene) = mean across
RNA_SEQ/ATAC/CHIP_TF of the Pearson correlation between the TSS-anchored,
strand-oriented predicted profiles of the human and mouse orthologous loci,
mapped to 0..1. Human is the R reference (1.0); mouse's R layer is the computed
concordance. When this replaces a curated prior and moves the ranking (Hippo:
0.90→0.614, mouse drops below human), that is a reported finding.

## Provenance is data, not decoration

Every axis value carries a provenance string
(`computed:{uniprot,ensembl,evo2,compara,alphagenome}` or `curated`). The
dashboards render it (legend + `*` on computed cells + tooltips). The rule: a
reader must always be able to tell a measurement from a prior. When you extend a
layer from curated to computed, flip the tag and, if it changes a conclusion,
say so in the footer.
