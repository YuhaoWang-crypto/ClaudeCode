"""
PSMOS — Pathway-aware model-organism selection.

Given a signalling pathway, recommend the right model organism(s) — per research
purpose — using a layered, provenance-labelled score. The sequence-constraint
layer is computed by the Evo2 genome foundation model (run on Modal H100); the
hard gate is an empirical ortholog search (UniProt + Ensembl); remaining axes
are explicitly-labelled comparative-genomics priors.

Modules: pathways, orthologs, cds, evo2_modal, run_evo2_scoring, scoring,
build_dashboard, run_all.
"""
