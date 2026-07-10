# perturbomics — perturbation-omics integration & combination analysis

A Claude Code **skill** (and a working Python package) that unifies drug, CRISPR,
single-cell, and foundation-model perturbations into one comparable object — a
signed ranked gene **signature** — then scores connectivity, mines drug+CRISPR
**combinations**, and fuses transcriptomic hits with pathway dynamics and live
drug-discovery data into a multi-scale **lead ranking**. Every claim is labelled
✅ rigorous (deterministic computation) vs ⚠️ hypothesis (needs validation).

## Install

Drop the `perturbomics/` folder into a project's `.claude/skills/` directory
(Claude Code auto-discovers it), or use the package standalone:

```bash
pip install numpy pandas scipy          # core; the two offline demos need only these
# optional: pip install pydeseq2 decoupler scanpy anndata   # real single-cell DGE
cd assets && python3 -m perturbomics.demo
```

## What's inside

```
SKILL.md                     the skill entry point (read this first)
README.md                    this file
reference/
  data-access.md             how to get each dataset (Repurposing Hub, CMap/clue.io,
                             GPP CRISPick, Geneformer) + the enrichment MCP servers
  dge-signatures.md          pseudobulk single-cell DGE → signature (the rigorous way)
  connectivity.md            GSEA enrichment score, CMap WTCS, NCS/τ, null models
  geneformer.md              embeddings + in-silico perturbation → signature
  integration.md             perturbomics × network-biomarker × drug-MCP funnel
assets/
  perturbomics/              the importable package
    signature.py             the common object + constructors from every source
    connectivity.py          enrichment_score, weighted_connectivity_score, NCS
    combine.py               rank_reversers, best_combinations (cross-modality)
    pseudobulk.py            raw counts → pseudobulk → signature (PyDESeq2 or fallback)
    enrichr.py               load REAL perturbation libraries (Enrichr, no login)
    integrate.py             NetworkContext, DrugEvidence, integrated_leads
    demo.py                  offline synthetic end-to-end (instant)
    demo_integrate.py        offline 4-axis funnel (instant)
    realdata_ipf.py          REAL run: Enrichr disease/drug/CRISPR → reversers+combos
  requirements.txt
examples/
  real_leads_ipf.py          REAL integrated leads (live ChEMBL + ClinicalTrials)
figures/
  connectivity_report.html   self-contained report of the runs (open in a browser)
```

## The five entry points

| Command | What it does | Needs network? |
|---|---|---|
| `python3 -m perturbomics.demo` | synthetic drug+CRISPR+DGE end-to-end | no |
| `python3 -m perturbomics.demo_integrate` | synthetic 4-axis lead funnel | no |
| `python3 -m perturbomics.realdata_ipf <cache> "<disease>"` | real Enrichr reversers + combos for any CREEDS disease | downloads ~45 MB |
| `python3 examples/real_leads_ipf.py` | real integrated leads (live ChEMBL/trials evidence) | live MCP |
| open `figures/connectivity_report.html` | the written report of all runs | no |

## The one idea

Every perturbation → a `Signature` (signed ranked genes). Then **one** score
(Connectivity-Map WTCS) compares any two — drug↔drug, drug↔CRISPR,
perturbation↔disease — so you can rank reversers, mine combinations, and fuse
with network control + druggability + clinical status. See `SKILL.md`.

## License / data provenance

Code: use freely within your project. Data sources are third-party and carry
their own terms — Drug Repurposing Hub (CC-BY 4.0), LINCS/CMap, Broad GPP,
Geneformer (Apache-2.0), Enrichr libraries, ChEMBL, ClinicalTrials.gov. Cite the
originals; see `reference/data-access.md`.
