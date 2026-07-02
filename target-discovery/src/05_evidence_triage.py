"""
Step 5 — Evidence triage + target dossier.

The ranked, novelty-gated targets from Step 4 are handed to the MCP evidence
layer. MCP tools (ChEMBL, PubMed, ClinicalTrials.gov, Boltz) are invoked by
Claude, not from inside this process, so this script does two things:

  1. Emits a machine-readable *triage worklist* (`results/triage_worklist.json`)
     naming, per gene, exactly which MCP calls to run.
  2. Writes a dossier skeleton (`results/target_dossier.md`) for Claude to fill,
     applying the SAME template already validated on the seed shortlist
     (data/seed_shortlist.csv).

The per-gene evidence protocol (what Claude runs):
  - ChEMBL  target_search(gene_symbol)      -> target class, PDB count, ligandability
  - ChEMBL  get_bioactivity(target_chembl_id)-> known chemical matter
  - PubMed  search_articles("<gene> pulmonary fibrosis") -> functional/genetic evidence
  - ClinicalTrials search_trials(intervention=<gene/drug>) -> confirm novelty
  - Boltz   structure_and_binding           -> pocket / binder tractability (top tier only)
Scoring: genetic_evidence x expression_specificity x druggability x novelty x safety.
"""
from __future__ import annotations
import argparse
import json
import os

import pandas as pd
import yaml


def load_config(path: str) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def build_worklist(genes: list[str], disease: str) -> list[dict]:
    work = []
    for g in genes:
        work.append(
            {
                "gene_symbol": g,
                "mcp_calls": [
                    {"tool": "ChEMBL.target_search", "args": {"gene_symbol": g}},
                    {"tool": "ChEMBL.get_bioactivity",
                     "args": {"note": "use target_chembl_id from target_search"}},
                    {"tool": "PubMed.search_articles",
                     "args": {"query": f"{g} {disease} fibroblast"}},
                    {"tool": "ClinicalTrials.search_trials",
                     "args": {"intervention": g, "condition": disease}},
                    {"tool": "Boltz.start_structure_and_binding",
                     "args": {"note": "top-tier only; assess pocket/binder"}},
                ],
                "score_axes": [
                    "genetic_evidence", "expression_specificity",
                    "druggability", "novelty", "safety",
                ],
            }
        )
    return work


DOSSIER_HEADER = """# {disease} — Target Dossier

Generated from Geneformer in-silico perturbation (Step 3), novelty-gated (Step 4),
and evidence-triaged via MCP (Step 5). One page per target; ranked most-promising
first. Fill each section from the MCP results named in triage_worklist.json.

---
"""

DOSSIER_ENTRY = """## {rank}. {gene}  (gated score: {score})

- **Mechanistic hypothesis** (why KO reverses the pathological state): _TODO from Step 3_
- **Cell state acted on**: {cell_state}
- **Novelty / competition**: {crowding}  _(confirm via ClinicalTrials)_
- **Druggability** (ChEMBL target class, PDBs, known ligands): _TODO_
- **IPF functional / genetic evidence** (PubMed, with PMIDs): _TODO_
- **Structural tractability** (Boltz pocket / binder): _TODO_
- **Suggested modality**: _small molecule / biologic / degrader / senolytic_
- **Key risks** (essentiality, on-target tox, expression breadth): _TODO_

---
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()
    cfg = load_config(args.config)

    ranked_csv = cfg["ranking"]["out_csv"]
    if os.path.exists(ranked_csv):
        ranked = pd.read_csv(ranked_csv).head(cfg["triage"]["top_n"])
    else:
        # Bootstrap from the validated seed shortlist so the template is runnable
        # before the GPU perturbation stage has produced ranked_targets.csv.
        print(f"[info] {ranked_csv} not found; bootstrapping from seed_shortlist.csv")
        ranked = pd.read_csv("data/seed_shortlist.csv")
        ranked["gated_score"] = ""
        ranked["crowding"] = "novel"

    disease = cfg["disease"]["short"]
    genes = ranked["gene_symbol"].tolist()

    os.makedirs("results", exist_ok=True)
    worklist = build_worklist(genes, cfg["disease"]["name"])
    with open("results/triage_worklist.json", "w") as fh:
        json.dump(worklist, fh, indent=2)

    with open(cfg["triage"]["out_dossier"], "w") as fh:
        fh.write(DOSSIER_HEADER.format(disease=disease))
        for i, row in enumerate(ranked.itertuples(), 1):
            fh.write(
                DOSSIER_ENTRY.format(
                    rank=i,
                    gene=row.gene_symbol,
                    score=getattr(row, "gated_score", ""),
                    cell_state=getattr(row, "cell_state", "—"),
                    crowding=getattr(row, "crowding", "novel"),
                )
            )
    print(f"[done] worklist -> results/triage_worklist.json")
    print(f"[done] dossier skeleton -> {cfg['triage']['out_dossier']}")


if __name__ == "__main__":
    main()
