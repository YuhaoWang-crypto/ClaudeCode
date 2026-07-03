# External resource inventory

Access status verified during Phase-0 scaffolding. "Problem" column refers to the
A/B split in [`PLAN.md`](PLAN.md).

| Resource | URL | Problem | Access | Notes |
|---|---|---|---|---|
| **LNP_ML / LiON** | github.com/jswitten/LNP_ML | A | ✅ public, **MIT** | Chemprop 1.7.0, py3.8. Ships `data/all_data.csv` (>9k pts) + `main_script.py`. No pretrained checkpoint; retrain via one command. → `scripts/fetch_lnp_ml.sh` |
| **Su 2026 spatial-conf** | doi 10.1038/s41551-026-01640-8 | A | paper only | 1,408-lipid library (14 heads × {ester,amide} × 16 tails); MD → 2D density → 28 feats → SISSO. Reproduce method only if shape features help. |
| **DrugCLIP** | github.com/bowen-gao/DrugCLIP | B | ✅ public code | NeurIPS 2023; contrastive pocket↔molecule retrieval. Uni-Mol based. Needs GPU. |
| **drug-the-whole-genome** | drug-the-whole-genome.yanyanlan.com | B | web UI | Genome-wide DrugCLIP (*Science* 2025, "Deep contrastive learning enables genome-wide virtual screening"). #1 = PPI-ish, #4 = ligands-for-proteins. Check for bulk export / API; the web tool screens a query protein against libraries. |
| **humanPPI** | prodata.swmed.edu/humanPPI | B | ⚠️ 503 at check | Human protein–protein interactions (Qian/Grishin lab). For peptide/PPI-derived targeting. Retry later. |
| **ChEMBL** | (MCP tool in this session) | B | ✅ live now | Known actives / IC50-Ki per target, decoys. Use `mcp__ChEMBL__*`. |
| **Boltz** | (MCP tool in this session) | B | ✅ live now | Structure + **binding-affinity** co-folding. Validate ligand↔receptor & peptide↔receptor without local GPU. Use `mcp__Boltz_API__*`. |
| **PubChem PUG-REST** | pubchem.ncbi.nlm.nih.gov/rest/pug | A/B | ✅ live now | Name→SMILES resolution for reference structures. → `scripts/fetch_reference_lipids.py` |
| **Delivery-kinetics models** | Mihaila 2017/2019; Müller 2024 | A | papers | ODE/stochastic uptake→escape→expression. Phase-4 mechanistic simulator. |
| **Modal** | modal.com | infra | account needed | Run Chemprop/DrugCLIP/Boltz-batch/DiffDock/MD. Phase-5 wrappers. |

## Immediate to-dos to confirm access
- [ ] Clone LNP_ML, inspect `all_data.csv` schema (columns for SMILES, delivery, formulation, target).
- [ ] Check `drug-the-whole-genome` for a documented API / bulk download vs. UI-only.
- [ ] Retry humanPPI; determine download format.
- [ ] Confirm whether heavy compute runs on user GPU or Modal (drives Phase-5).
