# External resource inventory

Access status verified during Phase-0 scaffolding. "Problem" column refers to the
A/B split in [`PLAN.md`](PLAN.md).

| Resource | URL | Problem | Access | Notes |
|---|---|---|---|---|
| **LNP_ML / LiON** | github.com/jswitten/LNP_ML | A | ✅ public, **MIT** | Chemprop 1.7.0, py3.8. Ships `data/all_data.csv` (>9k pts) + `main_script.py`. No pretrained checkpoint; retrain via one command. → `scripts/fetch_lnp_ml.sh` |
| **Su 2026 spatial-conf** | doi 10.1038/s41551-026-01640-8 | A | paper only | 1,408-lipid library (14 heads × {ester,amide} × 16 tails); MD → 2D density → 28 feats → SISSO. Reproduce method only if shape features help. |
| **DrugCLIP** | github.com/bowen-gao/DrugCLIP | B | ✅ public code | NeurIPS 2023; contrastive pocket↔molecule retrieval. Uni-Mol based. Needs GPU. |
| **drug-the-whole-genome** | drug-the-whole-genome.yanyanlan.com | B | ✅ **mined via API** | Genome-wide DrugCLIP (*Science* 2025). SPA backed by FastAPI `dtwgapi.yanyanlan.com`: `GET /complexes/{uniprot}` → predicted pocket-ligand hits (drugclip_score, docking_score, nearby_residues, source); `POST /get_smiles {complex_list}` → SMILES. → `scripts/fetch_drugclip.py`. |
| **humanPPI** | prodata.swmed.edu/humanPPI | B | ✅ **mined via API** | Predicted human PPIs (Grishin/Cong lab, *Science* 2025). Flask endpoint `GET /humanPPI/data/{uniprot}?filter=<bool>` → partners (AF_Score, RF_Score, DCA_Score, subcellular locality, confidence). Bulk: `conglab.swmed.edu/humanPPI/downloads/final_predictions.tar.gz`. → `scripts/fetch_humanppi.py`. |
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
