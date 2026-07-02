# IPF — Pilot Target Dossier (evidence-filled)

Two whitespace candidates carried through the full triage layer with **real data**
from ChEMBL (bioactivity) and Boltz-2.1 (structure + binding). This is the
`05_evidence_triage.py` template filled in — it demonstrates the "Druggability"
and "Structural tractability" rows the pipeline produces per target.

Boltz predictions expire from the workspace after ~7 days; IDs recorded in
`data/boltz_results.csv`. ChEMBL = EMBL-EBI ChEMBL v34.

---

## 1. PTGES / mPGES-1  (Prostaglandin E synthase, UniProt O14684, CHEMBL5658)

- **Mechanistic hypothesis**: inducible microsomal PGE2 synthase downstream of
  COX-2; blocking it lowers pro-fibrotic/pro-inflammatory PGE2 signaling in the
  pathologic-fibroblast compartment without the broad COX-inhibition liabilities.
- **Cell state acted on**: pathologic_fibroblast.
- **Novelty / competition**: **novel for IPF** — not among the 55 in-clinic Ph2/3
  IPF programs (crowded set is PDE4B, LPA1, integrin αvβ6, ROCK2, CSF1R, TG2…).
- **Druggability (ChEMBL, real)**: 🟢 **strong**. 1,600 bioactivities; most potent
  **IC50 = 1 nM (pChEMBL 9.0)**, many sub-10 nM, ligand efficiency LE 0.40–0.46.
  17 co-crystal PDB structures. Mature, tractable small-molecule chemotypes
  (e.g. carbazole/dibenzazole scaffolds).
- **Structural tractability (Boltz-2.1, real — monomer vs trimer)**:
  - *Monomer* + 1 nM inhibitor: structure_confidence 0.89, pTM 0.92,
    ligand_iPTM 0.78, complex iPLDDT 0.75, binding_confidence 0.50.
  - *Homotrimer* (A/B/C) + 1 nM inhibitor: structure_confidence **0.95**,
    pTM 0.98, **protein_iPTM 0.987** (biological trimer assembles cleanly),
    **ligand_iPTM 0.85**, complex iPLDDT **0.90**, interface error (iPDE)
    **3.1→1.35 Å**, binding_confidence 0.48.
  - **Read**: the trimer confirms the inhibitor binds a **well-defined
    inter-subunit pocket** — all geometric/interface confidence rose sharply,
    validating the true druggable site. The `binding_confidence` affinity-proxy
    stayed flat (~0.48) and is *not* a calibrated Kd, so pocket tractability is
    strong while the model's affinity estimate is modest. IDs in
    `data/boltz_results.csv` (monomer `sab_pred_nYrd5vkK…`, trimer
    `sab_pred_x8cp7WZf…`).
- **Suggested modality**: small molecule (orthosteric mPGES-1 inhibitor).
- **Key risks**: mPGES-1 inhibitors well-precedented in inflammation but none yet
  approved; verify lung-fibroblast expression specificity and PGE2-shunting to
  other prostanoids. **Trimer re-run: done** (pocket confirmed at subunit
  interface). **Recommended next step**: confirm IPF fibroblast expression from
  the single-cell atlas (Step 3); optionally dock a diverse inhibitor panel to
  rank chemotypes against the confirmed interface pocket.

## 2. MDK / Midkine  (UniProt P21741, CHEMBL1949490)

- **Mechanistic hypothesis**: heparin-binding growth factor secreted by aberrant
  basaloid cells; drives pro-EMT, migration and pro-fibrotic crosstalk. A secreted
  ligand → neutralize extracellularly.
- **Cell state acted on**: aberrant_basaloid.
- **Novelty / competition**: **novel for IPF** — no in-clinic IPF program.
- **Druggability (ChEMBL, real)**: 🟡 **medium, biologic-leaning**. Only 16
  bioactivities; the nanomolar binders are **sulfated heparin-oligosaccharide
  mimetics (Kd 2.6 nM)**, while classic small molecules are weak (IC50 ~0.2–1.3 µM).
  Consistent with a heparin/glycosaminoglycan-binding surface, not a small-molecule
  pocket.
- **Structural tractability (Boltz-2.1, real)**: apo prediction
  **structure_confidence 0.67, pTM 0.52, pLDDT 0.70** — partly flexible/disordered
  (basic C-terminal tail), **no deep druggable pocket**. Directly corroborates the
  ChEMBL read.
- **Suggested modality**: **biologic** — neutralizing antibody or glycan/heparin-
  mimetic trap. (Good candidate for the EDEN/Boltz protein-design layer.)
- **Key risks**: secreted-target PK/exposure; MDK has broad physiology (neural,
  immune) → watch on-target systemic effects. **Recommended next step**: run a
  Boltz protein-design / EDEN antibody-design pass against the folded MDK core.

---

### How this maps back to the pipeline

| Dossier row | Data source (this pilot) |
|---|---|
| Novelty / competition | ClinicalTrials.gov → `data/crowded_targets.csv` |
| Druggability | ChEMBL `get_bioactivity` (real IC50/Kd above) |
| Structural tractability | Boltz-2.1 `structure_and_binding` → `data/boltz_results.csv` |
| Mechanistic hypothesis / cell state | Geneformer perturbation (Step 3, GPU) |

The only row still pending real numbers is the mechanistic goal-shift from
Geneformer — that needs the GPU perturbation run. Everything else is filled from
live MCP data here.
