# allosteric-biosensor (self-contained skill)

Design single-component allosteric protein-switch biosensors by circularly
permuting a small ML/de-novo ligand-binding domain and grafting it into a
reporter enzyme (TEM-1 β-lactamase), then validate in silico with Boltz-2.1.
Reproduces Guo, Baker & Alexandrov, Nat. Biotechnol. 2026
(doi:10.1038/s41587-026-03081-9).

## Layout
- `SKILL.md` — the skill instructions (workflow, honesty labeling).
- `reference/` — methodology, Boltz validation, adding a system.
- `assets/system_template.py` — skeleton for a new receptor/reporter/analyte.
- `scripts/biosensor_pipeline/` — the working Python package (runnable).

## Quickstart
```bash
pip install -r requirements.txt
cd scripts
python3 -m biosensor_pipeline.test_design          # correctness + reproducibility
python3 -m biosensor_pipeline.run_repro            # build+verify libraries
python3 -m biosensor_pipeline.discover "<SMILES>" --name myanalyte   # mine a receptor
```
Boltz structure+binding runs go through the Boltz MCP tools (see
`reference/boltz-validation.md`); the deterministic design + mining half runs
fully offline.

## Honesty labeling
Every output is ✅ rigorous (deterministic construction, structure-derived
annotations, real geometric measurements, model-returned metrics) or
⚠️ hypothesis (the switch proxy, "active site intact ⇒ ON"). A wet-lab
kobs(+L)/kobs(−L) titration is the only ground-truth dynamic range.
