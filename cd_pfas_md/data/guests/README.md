# Guest structures

By default the pipeline builds guests from the SMILES in `config/system.yaml`
(RDKit embed → antechamber). Drop a pre-built `*.mol2` here and point the guest's
`mol2:` field at it if you'd rather supply your own geometry/charges.

CRITICAL — protonation / charge states:
- The **sulfonated dye** and **PFAS** head groups must be the **deprotonated
  anion** (net charge −1) at assay pH. The SMILES in `system.yaml` already are.
- Replace the placeholder dye SMILES with **your actual reporter dye** and set
  its measured `Ka_M_inv` (from ITC or fluorescence titration) so the APR run
  has a calibration anchor.
- PFAS fluorine parameters are the main accuracy risk — see the `fluorine_policy`
  handling in `src/parameterize.py`.
