"""Demo: fuse perturbomics + network-biomarker + drug-MCP into lead ranking.

Offline & deterministic. Shows the four-axis funnel re-ranking single-axis
reversers. The biology is planted so the winner is checkable:

  * disease reversed by several drugs / KOs (transcriptomic axis)
  * one target (EGFR) sits at the pathway's irreducible-core + switch node
    (network axis)                       -> stand-in for network-biomarker M1/M19
  * drug evidence (ChEMBL pIC50, Boltz binding, ADME, trial phase) per compound
    -> stand-in for the ChEMBL/Boltz/ClinicalTrials MCP servers

Expected: a drug with only middling connectivity but a core+switch target that
is potent, structurally bound, and clinical-phase rises to the TOP once all
four axes are fused — which no single axis would have surfaced.

Run:  python3 -m perturbomics.demo_integrate
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .integrate import (NetworkContext, DrugEvidence, integrated_leads,
                        monitoring_biomarker, DEFAULT_WEIGHTS)

GENES = [f"G{i}" for i in range(60)]
RNG = np.random.default_rng(3)


def _sig(hits, scale=0.3):
    import pandas as pd
    s = pd.Series(0.0, index=GENES)
    for g, v in hits.items():
        s[g] = v
    return s + RNG.normal(0, scale, len(GENES))


def main() -> pd.DataFrame:
    from .signature import Signature
    from .connectivity import weighted_connectivity_score

    print("=" * 70)
    print("Integrated funnel: transcriptome x network x structure x clinic")
    print("=" * 70)

    disease = Signature(_sig({**{g: 3.0 for g in GENES[0:12]},
                              **{g: -3.0 for g in GENES[12:24]}}),
                        name="disease", modality="disease")

    # candidate perturbations with the GENE TARGET each engages
    lib = [
        Signature(_sig({g: -3.0 for g in GENES[0:12]}), name="gefitinib_like",
                  modality="drug", meta={"target": "EGFR"}),      # partial reverser
        Signature(_sig({**{g: -3.0 for g in GENES[0:12]},
                        **{g: 3.0 for g in GENES[12:24]}}), name="broad_tox",
                  modality="drug", meta={"target": "TUBB"}),      # strong reverser, bad target
        Signature(_sig({g: -3.0 for g in GENES[0:8]}), name="KO_MYC",
                  modality="crispr_ko", meta={"target": "MYC"}),
        Signature(_sig({g: 3.0 for g in GENES[16:24]}), name="KO_PTEN",
                  modality="crispr_ko", meta={"target": "PTEN"}),
    ]
    rev = pd.DataFrame([{
        "name": s.name, "modality": s.modality,
        "connectivity": weighted_connectivity_score(disease, s, k=20),
        "target": s.meta.get("target", ""), "moa": ""} for s in lib]
    ).sort_values("connectivity").reset_index(drop=True)

    print("\n[axis 1 only] transcriptomic reversers (perturbomics):")
    print(rev[["name", "target", "connectivity"]].to_string(
        index=False, formatters={"connectivity": "{:+.3f}".format}))

    # axis 2: network-biomarker says EGFR is core+switch, MYC is switch
    net = NetworkContext(core_nodes={"EGFR", "RAF"}, switch_nodes={"EGFR", "MYC"},
                         biomarker=("ERK_ppERK",))
    # axes 3-4: drug MCP evidence (illustrative stand-ins for ChEMBL/Boltz/trials)
    evidence = {
        "gefitinib_like": DrugEvidence("gefitinib_like", "EGFR", pIC50=8.5,
                                       boltz_bind=0.91, admet_ok=True, clinical_phase=4),
        "broad_tox": DrugEvidence("broad_tox", "TUBB", pIC50=7.0,
                                  boltz_bind=0.6, admet_ok=False, clinical_phase=1),
        "KO_MYC": DrugEvidence("KO_MYC", "MYC"),   # genetic: no drug evidence
        "KO_PTEN": DrugEvidence("KO_PTEN", "PTEN"),
    }

    leads = integrated_leads(rev, network=net, evidence=evidence)
    print("\n[all 4 axes] integrated leads (weights "
          f"{DEFAULT_WEIGHTS}):")
    print(leads.to_string(index=False, na_rep="—"))

    print(f"\nMonitor response with early-warning biomarker: "
          f"{monitoring_biomarker(net)}  (network-biomarker M4)")
    top = leads.iloc[0]
    print(f"\nWinner: {top['name']} (target {top['target']}, role "
          f"{top['net_role']}) — middling connectivity, but a core+switch node "
          "that is potent, bound, ADME-clean and phase-4. No single axis ranks "
          "it first; the fusion does.")
    print("\nRIGOR: axis values & the weighted sum are ✅ deterministic; "
          "'control node => better target', the E->effect map, and any "
          "therapeutic claim are ⚠️ hypotheses.")
    return leads


if __name__ == "__main__":
    main()
