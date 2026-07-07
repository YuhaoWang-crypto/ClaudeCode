"""
Module 8 — Map the network-stability biomarker onto real trial endpoints.

Data pulled live from ClinicalTrials.gov this session:
  * 71 trials for KRAS G12C + sotorasib/adagrasib
  * pivotal CodeBreaK 100 (NCT03600883, sotorasib) and KRYSTAL/NCT04685135
    (adagrasib vs docetaxel) endpoint structure.

The point: each mathematical layer of this pipeline corresponds to a
measured clinical endpoint, so the abstract biomarkers are not free-floating
--- they map onto quantities trials already collect.
"""

# Endpoints of CodeBreaK 100 (NCT03600883), from analyze_endpoints (this session)
CODEBREAK100 = {
    "nct_id": "NCT03600883",
    "drug": "sotorasib (AMG-510)",
    "primary": [
        "Objective response rate (ORR, RECIST 1.1)",
        "Duration of response (DOR)",
        "Disease control rate",
        "Treatment-emergent / grade>=3 adverse events, DLTs",
    ],
    "secondary": [
        "Progression-free survival (PFS)",
        "Overall survival (OS)",
        "Cmax / Tmax / AUC of sotorasib (PK)",
        "sotorasib exposure vs QTc interval (cardiac safety)",
    ],
}

# How each pipeline layer maps to a real, measured endpoint.
MAPPING = [
    ("M4  DNB / LLE (network tips OFF the oncogenic ON-state)",
     "ORR + DOR (RECIST): tumour shrinks = KRAS-ON core tipped off",
     "primary efficacy"),
    ("M6  target engagement E (ChEMBL pIC50 + Boltz binding)",
     "Cmax / AUC (PK) — exposure sets fraction of node occupied",
     "secondary PK"),
    ("M9  occupancy(Cfree, IC50) -> perturbation mu",
     "Cmax/AUC + occupancy PK/PD modelling",
     "secondary PK/PD"),
    ("M4  LLE -> 0+/positive (loss of stability = toxicity hypothesis)",
     "grade>=3 AEs, DLTs, exposure-QTc relationship",
     "safety"),
    ("M1/M5 irreducible-core / master-node identity",
     "KRAS G12C mutation status = the enrolment biomarker itself",
     "patient selection"),
]


def report():
    print("=" * 68)
    print("MODULE 8 — Biomarker layers  <->  real clinical trial endpoints")
    print("=" * 68)
    print(f"pivotal trial: {CODEBREAK100['nct_id']}  "
          f"({CODEBREAK100['drug']})")
    print("  primary endpoints:")
    for e in CODEBREAK100["primary"]:
        print(f"    - {e}")
    print("  secondary endpoints:")
    for e in CODEBREAK100["secondary"]:
        print(f"    - {e}")

    print("\nmapping pipeline layer  ->  measured endpoint:")
    for layer, endpoint, kind in MAPPING:
        print(f"\n  [{kind}]")
        print(f"    {layer}")
        print(f"      -> {endpoint}")

    print("\nkey point: the two most abstract biomarkers in this pipeline")
    print("(M4 LLE/DNB for efficacy-tipping, and LLE->positive for toxicity)")
    print("map onto endpoints trials ALREADY measure (ORR/DOR and grade>=3")
    print("AEs / exposure-QTc). PK secondary endpoints (Cmax/AUC) are exactly")
    print("the inputs M9 needs to turn exposure into network perturbation mu.")
    return {"trial": CODEBREAK100, "mapping": MAPPING}


if __name__ == "__main__":
    report()
