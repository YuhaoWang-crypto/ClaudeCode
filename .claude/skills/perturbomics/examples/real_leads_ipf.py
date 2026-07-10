"""REAL integrated leads for IPF: perturbomics WTCS x (illustrative) network
x LIVE ChEMBL potency/ADMET x LIVE ClinicalTrials phase.

Every drug axis value below was fetched live from the ChEMBL and ClinicalTrials
MCP servers this session (provenance in comments). Fed through the real
perturbomics.integrate code.
"""
import sys, json
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__file__), "..", "assets"))
import pandas as pd
from perturbomics import (NetworkContext, DrugEvidence, integrated_leads,
                          monitoring_biomarker)

# --- perturbomics transcriptomic reversers (real IPF run, |WTCS|) -------------
#     targets are the LIVE ChEMBL mechanism targets fetched this session.
reversers = pd.DataFrame([
    {"name": "procaterol",     "modality": "drug",      "connectivity": -1.000, "target": "ADRB2"},
    {"name": "trichostatin a", "modality": "drug",      "connectivity": -0.996, "target": "HDAC1"},
    {"name": "mln4924",        "modality": "drug",      "connectivity": -0.991, "target": "NAE1"},   # pevonedistat
    {"name": "canertinib",     "modality": "drug",      "connectivity": -0.950, "target": "EGFR"},
    {"name": "KO_CDK13",       "modality": "crispr_ko", "connectivity": -0.674, "target": "CDK13"},
    {"name": "KO_NBAS",        "modality": "crispr_ko", "connectivity": -0.590, "target": "NBAS"},
])

# --- axis 2: NETWORK CONTROL (ILLUSTRATIVE IPF pathway) -----------------------
# NOTE: ⚠️ illustrative core/switch set from IPF pathway biology (TGF-β/SMAD +
# RTK/ERBB→ERK), NOT a computed M1 quotient this session. To make it rigorous,
# run network-biomarker m1_symmetry/m19 on an IPF GRN and use its report() sets.
net = NetworkContext(
    core_nodes={"EGFR", "ERBB2", "SMAD3", "RAF1"},
    switch_nodes={"MAPK1", "MYC", "SNAI2"},
    biomarker=("ppERK",))

# --- axes 3 & 4: LIVE ChEMBL potency/ADMET + LIVE ClinicalTrials phase ---------
# pIC50 = best pChEMBL fetched; admet_ok from ChEMBL Ro5 (0 violations => True);
# clinical_phase = ChEMBL max_phase (ClinicalTrials trial counts noted).
# boltz_bind left None: the Boltz structure+binding job (async, minutes) was NOT
# run this session — engagement therefore uses ChEMBL potency alone.
evidence = {
    "canertinib":     DrugEvidence("canertinib", "EGFR",  pIC50=8.82, boltz_bind=None, admet_ok=True,  clinical_phase=3),  # EGFR 1.5nM; 3 CT.gov trials
    "mln4924":        DrugEvidence("mln4924",    "NAE1",  pIC50=8.33, boltz_bind=None, admet_ok=True,  clinical_phase=3),  # NEDD8 4.7nM; 41 CT.gov trials
    "trichostatin a": DrugEvidence("trichostatin a", "HDAC1", pIC50=8.84, boltz_bind=None, admet_ok=True, clinical_phase=1),  # HDAC 1.4nM
    "procaterol":     DrugEvidence("procaterol", "ADRB2", pIC50=None, boltz_bind=None, admet_ok=True,  clinical_phase=3),  # agonist; EC50 not fetched
    # CRISPR KO partners: genetic, no compound -> score on transcriptomic+network only
}

leads = integrated_leads(reversers, network=net, evidence=evidence)
pd.set_option("display.width", 160)
print("REAL integrated leads — IPF (perturbomics x network x ChEMBL x trials)\n")
print(leads.to_string(index=False, na_rep="—"))
print(f"\nMonitor a responder via early-warning biomarker: {monitoring_biomarker(net)} "
      "(network-biomarker M4 DNB) ⚠️ illustrative")
print("\nProvenance: transcriptomic=perturbomics real IPF run; engagement/clinical"
      "=LIVE ChEMBL+ClinicalTrials this session; network=⚠️ illustrative IPF pathway;"
      " Boltz binding NOT run (async).")
leads.to_json(
    "real_leads.json",
    orient="records", indent=2)
