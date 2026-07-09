# Config example — retarget the whole campaign by editing this

```python
config = {
    "disease": {"name": "NSCLC", "mondo": "MONDO_0005138", "efo": "MONDO_0005233"},
    # Stage 1: indications to score against NSCLC (each scored 0-5 per dimension)
    "candidate_indications": [
        "NSCLC", "Colorectal", "HCC", "Pancreatic", "SLE", "TNBC",
        "Ulcerative colitis", "Obesity", "MASH", "RA", "Atopic dermatitis",
    ],
    # Stage 2 (optional): if the axis is known, seed it; else discover from OT top targets
    "seed_targets": ["KRAS", "KEAP1", "NFE2L2", "SLC7A11", "GLS"],
    # Stage 2b: signed regulatory edges for the network-biomarker layer
    "pathway_edges": [
        ("KEAP1","NRF2","-"), ("NRF2","SLC7A11","+"), ("NRF2","GLS","+"),
        ("NRF2","NRF2","+"), ("KRAS","NRF2","+"), ("SLC7A11","NRF2","+"),
    ],
    # cohorts for cBioPortal mutation frequencies (Stage 2)
    "cohorts": ["luad_tcga_pan_can_atlas_2018","lusc_tcga_pan_can_atlas_2018","nsclc_tcga_broad_2016"],
    "antibody_design_cap": 100,   # cost control
    "cofolds_per_job": 4,
}
```

To run on a NEW disease: change `disease`, `candidate_indications`, and either
`seed_targets`+`pathway_edges` (if known) or leave them empty to discover the axis
from Open Targets top-associated targets in Stage 2. Pick cohorts matching the new
tumor type from cBioPortal `cancer_studies`.

## The target-axis principle (Stage 2)

Pick 3-5 targets that form a **mechanistic chain across biological levels**
(genomic driver -> signaling/TF -> surface effector -> metabolic effector), not the
5 highest-scoring independent targets. The chain is what lets each target map to its
biology-appropriate modality in Stage 5, and what the network-biomarker layer (2b)
collapses to an irreducible core.
