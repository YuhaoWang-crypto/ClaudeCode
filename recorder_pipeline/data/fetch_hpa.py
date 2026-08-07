"""Fetch HPA per-gene records for every recorder marker parent/writer/carrier."""
import json, subprocess, sys, time

# marker label -> (gene symbol, Ensembl). Ensembl verified against returned "Gene".
GENES = [
    ("GFAP (TBI parent)",            "GFAP",     "ENSG00000131095"),
    ("Uromodulin (CKD reserve)",     "UMOD",     "ENSG00000169344"),
    ("NGAL (AKI)",                   "LCN2",     "ENSG00000148346"),
    ("KIM-1 (AKI)",                  "HAVCR1",   "ENSG00000113249"),
    ("Collagen III (PRO-C3/C3M)",    "COL3A1",   "ENSG00000168542"),
    ("Collagen VI a3 (PRO-C6/C6M)",  "COL6A3",   "ENSG00000163359"),
    ("Albumin (C-Alb carrier)",      "ALB",      "ENSG00000163631"),
    ("ApoA-I (Tyr192 parent)",       "APOA1",    "ENSG00000118137"),
    ("MPO (ApoA-I writer)",          "MPO",      "ENSG00000005381"),
    ("PADI4 (citrullination writer)","PADI4",    "ENSG00000159339"),
    ("Elastase (NET writer)",        "ELANE",    "ENSG00000197561"),
    ("CA125 (ovarian)",              "MUC16",    "ENSG00000181143"),
    ("uEGF (tubular mass)",          "EGF",      "ENSG00000138798"),
    ("CR1 (C4d normalizer)",         "CR1",      "ENSG00000203710"),
    ("C3 (complement)",              "C3",       "ENSG00000125730"),
    ("C4A (C4d parent)",             "C4A",      "ENSG00000244731"),
    ("Antithrombin (TAT)",           "SERPINC1", "ENSG00000117601"),
    ("Prothrombin (TAT)",            "F2",       "ENSG00000180210"),
    ("Plasminogen (PIC)",            "PLG",      "ENSG00000122194"),
    ("a2-antiplasmin (PIC)",         "SERPINF2", "ENSG00000167711"),
    ("PAI-1 (tPAIC)",                "SERPINE1", "ENSG00000106366"),
    ("Thrombomodulin",               "THBD",     "ENSG00000178726"),
    ("NT-proBNP",                    "NPPB",     "ENSG00000120937"),
    ("PSA",                          "KLK3",     "ENSG00000142515"),
    ("CEA",                          "CEACAM5",  "ENSG00000105388"),
    ("Tau (p-tau217)",               "MAPT",     "ENSG00000186868"),
    ("APP (Abeta42/40)",             "APP",      "ENSG00000142192"),
    ("a-synuclein (SAA)",            "SNCA",     "ENSG00000145335"),
    ("Transferrin (%CDT)",           "TF",       "ENSG00000091513"),
    ("AFP (AFP-L3)",                 "AFP",      "ENSG00000081051"),
    ("Haemoglobin B (HbA1c)",        "HBB",      "ENSG00000244734"),
    ("HER2 (sHER2)",                 "ERBB2",    "ENSG00000141736"),
    ("ADAM10 (HER2 sheddase)",       "ADAM10",   "ENSG00000137845"),
    ("MMP9",                         "MMP9",     "ENSG00000100985"),
    ("TIMP1",                        "TIMP1",    "ENSG00000102265"),
]

KEEP = [
    "Gene", "Ensembl",
    "RNA tissue specificity", "RNA tissue specificity score",
    "RNA tissue distribution", "RNA tissue specific nTPM",
    "RNA single cell type specificity", "RNA single cell type specificity score",
    "RNA single cell type distribution", "RNA single cell type specific nCPM",
    "RNA single cell type group specific nCPM",
    "RNA blood cell specificity", "RNA blood cell distribution",
    "RNA blood cell specific nTPM",
    "RNA tissue cell type enrichment",
    "Blood concentration - Conc. blood MS [pg/L]",
    "Blood concentration - Conc. blood IM [pg/L]",
    "Secretome location", "Protein class",
    "Single cell expression cluster",
]

out = {}
for label, sym, ens in GENES:
    url = f"https://www.proteinatlas.org/{ens}.json"
    try:
        raw = subprocess.run(["curl", "-sSL", "--compressed", "--max-time", "60", url],
                             capture_output=True, text=True, check=True).stdout
        d = json.loads(raw)
        if isinstance(d, list):
            d = d[0]
    except Exception as e:
        print(f"FAIL {label} {ens}: {e}", file=sys.stderr)
        continue
    got = d.get("Gene")
    if got != sym:
        print(f"MISMATCH {label}: asked {sym} got {got} ({ens})", file=sys.stderr)
    rec = {k: d.get(k) for k in KEEP}
    rec["_label"] = label
    rec["_verified"] = (got == sym)
    out[label] = rec
    print(f"ok {label:32s} {sym:10s} tissue={d.get('RNA tissue specificity')}"
          f" | blood={d.get('RNA blood cell specificity')}")
    time.sleep(0.25)

json.dump(out, open("hpa_markers.json", "w"), indent=1)
print(f"\nwrote {len(out)}/{len(GENES)} records")
