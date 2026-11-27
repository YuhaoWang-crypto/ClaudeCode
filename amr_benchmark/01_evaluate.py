"""AMR genotype->phenotype concordance on Klebsiella pneumoniae (BV-BRC lab phenotypes + NCBI
AMRFinder/NDARO acquired resistance genes). Resistance determinant -> predicted R; per drug report
sensitivity, specificity, concordance, and very-major errors (missed resistance = the clinical warning)."""
import os, json, re, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
A = os.environ["AMR_DATA"]; OUT = "results"
pheno = pd.read_csv(f"{A}/pheno_sel.tsv", sep="\t")
geno_all = pd.read_csv(f"{A}/genotype.tsv", sep="\t")
queried = set(geno_all.genome_id.astype(str))                               # genomes we successfully retrieved
geno = geno_all[geno_all.source == "NDARO"].copy(); geno["gene"] = geno.gene.astype(str); geno["antibiotics"] = geno.antibiotics.fillna("")

# gene-name -> resistance class (standard acquired-determinant knowledge)
def gene_classes(g, ab):
    c = set(); gl = g.lower()
    if re.search(r"kpc|ndm|vim|imp|oxa-?48|oxa-?181|oxa-?232|ges-?5|carbapenem", gl): c |= {"carbapenem", "cephalosporin", "penicillin"}
    if re.search(r"ctx-m|^blacmy|dha|fox|act|^blaoxa-?10|ges|per|veb|^blashv|^blatem", gl): c |= {"cephalosporin", "penicillin"}   # ESBL/AmpC broad (approx)
    if re.search(r"ctx-m|cmy|dha|^blaoxa-?48", gl): c |= {"monobactam"}
    if gl.startswith(("blatem", "blashv", "blaoxa-1", "blaoxa-9", "blaz")): c |= {"penicillin"}
    if re.search(r"aac\(|aph\(|ant\(|aad|rmt|arma|str[ab]|aadA", gl): c |= {"aminoglycoside"}
    if re.search(r"qnr|aac\(6'\)-ib-cr|qep|oqx|gyra|parc", gl): c |= {"quinolone"}
    if re.search(r"\bsul\d|^sul", gl): c |= {"sulfonamide"}
    if re.search(r"dfr", gl): c |= {"trimethoprim"}
    # supplement from the antibiotics annotation
    for tok in ab.lower().replace("|", " ").split():
        for k, ws in [("carbapenem",["carbapenem"]),("cephalosporin",["cephalosporin","cefiderocol"]),("penicillin",["penicillin","beta-lactam"]),
                      ("aminoglycoside",["aminoglycoside","gentamicin","tobramycin","amikacin","kanamycin"]),("quinolone",["quinolone","fluoroquinolone"]),
                      ("sulfonamide",["sulfonamide","sulfamethoxazole"]),("trimethoprim",["trimethoprim"]),("monobactam",["aztreonam","monobactam"])]:
            if any(w in tok for w in ws): c.add(k)
    return c
g2cls = {}
for gid, sub in geno.groupby("genome_id"):
    cls = set()
    for g, ab in zip(sub.gene, sub.antibiotics): cls |= gene_classes(g, ab)
    g2cls[str(gid)] = cls
drug_class = {"meropenem": {"carbapenem"}, "imipenem": {"carbapenem"},
    "ceftazidime": {"cephalosporin"}, "ceftriaxone": {"cephalosporin"}, "cefepime": {"cephalosporin"},
    "ciprofloxacin": {"quinolone"}, "levofloxacin": {"quinolone"},
    "gentamicin": {"aminoglycoside"}, "tobramycin": {"aminoglycoside"}, "amikacin": {"aminoglycoside"},
    "trimethoprim/sulfamethoxazole": {"sulfonamide", "trimethoprim"}, "aztreonam": {"monobactam", "cephalosporin"}}
pheno = pheno[pheno.genome_id.astype(str).isin(queried)]                    # only successfully-genotyped genomes
rows = []
for drug, cls in drug_class.items():
    d = pheno[pheno.antibiotic == drug]
    yt = d.R.values; yp = np.array([int(bool(g2cls.get(str(g), set()) & cls)) for g in d.genome_id])
    if len(yt) < 50 or len(set(yt)) < 2: continue
    R, Sm = yt == 1, yt == 0
    rows.append({"drug": drug, "n": int(len(yt)), "pct_R": round(float(yt.mean()), 2),
        "sensitivity": round(float(yp[R].mean()), 3), "specificity": round(float((1 - yp[Sm]).mean()), 3),
        "PPV": round(float(yt[yp == 1].mean()), 3) if (yp == 1).any() else None,
        "concordance": round(float((yp == yt).mean()), 3), "AUROC": round(float(roc_auc_score(yt, yp)), 3),
        "very_major_error_missedR": round(float((yp[R] == 0).mean()), 3), "major_error_falseR": round(float((yp[Sm] == 1).mean()), 3)})
res = pd.DataFrame(rows); res.to_csv(f"{OUT}/amr_concordance.csv", index=False)
acq = res[res.drug.isin(["meropenem","imipenem","ceftazidime","ceftriaxone","gentamicin","tobramycin","amikacin","trimethoprim/sulfamethoxazole"])]
S = {"n_genomes_genotyped": len(queried), "per_drug": res.to_dict("records"),
     "mean_AUROC_acquired_gene_drugs": round(float(acq.AUROC.mean()), 3),
     "mean_sensitivity_acquired": round(float(acq.sensitivity.mean()), 3),
     "fluoroquinolone_AUROC": float(res[res.drug=="ciprofloxacin"].AUROC.iloc[0]) if (res.drug=="ciprofloxacin").any() else None}
json.dump(S, open(f"{OUT}/summary.json", "w"), indent=1)
print(res.to_string(index=False)); print("\nmean AUROC (acquired-gene drugs):", S["mean_AUROC_acquired_gene_drugs"])
