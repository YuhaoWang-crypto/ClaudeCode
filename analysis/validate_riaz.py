#!/usr/bin/env python3
"""
External validation of the CD8 memory-vs-exhaustion responder signature on the
Riaz et al. (Cell 2017) anti-PD-1 (nivolumab) melanoma cohort.

Data (referenced by name, not committed): GEO GSE91061
  - GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz  (bulk RNA-seq, Entrez IDs)
  - GSE91061_series_matrix.txt.gz                        (per-sample response/visit)

The signature learned on GSE120575 single cells is applied to bulk tumor
transcriptomes: signature = mean(z[memory genes]) - mean(z[exhaustion genes]).
"""
import os
import gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

DATA = "/home/user/ClaudeCode/data"
OUT = "/home/user/ClaudeCode/results"
FIG = os.path.join(OUT, "figures")

FPKM = os.path.join(DATA, "GSE91061_fpkm.csv.gz")
SMAT = os.path.join(DATA, "GSE91061_series_matrix.txt.gz")

# signature genes -> Entrez GeneID (hg19KnownGene rows are Entrez IDs)
MEM_ENTREZ = {  # CD8_G memory-like
    "TCF7": 6932, "IL7R": 3575, "SELL": 6402, "CCR7": 1236, "CD28": 940,
    "LEF1": 51176, "CD27": 939, "GZMK": 3003,
}
EXH_ENTREZ = {  # CD8_B exhaustion
    "PDCD1": 5133, "HAVCR2": 84868, "LAG3": 3902, "TIGIT": 201633,
    "CTLA4": 1493, "ENTPD1": 953, "CD38": 952, "TOX": 9760,
}
CD8_ENTREZ = {  # CD8 T-cell abundance / infiltration
    "CD8A": 925, "CD8B": 926, "GZMA": 3001, "NKG7": 4818, "CCL5": 6352,
}


def parse_series_matrix():
    titles, char_rows = None, []
    with gzip.open(SMAT, "rt", encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                char_rows.append(vals)
    meta = pd.DataFrame({"title": titles})
    for vals in char_rows:
        joined = " ".join(vals).lower()
        if "response:" in joined:
            meta["response"] = [v.split(":", 1)[1].strip() for v in vals]
        elif "visit" in joined or "pre or on" in joined:
            meta["visit"] = [v.split(":", 1)[1].strip() for v in vals]
    return meta.set_index("title")


def main():
    meta = parse_series_matrix()
    print("[riaz] samples:", len(meta),
          "| response counts:\n", meta["response"].value_counts())

    fpkm = pd.read_csv(FPKM, index_col=0)
    fpkm.index = fpkm.index.astype(str)
    print("[riaz] fpkm genes x samples =", fpkm.shape)

    def zscores(entrez_map):
        ids = [str(e) for e in entrez_map.values() if str(e) in fpkm.index]
        missing = [g for g, e in entrez_map.items() if str(e) not in fpkm.index]
        sub = np.log2(fpkm.loc[ids] + 1.0)
        z = sub.sub(sub.mean(1), axis=0).div(sub.std(1) + 1e-9, axis=0)
        return z.mean(0), missing  # per-sample mean z

    mem_z, miss_m = zscores(MEM_ENTREZ)
    exh_z, miss_e = zscores(EXH_ENTREZ)
    cd8_z, miss_c = zscores(CD8_ENTREZ)
    print("[riaz] missing memory:", miss_m, "| exhaustion:", miss_e, "| CD8:", miss_c)

    df = meta.copy()
    df["memory"] = mem_z
    df["exhaustion"] = exh_z
    df["CD8_abundance"] = cd8_z
    df["state_diff"] = mem_z - exh_z          # the single-cell state signature
    # Responder = PRCR (CR/PR); Non-responder = PD; SD excluded (ambiguous)
    df["group"] = df["response"].map({"PRCR": "R", "PD": "NR"})
    df.to_csv(os.path.join(OUT, "riaz_signature_scores.csv"))

    scores = ["memory", "exhaustion", "CD8_abundance", "state_diff"]
    rows = []
    for visit in ["Pre", "On"]:
        for col in scores:
            d = df[(df["visit"] == visit) & df["group"].notna()]
            r = d.loc[d["group"] == "R", col]
            n = d.loc[d["group"] == "NR", col]
            if len(r) >= 3 and len(n) >= 3:
                y = (d["group"] == "R").astype(int)
                auc = roc_auc_score(y, d[col])
                p = mannwhitneyu(r, n, alternative="two-sided").pvalue
            else:
                auc, p = np.nan, np.nan
            rows.append(dict(visit=visit, score=col, n_R=len(r), n_NR=len(n),
                             AUC=round(auc, 3), pval=round(p, 4)))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "riaz_validation_stats.csv"), index=False)
    print("[riaz] validation (Responder=CR/PR vs Non-responder=PD):\n",
          res.to_string(index=False))

    # figure: AUC of each score at Pre and On
    fig, ax = plt.subplots(figsize=(7, 4.2))
    piv = res.pivot(index="score", columns="visit", values="AUC").loc[scores]
    piv.plot(kind="bar", ax=ax, color=["#4C78A8", "#F58518"])
    ax.axhline(0.5, color="grey", ls="--", lw=0.8)
    ax.set_ylabel("ROC AUC (predict CR/PR vs PD)")
    ax.set_title("Riaz anti-PD-1 (GSE91061): bulk validation of signature components")
    ax.set_ylim(0, 0.85)
    for c in ax.containers:
        ax.bar_label(c, fmt="%.2f", fontsize=7)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, "riaz_validation.png"), dpi=130)
    plt.close()
    print("[done] riaz validation")


if __name__ == "__main__":
    main()
