"""Score every IEDB peptide-allele pair with MHCflurry 2.0 (presentation, affinity, processing),
and build length/allele-matched decoy peptides from the human proteome for the presentation test.

MHCflurry predicts binding/presentation; it is NOT trained on T-cell immunogenicity labels, so the
IEDB Positive/Negative outcome is an independent test of the 'binding != immunogenicity' boundary.
"""
import os, numpy as np, pandas as pd
from mhcflurry import Class1PresentationPredictor
WORK = os.environ["IB_WORK"]; OUT = "results"; rng = np.random.default_rng(0)
pred = Class1PresentationPredictor.load()
supported = set(pred.supported_alleles)

pairs = pd.read_parquet(f"{WORK}/tcell_pairs.parquet")
pairs = pairs[pairs.allele.isin(supported)].reset_index(drop=True)

def score(peptides, allele):
    uniq = list(dict.fromkeys(peptides))
    r = pred.predict(peptides=uniq, alleles=[allele], verbose=0)
    m = r.drop_duplicates("peptide").set_index("peptide")[["presentation_score", "affinity", "processing_score"]]
    return m.loc[list(peptides)].values

# ---- score real pairs, grouped by allele
sc = np.full((len(pairs), 3), np.nan)
for al, idx in pairs.groupby("allele").groups.items():
    idx = np.array(idx); sc[idx] = score(pairs.seq.values[idx], al)
pairs[["presentation_score", "affinity", "processing_score"]] = sc
pairs["pred_binder_500nM"] = pairs.affinity < 500
pairs.to_parquet(f"{WORK}/pairs_scored.parquet")

# ---- decoys: for each positive peptide, draw proteome peptides of the same length for its allele
prot = []
cur = ""
for l in open(f"{WORK}/human_sp.fasta"):
    if l.startswith(">"):
        if cur: prot.append(cur)
        cur = ""
    else: cur += l.strip()
if cur: prot.append(cur)
prot = [p for p in prot if len(p) >= 12 and set(p) <= set("ACDEFGHIKLMNPQRSTVWY")]
def random_peptides(n, L):
    out = []
    while len(out) < n:
        p = prot[rng.integers(len(prot))]; i = rng.integers(0, len(p) - L + 1); out.append(p[i:i + L])
    return out
pos = pairs[pairs.pos == 1]
dec_rows = []
for (al, L), g in pos.groupby(["allele", "len"]):
    peps = random_peptides(len(g) * 3, L)                        # 3 decoys per positive
    s = score(peps, al)
    dec_rows.append(pd.DataFrame({"allele": al, "len": L, "seq": peps,
        "presentation_score": s[:, 0], "affinity": s[:, 1], "processing_score": s[:, 2]}))
dec = pd.concat(dec_rows, ignore_index=True); dec.to_parquet(f"{WORK}/decoys_scored.parquet")
print("pairs scored:", len(pairs), "| alleles:", pairs.allele.nunique(), "| decoys:", len(dec))
print("pos median presentation:", round(pos.presentation_score.median(), 3),
      "| decoy median:", round(dec.presentation_score.median(), 3))
