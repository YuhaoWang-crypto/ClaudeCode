"""Known-answer validation: IFN/viral-response reQTLs enrich for interferon-response genes;
reQTL sets are stimulus-coherent (IFNg ~ IAV type-I/II IFN, distinct from Pam3). Plus figures."""
import os, json, csv, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
R = os.environ["RQ_DATA"]; Q = os.environ["QB_DATA"]; OUT = "results"
# ENSG -> symbol from HGNC
ens2sym = {}
for x in csv.DictReader(open(f"{Q}/hgnc.tsv"), delimiter="\t"):
    e = x.get("ensembl_gene_id", "")
    if e: ens2sym[e] = x["symbol"]
def syms(ensids): return {ens2sym.get(str(e).split(".")[0]) for e in ensids if isinstance(e,str)} - {None}
def gset(name):
    d = json.load(open(f"{R}/{name}.json")); return set(d[list(d)[0]]["geneSymbols"])
IFNg = gset("HALLMARK_INTERFERON_GAMMA_RESPONSE"); IFNa = gset("HALLMARK_INTERFERON_ALPHA_RESPONSE"); INFL = gset("HALLMARK_INFLAMMATORY_RESPONSE")
IFN = IFNg | IFNa
ov = json.load(open(f"{OUT}/egene_overlap.json"))

def enrich(comp, gs, label):
    g = pd.read_csv(f"{OUT}/gained_{comp}.csv"); gained = syms(g.molecular_trait_id)
    base = comp.split("_to_")[0]
    allgenes = syms(pd.read_csv(f"{R}/{base}.permuted.tsv.gz", sep="\t").molecular_trait_id)
    a = len(gained & gs); b = len(gained - gs); c = len((allgenes & gs) - gained); d = len((allgenes - gs) - gained)
    orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return {"set": label, "gained_in_set": a, "gained_total": len(gained), "pct_in_set": round(100*a/len(gained),1),
            "background_pct": round(100*len(allgenes & gs)/len(allgenes),1), "odds_ratio": round(float(orr),2), "p": float(p)}
S = {"enrichment": {}}
for comp in ["naive_to_IAV", "naive_to_R848", "naive_to_LPS", "naive_to_Pam3", "mac_naive_to_mac_IFNg"]:
    S["enrichment"][comp] = {k: enrich(comp, s, k) for k, s in [("IFN response", IFN), ("inflammatory", INFL)]}
# stimulus coherence: Jaccard of gained-eGene sets between stimuli
gsets = {c: syms(pd.read_csv(f"{OUT}/gained_naive_to_{c}.csv").molecular_trait_id) for c in ["IAV","R848","LPS","Pam3"]}
def jac(a,b): return round(len(a&b)/len(a|b),3)
S["stimulus_coherence_jaccard"] = {f"{a}~{b}": jac(gsets[a],gsets[b]) for i,a in enumerate(["IAV","R848","LPS","Pam3"]) for b in ["IAV","R848","LPS","Pam3"][i+1:]}
conf = {c: json.load(open(f"{OUT}/confirm_{c}.json")) for c in ["naive_to_IAV","mac_naive_to_mac_IFNg"]}
S["tabix_confirmation"] = conf
json.dump(S, open(f"{OUT}/validation.json", "w"), indent=1, default=str)

# ---- figures
fig, ax = plt.subplots(1, 3, figsize=(18, 4.8))
comps = ["naive->LPS","naive->R848","naive->IAV","naive->Pam3","mac_naive->mac_IFNg"]
labels = ["monocyte\n+LPS","monocyte\n+R848","monocyte\n+influenza","monocyte\n+Pam3","macrophage\n+IFNγ"]
gained = [ov["comparisons"][c]["response_gained_eGenes"] for c in comps]
pct = [ov["comparisons"][c]["pct_of_stim_eGenes_baseline_invisible"] for c in comps]
x = np.arange(len(comps))
b1 = ax[0].bar(x, pct, color="#d73027")
for i,v in enumerate(pct): ax[0].text(i, v+0.6, f"{v:.0f}%\n({gained[i]})", ha="center", fontsize=8)
ax[0].set_xticks(x, labels, fontsize=8); ax[0].set(ylabel="% of stimulated eGenes invisible at baseline", title="A. Response-eQTLs a resting-state scan misses\n(gained eGenes per stimulus)", ylim=(0,45))
# B: confirmation - beta stim vs naive for IAV
d = pd.read_csv(f"{OUT}/confirm_naive_to_IAV.csv")
ax[1].scatter(d.abs_beta_naive, d.abs_beta_stim, s=14, c=np.where(d.p_naive>0.01,"#d73027","#999"), alpha=.7)
mx = max(d.abs_beta_stim.max(), d.abs_beta_naive.max())
ax[1].plot([0,mx],[0,mx],"k:",lw=.8); ax[1].set(xlabel="|effect| at baseline (naive)", ylabel="|effect| under influenza",
    title=f"B. Confirmed reQTLs (monocyte+influenza)\nred = null at baseline (p>0.01), {(d.p_naive>0.01).mean()*100:.0f}% of gained")
# C: IFN enrichment
e = S["enrichment"]; ors = [e[c]["IFN response"]["odds_ratio"] for c in ["naive_to_IAV","naive_to_R848","naive_to_LPS","naive_to_Pam3","mac_naive_to_mac_IFNg"]]
ax[2].bar(x, ors, color=["#d73027" if o>1.5 else "#999" for o in ors]); ax[2].axhline(1,ls=":",c="k",lw=.8)
for i,o in enumerate(ors): ax[2].text(i,o+0.03,f"{o:.1f}",ha="center",fontsize=9)
ax[2].set_xticks(x, labels, fontsize=8); ax[2].set(ylabel="odds ratio: reQTL genes in IFN-response set", title="C. Known-answer: IFN/viral reQTLs enrich\nfor interferon-response genes")
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_reqtl.png", dpi=130)
print(json.dumps(S, indent=1, default=str))
