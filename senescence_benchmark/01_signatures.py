"""Build senescence signatures (public) as gene-symbol sets."""
import os, json, csv
S = os.environ["SEN_DIR"]; Q = os.environ["QB_DATA"]
def msig(name): d = json.load(open(f"{S}/../sen_{name}.json")); return set(d[list(d)[0]]["geneSymbols"])
sig = {}
sig["SenMayo"] = msig("SAUL_SEN_MAYO")
sig["Fridman_up"] = msig("FRIDMAN_SENESCENCE_UP")
sig["Reactome_senescence"] = msig("REACTOME_CELLULAR_SENESCENCE")
# CellAge genes that induce senescence
ca = list(csv.DictReader(open(f"{S}/cellage3.tsv"), delimiter="\t"))
sig["CellAge_induces"] = {x["Gene symbol"] for x in ca if x.get("Senescence Effect") == "Induces"}
# SenCID feature set (seneset) -> symbols (already symbols)
sig["SenCID_geneset"] = {l.strip() for l in open(f"{S}/../ext/sencid/SenCID/resource/seneset.txt") if l.strip()}
# simple curated up-markers (arrest + SASP + lysosomal)
sig["Core_markers_up"] = {"CDKN1A", "CDKN2A", "SERPINE1", "GLB1", "TP53", "IL6", "CXCL8", "IL1B", "MMP3", "MMP1", "TNF", "IGFBP3", "GDF15"}
# proliferation-loss (senescent cells lose these; used as negative -> we invert at scoring)
sig["Proliferation_down"] = {"MKI67", "TOP2A", "CCNB1", "CCNB2", "CDK1", "LMNB1", "PCNA", "CENPA", "BUB1", "FOXM1"}
json.dump({k: sorted(v) for k, v in sig.items()}, open("results/signatures.json", "w"))
# cross-type annotation from CellAge (which senescence type each inducing gene is linked to)
typ = {}
for x in ca:
    if x.get("Senescence Effect") == "Induces":
        typ.setdefault(x.get("Type of senescence", "Unclear"), []).append(x["Gene symbol"])
json.dump({k: len(v) for k, v in typ.items()}, open("results/cellage_types.json", "w"))
print({k: len(v) for k, v in sig.items()})
print("CellAge senescence types:", {k: len(v) for k, v in typ.items()})
