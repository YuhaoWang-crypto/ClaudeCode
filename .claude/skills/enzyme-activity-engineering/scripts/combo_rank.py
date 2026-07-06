import json, itertools
import numpy as np

combo = json.load(open("data/combo_scores.json"))
rows = json.load(open("data/ranked_mutations.json"))
by = {r["mutation"]: r for r in rows}

E1 = combo["esm1v_t33_650M_UR90S_1"]["scores"]
E2 = combo["esm2_t33_650M_UR50D"]["scores"]

SINGLE = {168: "L168K", 193: "H193R", 231: "F231Y", 248: "A248K", 257: "V257A"}
CLASS = {168: "fold", 193: "fold", 231: "RNA", 248: "DNA", 257: "turnover"}
# single-mutant structural LLRs (no epistasis available from per-position marginals)
mpnn_s = {p: by[m]["raw_llr"]["ProteinMPNN"] for p, m in SINGLE.items()}
esmif_s = {p: by[m]["raw_llr"]["ESM-IF"] for p, m in SINGLE.items()}

recs = []
for key, d in E1.items():
    combo_pos = d["combo"]
    n = len(combo_pos)
    esm1_epi, esm1_add = E1[key]["epistatic"], E1[key]["additive"]
    esm2_epi, esm2_add = E2[key]["epistatic"], E2[key]["additive"]
    esm_joint = (esm1_epi + esm2_epi) / 2          # epistasis-aware ESM (both models)
    esm_add = (esm1_add + esm2_add) / 2
    epistasis = esm_joint - esm_add                 # + = synergy, - = interference
    struct_add = sum(mpnn_s[p] for p in combo_pos) + sum(esmif_s[p] for p in combo_pos)
    classes = sorted(set(CLASS[p] for p in combo_pos))
    # mechanism diversity: distinct activity routes covered
    recs.append({"key": key, "n": n, "pos": combo_pos,
                 "esm_joint": esm_joint, "esm_add": esm_add, "epistasis": epistasis,
                 "esm1_epi": esm1_epi, "esm2_epi": esm2_epi,
                 "struct_add": struct_add, "classes": classes, "ndiv": len(classes)})

singles = [r for r in recs if r["n"] == 1]
multis = [r for r in recs if r["n"] >= 2]
multis.sort(key=lambda r: r["esm_joint"], reverse=True)
json.dump(recs, open("data/combo_ranked.json", "w"))

def fmt(r):
    return (f"{r['key']:<28} n={r['n']} | ESM-joint {r['esm_joint']:+6.2f} | add {r['esm_add']:+6.2f} "
            f"| epi {r['epistasis']:+5.2f} | struct+ {r['struct_add']:+6.2f} | {'/'.join(r['classes'])}")

print("=== SINGLES (reference) ===")
for r in sorted(singles, key=lambda r: r["esm_joint"], reverse=True):
    print("  ", fmt(r))

print("\n=== MULTI-MUTANTS ranked by epistasis-aware ESM joint score ===")
for r in multis:
    print("  ", fmt(r))

print("\n=== Largest epistasis magnitudes (|ESM joint - additive|) ===")
for r in sorted(multis, key=lambda r: abs(r["epistasis"]), reverse=True)[:8]:
    tag = "SYNERGY" if r["epistasis"] > 0 else "INTERFERENCE"
    print(f"   {r['key']:<30} epistasis {r['epistasis']:+.2f}  [{tag}]  classes={'/'.join(r['classes'])}")

print("\n=== F231Y / V257A interaction (the only spatially-close pair, 6.5A) ===")
for r in recs:
    if set(r["pos"]) == {231, 257}:
        print(f"   F231Y+V257A: ESM-1v epi {r['esm1_epi']:+.2f}, ESM2 epi {r['esm2_epi']:+.2f}, "
              f"joint {r['esm_joint']:+.2f} vs additive {r['esm_add']:+.2f} -> epistasis {r['epistasis']:+.2f}")
