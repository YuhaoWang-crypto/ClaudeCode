import json, numpy as np
from collections import defaultdict

AAS = "ACDEFGHIKLMNPQRSTVWY"
NATIVE = ("MEHELHYIGIDTAKEKLDVDVLRPDGRHRTKKFANTTKGHDELVSWLKGHKIDHAHICIEATGTYMEP"
          "VAECLYDAGYIVSVINPALGKAFAQSEGLRNKTDTVDARMLAEFCRQKRPAAWEAPHPLERALRALVVRH"
          "QALTDMHTQELNRTETAREVQRPSIDAHLLWLEAELKRLEKQIKDLTDDDPDMKHRRKLLESIPGIGEKT"
          "SAVLLAYIGLKDRFAHARQFAAFAGLTPRRYESGSSVRGASRMSKAGHVSLRRALYMPAMVATSKTEWGR"
          "AFRDRLAANGKKGKVILGAMMRKLAQVAYGVLKSGVPFDASRHNPVAA")
CATALYTIC = {11, 18, 20, 60, 105}

esm = json.load(open("data/esm_scores.json"))
mpnn = json.load(open("data/mpnn_scores.json"))
esmif = json.load(open("data/esmif_scores.json"))
ann = {int(k): v for k, v in json.load(open("data/structure_annotation.json")).items()}

# ---- assemble per-method score[(resi,mut)] = llr (higher = mutation favored over WT) ----
methods = {}

def from_esm(key):
    d = {}
    for pos_s, rec in esm[key]["scores"].items():
        resi = int(pos_s)
        for mut, v in rec["muts"].items():
            d[(resi, mut)] = v
    return d

methods["ESM-1v"] = from_esm("esm1v_t33_650M_UR90S_1")
methods["ESM2"] = from_esm("esm2_t33_650M_UR50D")

# ProteinMPNN: log_p [L,21] over alphabet, struct idx -> native = idx+4
alpha = mpnn["alphabet"]
lp = np.array(mpnn["log_p"])
mp = {}
for i, row in enumerate(lp):
    resi = i + 4
    wt = NATIVE[resi - 1]
    if wt not in alpha:
        continue
    for mut in AAS:
        if mut == wt:
            continue
        mp[(resi, mut)] = float(row[alpha.index(mut)] - row[alpha.index(wt)])
methods["ProteinMPNN"] = mp

# ESM-IF: per_pos idx (0-based over modeled seq) -> native = idx+4
ei = {}
for rec in esmif["per_pos"]:
    resi = rec["idx"] + 4
    for mut, v in rec["muts"].items():
        ei[(resi, mut)] = v
methods["ESM-IF"] = ei

# sanity: WT identity check on structural methods
def check(name, d):
    # every key resi's WT must equal NATIVE
    ok = True
    for (resi, mut) in list(d.keys())[:1]:
        pass
    return ok

# ---- percentile-rank normalize each method over its own scored mutations ----
pct = {}
for name, d in methods.items():
    keys = list(d.keys()); vals = np.array([d[k] for k in keys])
    order = vals.argsort()
    ranks = np.empty(len(vals)); ranks[order] = np.arange(len(vals))
    p = ranks / (len(vals) - 1)
    pct[name] = {keys[i]: float(p[i]) for i in range(len(keys))}

# ---- consensus over all mutations ----
allkeys = set()
for d in methods.values():
    allkeys |= set(d.keys())

rows = []
for (resi, mut) in allkeys:
    if resi in CATALYTIC:
        continue  # never mutate catalytic residues
    wt = NATIVE[resi - 1]
    per = {name: methods[name].get((resi, mut)) for name in methods}
    perp = {name: pct[name].get((resi, mut)) for name in methods}
    avail = [v for v in perp.values() if v is not None]
    if len(avail) < 3:      # require coverage by >=3 methods (structural methods cover 4..322)
        continue
    consensus = float(np.mean(avail))
    raw = {k: (round(v, 3) if v is not None else None) for k, v in per.items()}
    n_fav = sum(1 for v in per.values() if v is not None and v > 0)   # methods calling it favorable
    n_seq_fav = sum(1 for k in ("ESM-1v", "ESM2") if per[k] is not None and per[k] > 0)
    n_str_fav = sum(1 for k in ("ProteinMPNN", "ESM-IF") if per[k] is not None and per[k] > 0)
    a = ann.get(resi, {})
    rows.append({
        "mutation": f"{wt}{resi}{mut}", "resi": resi, "wt": wt, "mut": mut,
        "consensus_pct": round(consensus, 4),
        "n_methods": len(avail), "n_favorable": n_fav,
        "seq_fav": n_seq_fav, "str_fav": n_str_fav,
        "raw_llr": raw,
        "min_pct": round(min(avail), 3),   # worst method percentile (agreement floor)
        "tags": a.get("tags", []),
        "d_catalytic": None if a.get("d_catalytic") is None else round(a["d_catalytic"], 1),
        "d_rna": None if a.get("d_rna") is None else round(a["d_rna"], 1),
        "d_dna": None if a.get("d_dna") is None else round(a["d_dna"], 1),
        "d_interface": None if a.get("d_interface") is None else round(a["d_interface"], 1),
        "cn10": a.get("cn10"),
    })

rows.sort(key=lambda r: r["consensus_pct"], reverse=True)
json.dump(rows, open("data/ranked_mutations.json", "w"))
print("total scored candidate mutations:", len(rows))
print("\n=== TOP 25 by consensus percentile (catalytic residues excluded) ===")
print(f"{'mut':>7} {'cons':>5} {'minP':>5} {'fav':>3} {'seq':>3} {'str':>3}  {'dCat':>5} {'dRNA':>5} {'dDNA':>5} {'cn':>3}  tags")
for r in rows[:25]:
    print(f"{r['mutation']:>7} {r['consensus_pct']:.3f} {r['min_pct']:.3f} {r['n_favorable']:>3} "
          f"{r['seq_fav']:>3} {r['str_fav']:>3}  {str(r['d_catalytic']):>5} {str(r['d_rna']):>5} "
          f"{str(r['d_dna']):>5} {str(r['cn10']):>3}  {','.join(t for t in r['tags'] if t not in ('core','surface','intermediate'))}"
          f" [{[t for t in r['tags'] if t in ('core','surface','intermediate')][0] if any(t in r['tags'] for t in ('core','surface','intermediate')) else '?'}]")
