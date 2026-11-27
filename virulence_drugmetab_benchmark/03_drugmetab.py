"""Microbial drug-biotransformation potential from genomes: survey canonical gut drug-metabolizing
enzymes across representative gut species, and test whether capability is predictable from species
(conserved) or strain-specific (variable). Activity validation needs wet-lab (Zimmermann) data - gated."""
import os, json, time, urllib.request, urllib.parse
OUT = "results"; HDR = {"User-Agent": "research-client/1.0", "Accept": "application/json"}
def api(path, q):
    u = "https://www.bv-brc.org/api/" + path + "/?" + urllib.parse.quote(q, safe="(),=&") + "&http_accept=application/json&limit(1)"
    for _ in range(4):
        try:
            r = urllib.request.urlopen(urllib.request.Request(u, headers=HDR), timeout=60)
            return int(r.headers.get("Content-Range", "items 0-0/0").split("/")[-1])
        except Exception: time.sleep(3)
    return 0
# drug-metabolizing enzymes (drug they act on)
enz = {"beta-glucuronidase": ("GUS", "SN-38/irinotecan, NSAID reactivation"),
       "azoreductase": ("azoR", "sulfasalazine -> 5-ASA"),
       "nitroreductase": ("nfs", "metronidazole, nitrazepam"),
       "choloylglycine hydrolase": ("BSH", "bile salt hydrolase; drug deconjugation"),
       "tyrosine decarboxylase": ("tyrDC", "L-DOPA inactivation"),
       "cardiac glycoside": ("cgr", "digoxin inactivation")}
# representative gut taxa (species taxon_id)
taxa = {"Bacteroides fragilis": 817, "Bacteroides thetaiotaomicron": 818, "Escherichia coli": 562,
        "Eggerthella lenta": 84112, "Enterococcus faecalis": 1351, "Enterococcus faecium": 1352,
        "Bifidobacterium longum": 216816, "Lactobacillus plantarum": 1590, "Clostridioides difficile": 1496,
        "Parabacteroides distasonis": 823, "Akkermansia muciniphila": 239935, "Clostridium perfringens": 1502}
res = {}
tot = {}
for sp, tid in taxa.items():
    tot[sp] = api("genome", f"eq(taxon_id,{tid})")
    res[sp] = {}
    for e in enz:
        n = api("genome_feature", f'and(eq(taxon_id,{tid}),keyword("{e}"))')
        res[sp][e] = n
# genomes carrying enzyme (approx by feature count presence); report prevalence = genomes_with_feature/total
S = {"enzymes": {e: {"short": enz[e][0], "drug": enz[e][1]} for e in enz}, "total_genomes": tot, "feature_counts": res,
     "prevalence_fraction": {}}
for sp in taxa:
    T = max(tot[sp], 1)
    S["prevalence_fraction"][sp] = {enz[e][0]: round(res[sp][e] / T, 2) for e in enz}
json.dump(S, open(f"{OUT}/drugmetab_summary.json", "w"), indent=1)
print("total genomes per species:", tot)
import pandas as pd
P = pd.DataFrame(S["prevalence_fraction"]).T
print(P.to_string())
