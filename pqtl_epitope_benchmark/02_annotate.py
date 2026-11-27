"""Annotate each cis-pQTL lead with VEP (GRCh38): most-severe consequence + affected genes.
Flag = candidate epitope artifact if the lead is a protein-altering variant in the COGNATE gene."""
import os, json, time, urllib.request, csv
OUT = "results"
rows = list(csv.DictReader(open(f"{OUT}/cis_leads.tsv"), delimiter="\t"))
rows = [r for r in rows if r["lead_rsid"].startswith("rs")]
PAV = {"missense_variant", "stop_gained", "stop_lost", "start_lost", "inframe_insertion",
       "inframe_deletion", "splice_acceptor_variant", "splice_donor_variant", "protein_altering_variant",
       "frameshift_variant"}
def vep(ids):
    req = urllib.request.Request("https://rest.ensembl.org/vep/human/id",
        data=json.dumps({"ids": ids}).encode(), headers={"Content-Type": "application/json", "Accept": "application/json"})
    for _ in range(4):
        try: return json.load(urllib.request.urlopen(req, timeout=60))
        except Exception: time.sleep(3)
    return []
byrs = {}
ids = list({r["lead_rsid"] for r in rows})
for i in range(0, len(ids), 180):
    for res in vep(ids[i:i+180]):
        byrs[res["input"]] = res
    time.sleep(0.3)
out = []
for r in rows:
    v = byrs.get(r["lead_rsid"])
    if not v: continue
    msc = v.get("most_severe_consequence", "")
    tc = v.get("transcript_consequences", [])
    cognate_cons = set()
    for c in tc:
        if c.get("gene_symbol") == r["gene"]: cognate_cons |= set(c.get("consequence_terms", []))
    is_pav_cognate = bool(cognate_cons & PAV)
    r.update({"most_severe": msc, "cognate_consequences": ";".join(sorted(cognate_cons)) or "none",
              "pav_in_cognate": int(is_pav_cognate), "pav_any": int(msc in PAV)})
    out.append(r)
cols = list(out[0].keys())
with open(f"{OUT}/cis_leads_annotated.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in out: f.write("\t".join(str(r[c]) for c in cols) + "\n")
print("annotated", len(out), "| PAV in cognate gene:", sum(r["pav_in_cognate"] for r in out))
