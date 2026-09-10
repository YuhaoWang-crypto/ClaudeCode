"""Fetch activity data from ChEMBL for independent re-evaluation of the VC suite models."""
import json, urllib.request, urllib.parse, time, os, sys

B = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
OUT = os.path.dirname(os.path.abspath(__file__)) + "/data"
os.makedirs(OUT, exist_ok=True)

FIELDS = ["molecule_chembl_id", "canonical_smiles", "standard_type", "standard_relation",
          "standard_value", "standard_units", "assay_chembl_id", "assay_description",
          "document_chembl_id", "document_year", "target_chembl_id", "target_organism",
          "assay_type", "pchembl_value", "data_validity_comment"]


def get(url, tries=8):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.load(r)
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(2 ** i)


def fetch(name, **q):
    path = f"{OUT}/{name}.jsonl"
    if os.path.exists(path) and os.path.getsize(path) > 0:
        print(f"[skip] {name} already present")
        return
    q.setdefault("limit", 1000)
    rows, offset = [], 0
    while True:
        qq = dict(q); qq["offset"] = offset
        d = get(B + "?" + urllib.parse.urlencode(qq))
        acts = d["activities"]
        if not acts:
            break
        for a in acts:
            rows.append({k: a.get(k) for k in FIELDS})
        total = d["page_meta"]["total_count"]
        offset += len(acts)
        print(f"  {name}: {offset}/{total}", flush=True)
        if offset >= total:
            break
    with open(path, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    print(f"[done] {name}: {len(rows)} rows -> {path}")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    jobs = {
        "na1": dict(name="na_CHEMBL6135", target_chembl_id="CHEMBL6135"),
        "na2": dict(name="na_CHEMBL2051", target_chembl_id="CHEMBL2051"),
        "pa":  dict(name="pa_CHEMBL1169598", target_chembl_id="CHEMBL1169598"),
        "hinf": dict(name="mic_h_influenzae", target_organism="Haemophilus influenzae", standard_type="MIC"),
        "eclo": dict(name="mic_e_cloacae", target_organism="Enterobacter cloacae", standard_type="MIC"),
        "herg": dict(name="herg_CHEMBL240", target_chembl_id="CHEMBL240", standard_type="IC50"),
    }
    for k, j in jobs.items():
        if which in ("all", k):
            fetch(**j)
