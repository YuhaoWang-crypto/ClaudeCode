"""
Functional annotation of a gene list with g:Profiler g:GOSt
(https://biit.cs.ut.ee/gprofiler, open API), restricted to a custom
background (the measured genes). Results are cached by query hash under
figures/_data_cache/gprofiler/.
"""
import hashlib
import json
import os
import subprocess
import time

from netjac import DATA_CACHE

API = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
DIR = os.path.join(DATA_CACHE, "gprofiler")


def enrich(genes, background=None, sources=("GO:BP", "REAC", "CORUM", "KEGG"),
           threshold=0.05, tries=4, organism="hsapiens"):
    body = {"organism": organism, "query": sorted(set(genes)),
            "sources": list(sources), "user_threshold": threshold,
            "no_evidences": True}
    if background is not None:
        body["domain_scope"] = "custom"
        body["background"] = sorted(set(background))
    s = json.dumps(body, sort_keys=True)
    key = hashlib.sha1(s.encode()).hexdigest()[:16]
    os.makedirs(DIR, exist_ok=True)
    fn = os.path.join(DIR, key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    for k in range(tries):
        r = subprocess.run(["curl", "-sf", "--max-time", "300", "-X", "POST", API,
                            "-H", "Content-Type: application/json", "-d", s],
                           capture_output=True)
        if r.returncode == 0 and r.stdout:
            res = [{"source": x["source"], "name": x["name"], "p": x["p_value"],
                    "k": x["intersection_size"], "term_size": x["term_size"]}
                   for x in json.loads(r.stdout)["result"]]
            res.sort(key=lambda x: x["p"])
            json.dump(res, open(fn, "w"))
            return res
        time.sleep(2 ** (k + 1))
    return []


def top_term(genes, background=None, prefer=("REAC", "GO:BP", "CORUM", "KEGG"),
             max_term=1500):
    """Best (smallest p) term of manageable size, as 'name (source, p)'."""
    res = [x for x in enrich(genes, background) if x["term_size"] <= max_term]
    if not res:
        return "n.s.", 1.0
    x = res[0]
    return f"{x['name']} [{x['source']}]", x["p"]
