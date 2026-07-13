"""
reactome.py — pathway enrichment / lookup via the Reactome AnalysisService and
ContentService (https://reactome.org/).

Reactome is an expert-curated pathway database. Its AnalysisService takes a gene
list and returns over-represented pathways with FDR, projected onto human.

IMPORTANT — environment note
----------------------------
reactome.org sits behind a Cloudflare bot-challenge that blocks server-side
(datacenter) requests: from this cloud container the endpoints return an
interstitial "Just a moment..." page instead of JSON. The functions below detect
that and raise ``ReactomeBlocked`` with a clear message. They work unchanged from
an un-challenged network (e.g. your laptop).

Good news: you usually don't need this. **STRING's enrichment API already
returns Reactome pathways** (category ``"Reactome Pathways"``) — see
``string_api.top_pathways`` — so Reactome coverage is available through a route
that is not Cloudflare-blocked. Use this module when you specifically need
Reactome's own analysis token, pathway hierarchy, or diagram links.
"""
from __future__ import annotations

from typing import Iterable

import requests

ANALYSIS = "https://reactome.org/AnalysisService"
CONTENT = "https://reactome.org/ContentService"
_UA = {"User-Agent": "Mozilla/5.0 (ppi_pathway skill)"}


class ReactomeBlocked(RuntimeError):
    """Raised when Cloudflare returns a challenge page instead of data."""


def _looks_blocked(resp: requests.Response) -> bool:
    ct = resp.headers.get("Content-Type", "")
    if "text/html" in ct and ("Just a moment" in resp.text
                              or "challenge" in resp.text.lower()):
        return True
    return False


def _check(resp: requests.Response):
    if _looks_blocked(resp):
        raise ReactomeBlocked(
            "reactome.org returned a Cloudflare challenge (blocked for "
            "server-side requests from this environment). Use "
            "string_api.top_pathways() for Reactome-category enrichment instead, "
            "or run this module from an un-challenged network.")
    resp.raise_for_status()


def analyse(genes: Iterable[str], projection: bool = True,
            interactors: bool = False) -> dict:
    """Submit a gene list to the Reactome AnalysisService. Returns the parsed
    JSON summary (includes a ``summary.token`` you can reuse for details/plots
    and a ``pathways`` list with p-values / FDR). Raises ``ReactomeBlocked`` if
    Cloudflare intercepts the request."""
    ids = "\n".join(str(g).strip() for g in genes if str(g).strip())
    endpoint = "identifiers/projection" if projection else "identifiers/"
    params = {"interactors": str(interactors).lower(), "pageSize": 100, "page": 1,
              "sortBy": "ENTITIES_PVALUE", "order": "ASC"}
    resp = requests.post(f"{ANALYSIS}/{endpoint}", params=params,
                         data=ids, headers={"Content-Type": "text/plain", **_UA},
                         timeout=90)
    _check(resp)
    return resp.json()


def top_pathways(genes: Iterable[str], fdr_max: float = 0.05,
                 n: int = 25) -> list[dict]:
    """Convenience: enriched Reactome pathways as flat rows (name, stId, p, fdr,
    entities found/total), FDR-filtered and sorted."""
    data = analyse(genes)
    rows = []
    for p in data.get("pathways", []):
        ent = p.get("entities", {})
        fdr = ent.get("fdr", 1.0)
        if fdr <= fdr_max:
            rows.append({
                "stId": p.get("stId"),
                "name": p.get("name"),
                "p_value": ent.get("pValue"),
                "fdr": fdr,
                "found": ent.get("found"),
                "total": ent.get("total"),
                "url": f"https://reactome.org/PathwayBrowser/#/{p.get('stId')}",
            })
    rows.sort(key=lambda r: r["fdr"])
    return rows[:n]


def version() -> str:
    resp = requests.get(f"{CONTENT}/data/database/version", headers=_UA, timeout=30)
    _check(resp)
    return resp.text.strip()


if __name__ == "__main__":
    demo = ["TP53", "BRCA1", "EGFR", "MYC", "CDK2"]
    try:
        for r in top_pathways(demo, n=8):
            print(f"{r['fdr']:.2e}  {r['found']}/{r['total']}  {r['name']}")
    except ReactomeBlocked as e:
        print("BLOCKED:", e)
