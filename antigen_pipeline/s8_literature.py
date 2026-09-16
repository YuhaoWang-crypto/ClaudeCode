"""Step 8 -- literature evidence for the head of the ranking.

Each of the top candidates is queried against Europe PMC with three fixed
query shapes: indication biology, surface-targeting modality, and clinical
stage.  Combined with the Open Targets drug record from step 4, that places
every candidate in one of four evidence classes, so a reader can tell an
already-clinical target from a genuinely novel one without trusting the score
alone.

Hit counts and the titles behind them are written to disk; nothing in the
report is recalled from memory.
"""

from __future__ import annotations

import time
import urllib.parse

import pandas as pd
import requests

from . import config as C

MODALITY_TERMS = ('"antibody-drug conjugate" OR "CAR-T" OR "CAR T cell" OR '
                  '"bispecific antibody" OR "monoclonal antibody" OR "immunotoxin"')
INDICATION_TERMS = '"lung adenocarcinoma" OR "non-small cell lung cancer" OR NSCLC'
# The gene symbol is pinned to title/abstract. Unrestricted full-text matching
# over-counts by one to two orders of magnitude -- a paper that merely cites a
# sequence accession is not evidence about the target.
GENE_FIELD = "TITLE_ABS"
# A symbol that cannot exist, queried the same way, calibrates the floor.
CONTROL_SYMBOL = "ZZZFAKE1"


def _search(query: str, page_size: int = 5, retries: int = 3) -> dict:
    # No sort parameter: the API's CITED_BY_COUNT sort silently returns an
    # empty body, and relevance order is what we want for the top titles.
    params = {"query": query, "format": "json", "pageSize": page_size,
              "resultType": "lite"}
    url = f"{C.EUROPEPMC_REST}?{urllib.parse.urlencode(params)}"
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=90)
            r.raise_for_status()
            payload = r.json()
            if "hitCount" not in payload:
                raise RuntimeError(f"malformed Europe PMC response: {list(payload)}")
            return payload
        except Exception as exc:
            if attempt == retries - 1:
                print(f"    europe pmc query failed: {exc}")
                return {"hitCount": None, "resultList": {"result": []}}
            time.sleep(2 ** (attempt + 1))
    return {"hitCount": None, "resultList": {"result": []}}


def _titles(payload: dict, k: int = 3) -> list[dict]:
    out = []
    for r in (payload.get("resultList", {}).get("result") or [])[:k]:
        out.append({"pmid": r.get("pmid") or r.get("id"), "year": r.get("pubYear"),
                    "journal": r.get("journalTitle"), "title": r.get("title")})
    return out


def classify(row, modality_hits, indication_hits, clinical_hits) -> tuple[str, str]:
    stage = str(row.get("ot_max_clinical_stage") or "")
    n_drugs = int(row.get("ot_n_drugs") or 0)
    if n_drugs > 0 and modality_hits and modality_hits >= 5:
        return ("clinically validated target",
                f"{n_drugs} drug programmes in Open Targets (max stage {stage or 'n/a'}), "
                f"{modality_hits} surface-modality papers")
    if modality_hits and modality_hits >= 3:
        return ("preclinical surface-targeting evidence",
                f"{modality_hits} papers pairing the gene with an antibody/CAR modality")
    if indication_hits and indication_hits >= 10:
        return ("tumour-biology evidence only",
                f"{indication_hits} lung-cancer papers but no established surface-targeting programme")
    return ("novel -- no surface-targeting literature",
            f"{indication_hits or 0} indication papers, {modality_hits or 0} modality papers")


def run(ranked: pd.DataFrame | None = None, n: int | None = None) -> pd.DataFrame:
    if ranked is None:
        ranked = pd.read_csv(C.RESULTS / "s6_ranked_candidates.csv")
    n = n or C.N_LITERATURE_CANDIDATES
    top = ranked.head(n)
    print(f"[S8] literature triage of the top {n} candidates (Europe PMC)")

    control = _search(f'({GENE_FIELD}:"{CONTROL_SYMBOL}") AND ({MODALITY_TERMS})', 1)
    print(f"  query calibration: the impossible symbol {CONTROL_SYMBOL} returns "
          f"{control.get('hitCount')} hits")

    rows = []
    for r in top.itertuples():
        gene = r.gene
        q_ind = f'({GENE_FIELD}:"{gene}") AND ({INDICATION_TERMS})'
        q_mod = f'({GENE_FIELD}:"{gene}") AND ({MODALITY_TERMS})'
        q_clin = (f'({GENE_FIELD}:"{gene}") AND ({MODALITY_TERMS}) AND '
                  f'(PUB_TYPE:"Clinical Trial" OR "phase 1" OR "phase I" OR "phase 2")')
        p_ind, p_mod, p_clin = _search(q_ind), _search(q_mod), _search(q_clin)
        h_ind = p_ind.get("hitCount")
        h_mod = p_mod.get("hitCount")
        h_clin = p_clin.get("hitCount")
        klass, note = classify(r._asdict(), h_mod, h_ind, h_clin)
        rows.append({
            "rank": r.rank, "gene": gene, "final_score": r.final_score, "tier": r.tier,
            "evidence_class": klass, "evidence_note": note,
            "hits_indication": h_ind, "hits_modality": h_mod, "hits_clinical": h_clin,
            "ot_n_drugs": getattr(r, "ot_n_drugs", None),
            "ot_drugs": getattr(r, "ot_drugs", ""),
            "top_modality_papers": "; ".join(
                f"{t['title']} ({t['journal']}, {t['year']}, PMID {t['pmid']})"
                for t in _titles(p_mod)),
            "top_indication_papers": "; ".join(
                f"{t['title']} ({t['journal']}, {t['year']}, PMID {t['pmid']})"
                for t in _titles(p_ind, 2)),
        })
        print(f"  {r.rank:>3} {gene:<10} {klass:<42} "
              f"indication={h_ind} modality={h_mod} clinical={h_clin}")

    lit = pd.DataFrame(rows)
    lit.to_csv(C.RESULTS / "s8_literature_evidence.csv", index=False)

    counts = lit.evidence_class.value_counts().to_dict()
    C.write_provenance("s8", {
        "candidates_reviewed": int(len(lit)),
        "source": C.EUROPEPMC_REST,
        "gene_field": GENE_FIELD,
        "control_symbol": CONTROL_SYMBOL,
        "control_hits": control.get("hitCount"),
        "evidence_class_counts": {str(k): int(v) for k, v in counts.items()},
        "per_gene": lit[["rank", "gene", "evidence_class", "hits_indication",
                         "hits_modality", "hits_clinical"]].to_dict("records"),
    })
    return lit


if __name__ == "__main__":
    run()
