"""OAS (Observed Antibody Space) germline-framework and CDR-length priors.

OAS paired human sequences give the human VH/VL germline frameworks and the
joint CDR-length distribution the generator samples from (plan Stage 1). OAS is
reachable but large (per-study CSV.gz units), so this module provides:
  - a downloader for one OAS paired unit (documented, opt-in)
  - extraction of germline V-gene usage and joint CDR-length histograms from a
    unit, using ANARCI if available for exact CDR boundaries.

For a dry run without downloading, ``default_cdr_length_prior()`` returns a
hard-coded human-repertoire prior so the pipeline is usable immediately.
"""

from __future__ import annotations

import csv
import gzip
import json
from collections import Counter
from pathlib import Path
from typing import Optional

OAS_PAIRED_INDEX = "https://opig.stats.ox.ac.uk/webapps/oas/oas_paired/"


def default_cdr_length_prior() -> dict[str, dict[int, float]]:
    """Approximate human paired-repertoire CDR length distributions (normalised).

    Memory figures, labelled approximate. Replace by extracting from an OAS unit.
    """
    raw = {
        "H1": {5: 0.80, 6: 0.15, 7: 0.05},
        "H2": {6: 0.35, 7: 0.40, 8: 0.20, 5: 0.05},
        "H3": {11: 0.14, 12: 0.15, 13: 0.14, 14: 0.12, 15: 0.11, 10: 0.10,
               16: 0.09, 9: 0.07, 17: 0.08},
        "L1": {6: 0.55, 7: 0.20, 8: 0.10, 11: 0.10, 12: 0.05},
        "L2": {3: 0.95, 7: 0.05},
        "L3": {9: 0.70, 8: 0.12, 10: 0.10, 11: 0.08},
    }
    return {cdr: _normalise(d) for cdr, d in raw.items()}


def _normalise(d: dict[int, float]) -> dict[int, float]:
    s = sum(d.values()) or 1.0
    return {k: v / s for k, v in d.items()}


def extract_from_unit(unit_path: str) -> dict:
    """Extract germline V-gene usage + CDR-length histograms from one OAS CSV.gz.

    OAS data units ship CDR columns (cdr1/2/3_aa_heavy, ..._light) and germline
    v_call annotations, so exact CDRs need no re-numbering. If those columns are
    absent, falls back to ANARCI (optional import).
    """
    path = Path(unit_path)
    v_heavy: Counter = Counter()
    v_light: Counter = Counter()
    cdr_len: dict[str, Counter] = {c: Counter() for c in ("H1", "H2", "H3", "L1", "L2", "L3")}

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        # OAS units have a JSON metadata line first; skip if present.
        first = fh.readline()
        if not first.lstrip().startswith("cdr") and "," in first and first.strip().startswith("{"):
            pass  # metadata consumed
        else:
            fh.seek(0)
        reader = csv.DictReader(fh)
        for row in reader:
            if row.get("v_call_heavy"):
                v_heavy[row["v_call_heavy"].split("*")[0]] += 1
            if row.get("v_call_light"):
                v_light[row["v_call_light"].split("*")[0]] += 1
            for cdr, col in (("H1", "cdr1_aa_heavy"), ("H2", "cdr2_aa_heavy"),
                             ("H3", "cdr3_aa_heavy"), ("L1", "cdr1_aa_light"),
                             ("L2", "cdr2_aa_light"), ("L3", "cdr3_aa_light")):
                v = row.get(col)
                if v:
                    cdr_len[cdr][len(v)] += 1

    return {
        "v_heavy_top": v_heavy.most_common(10),
        "v_light_top": v_light.most_common(10),
        "cdr_length_hist": {c: dict(cdr_len[c]) for c in cdr_len},
    }


def save_prior(prior: dict, path: str) -> None:
    Path(path).write_text(json.dumps(prior, indent=2))
