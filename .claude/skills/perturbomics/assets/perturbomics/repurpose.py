"""Drug repurposing — reuse an EXISTING drug for a NEW indication.

Discovering a target asks "what gene reverses the disease?". Repurposing asks a
tighter question: "which ALREADY-APPROVED drug reverses the disease **and is not
already used for it**?" — the Connectivity Map's founding use case.

Two ingredients on top of the ordinary reverser ranking:
  1. EXISTING-DRUG filter — keep only launched / late-clinical compounds
     (a repurposing candidate must already be a real drug).
  2. INDICATION NOVELTY — compare the drug's CURRENT indication / disease area
     to the target disease. A launched drug that reverses the disease but is
     approved for something ELSE is an indication-shift opportunity; one already
     used for this disease is not "repurposing".

Data source: the **Broad Drug Repurposing Hub** annotation file
(``repurposing_drugs_YYYYMMDD.txt``), columns
``pert_iname, clinical_phase, moa, target, disease_area, indication`` — so the
approval status and current indication are REAL, not guessed.

    reversers (perturbomics)  x  Repurposing Hub annotations  ->  ranked
    launched drugs that reverse the disease from a DIFFERENT indication.

RIGOR
  ✅ rigorous : the |WTCS| reversal, the launched/phase status and current
                indication (read from the Hub), and the score arithmetic.
  ⚠️ hypothesis: "reverses the signature => will treat the new indication".
                 A repurposing HYPOTHESIS to validate — never a clinical claim.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

# clinical_phase (Hub) -> a de-risking weight; launched drugs are the prize.
_PHASE_WEIGHT = {
    "launched": 1.00, "phase 3": 0.85, "phase 2/phase 3": 0.80,
    "phase 2": 0.65, "phase 1/phase 2": 0.55, "phase 1": 0.40,
    "preclinical": 0.15, "withdrawn": 0.30,  # withdrawn: known safety, still informative
}

_SALT_WORDS = ("hydrochloride", "dihydrochloride", "sulfate", "besylate",
               "mesylate", "maleate", "sodium", "hcl", "citrate", "tartrate",
               "phosphate", "acetate", "fumarate", "succinate")


def _norm(name: str) -> str:
    n = str(name).strip().lower()
    for s in _SALT_WORDS:
        n = n.replace(" " + s, "")
    return n.strip()


def load_repurposing_hub(path: str) -> dict[str, dict]:
    """Parse a Drug Repurposing Hub ``repurposing_drugs_*.txt`` file.

    Returns ``{normalized_pert_iname: {clinical_phase, moa, target, disease_area,
    indication}}``. Leading ``!`` comment lines and the header are skipped.
    """
    out: dict[str, dict] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!") or line.startswith("pert_iname\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[0].strip():
                continue
            parts += [""] * (6 - len(parts))
            name, phase, moa, target, area, indication = parts[:6]
            out[_norm(name)] = {
                "pert_iname": name, "clinical_phase": phase.strip(),
                "moa": moa.strip(), "target": target.strip(),
                "disease_area": area.strip(), "indication": indication.strip(),
            }
    return out


def _novelty(ann: dict, disease_terms: list[str]) -> float:
    """Indication novelty in [0,1] for a launched/clinical drug.

    1.0  = different disease area entirely (true indication shift)
    0.5  = same disease AREA but a different specific indication (adjacent)
    0.0  = already indicated for the target disease (NOT a repurposing gap)
    """
    text_ind = ann.get("indication", "").lower()
    text_area = ann.get("disease_area", "").lower()
    terms = [t.lower() for t in disease_terms]
    if any(t in text_ind for t in terms):
        return 0.0                                   # already used for it
    if any(t in text_area for t in terms):
        return 0.5                                   # same area, new indication
    return 1.0                                       # different area


def screen_repurposing(
    reversers: pd.DataFrame,
    hub: dict[str, dict],
    disease_terms: list[str],
    min_weight: float = 0.4,
    require_launched: bool = False,
    top_n: int = 25,
) -> pd.DataFrame:
    """Rank existing drugs that reverse the disease from a new indication.

    Parameters
    ----------
    reversers : output of ``rank_reversers`` (needs ``name``, ``connectivity``;
        drug ``name`` should match a Hub ``pert_iname``).
    hub : output of ``load_repurposing_hub``.
    disease_terms : words identifying the TARGET disease / its area, e.g.
        ``["pulmonary", "fibrosis", "respiratory"]`` — used for the novelty check.
    min_weight : drop drugs below this phase weight (0.4 ≈ Phase 1+; use higher
        to demand later-stage). ``require_launched`` forces launched-only.
    top_n : rows to return.

    Score = |WTCS| (reversal) × phase_weight (de-risking) × indication_novelty.
    Only drugs found in the Hub AND reversing (connectivity<0) are considered.
    Returns a ranked DataFrame with the drug's real current indication + moa.
    ✅ arithmetic / annotations; ⚠️ the therapeutic hypothesis.
    """
    rows = []
    for _, r in reversers.iterrows():
        if r["connectivity"] >= 0:                   # must reverse the disease
            continue
        ann = hub.get(_norm(r["name"]))
        if ann is None:
            continue
        phase = ann["clinical_phase"].lower()
        w = _PHASE_WEIGHT.get(phase, 0.0)
        if require_launched and phase != "launched":
            continue
        if w < min_weight:
            continue
        nov = _novelty(ann, disease_terms)
        if nov == 0.0:                               # already indicated -> not repurposing
            continue
        reversal = min(abs(r["connectivity"]), 1.0)
        score = reversal * w * nov
        rows.append({
            "drug": ann["pert_iname"],
            "clinical_phase": ann["clinical_phase"],
            "moa": ann["moa"],
            "target": ann["target"].split("|")[0] if ann["target"] else "",
            "current_indication": ann["indication"] or ann["disease_area"],
            "reversal": round(reversal, 3),
            "novelty": nov,
            "repurpose_score": round(score, 4),
        })
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows).drop_duplicates("drug")
    return df.sort_values("repurpose_score", ascending=False).head(top_n).reset_index(drop=True)
