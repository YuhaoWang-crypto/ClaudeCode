"""
calibration.py -- turn a raw Boltz interface score into a CALIBRATED triage signal,
grounded on real two-CRO wet-lab labels from the Anthropic claude-protein-binder-design
release (see reference/boltz-score-calibration.md; reproduced in binder_audit/run_audit.py).

Key facts this encodes (all measured, not assumed):
  * Boltz-2 ipSAE_min triages binders WITHIN a target at ~1.5x enrichment (mAP 0.462,
    beats chance on 13/13 targets). Use WITHIN-PANEL PERCENTILES, not absolute cutoffs.
  * The score predicts WHETHER a design binds, NOT how tightly (score vs logKD rho=-0.05).
    So it is a feasibility / ordering signal, never an affinity or dynamic-range proxy.
"""
from __future__ import annotations

# Reproduced Boltz-2 within-target enrichment (binder_audit/audit_boltz2_calibration.json):
# keep top X% of a target's designs by score -> hit rate. Baseline (no filter) = 26.8%.
BOLTZ2_ENRICHMENT = [           # (keep_top_fraction, hit_rate_pct)
    (1.00, 26.8), (0.50, 35.3), (0.30, 38.9), (0.10, 40.2), (0.05, 42.6),
]
BASELINE_HIT_RATE = 26.8
FEASIBILITY_NOTE = ("feasibility/triage only — Boltz score predicts WHETHER it binds, "
                    "NOT affinity or dynamic range (wet-lab titration is ground truth)")


def expected_hit_rate(top_fraction: float) -> float:
    """Interpolate expected binder hit-rate (%) for keeping the top `top_fraction` of a
    panel by Boltz score, from the reproduced enrichment curve. top_fraction in (0,1]."""
    if not (0 < top_fraction <= 1):
        raise ValueError("top_fraction must be in (0, 1]")
    pts = sorted(BOLTZ2_ENRICHMENT)                      # ascending fraction
    for (f0, h0), (f1, h1) in zip(pts, pts[1:]):
        if f0 <= top_fraction <= f1:
            w = (top_fraction - f0) / (f1 - f0) if f1 > f0 else 0.0
            return round(h0 + w * (h1 - h0), 1)
    return round(pts[-1][1] if top_fraction > pts[-1][0] else pts[0][1], 1)


def calibrated_rank(candidates: list, score_key: str = "ligand_iptm") -> list:
    """Rank candidates within a panel and annotate each with a calibrated triage read.

    Returns a new list (highest score first). Each item gains:
      rank, within_panel_percentile (1.0=best), keep_top_fraction (= 1 - pct + 1/n, i.e.
      the smallest top-fraction that still includes this candidate),
      expected_hit_rate_pct (from the reproduced curve), triage_note.
    Scores that are None sort last.
    """
    scored = [c for c in candidates if c.get(score_key) is not None]
    unscored = [c for c in candidates if c.get(score_key) is None]
    scored.sort(key=lambda c: c[score_key], reverse=True)
    n = len(scored)
    out = []
    for i, c in enumerate(scored):
        keep_frac = (i + 1) / n                          # top-fraction that includes rank i+1
        item = dict(c)
        item["rank"] = i + 1
        item["within_panel_percentile"] = round(1 - i / n, 3) if n > 1 else 1.0
        item["keep_top_fraction"] = round(keep_frac, 3)
        item["expected_hit_rate_pct"] = expected_hit_rate(keep_frac)
        item["enrichment_vs_baseline"] = round(item["expected_hit_rate_pct"] / BASELINE_HIT_RATE, 2)
        item["triage_note"] = FEASIBILITY_NOTE
        out.append(item)
    for c in unscored:
        item = dict(c); item["rank"] = None; item["triage_note"] = "no score"
        out.append(item)
    return out
