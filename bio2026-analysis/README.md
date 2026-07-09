# BIO2026 Next-Wave Signal Analysis

Seller's-market bias-corrected read of the BIO2026 asset list (11,976 assets)
to infer the next licensing/BD wave, weighting toward early-stage and academic
(university Phase I/II) signals.

## Method

The BIO asset list over-represents whatever is currently fundable/licensable, so
raw volume is a *lagging/coincident* indicator. Each theme's Next-Wave score is a
weighted sum of standardized terms:

- **+ immaturity gap** (early − late phase share) — rising vs peaking shape
- **+ early-phase enrichment** — science leads commerce ~2–4 years
- **+ academic-origin share** — least contaminated by licensing fashion
- **− crowding penalty** log(N) — a full shelf = today's crowded trade
- **− partnered rate** — already buyer-validated = current wave, not next

## Scripts

Run in order (each writes `df.pkl` / `scored.pkl` intermediates, git-ignored):

1. `analyze.py` — load Excel, classify academic origin, mine modality keywords → `df.pkl`
2. `score.py` — compute the bias-corrected Next-Wave score per theme → `scored.pkl`
3. `deep.py` — convergence intersections + early-academic target mining

```bash
pip install openpyxl pandas numpy
python3 analyze.py && python3 score.py && python3 deep.py
```

Input Excel (`AssetsExport.xlsx`) is not committed — point the `F=`/`read_excel`
path in the scripts at your local copy.

## Report

`bio2026-nextwave.html` — self-contained interactive report (leaderboard,
convergence plays, early-academic target tags, what's peaking, caveats).

## Caveats

Keyword-regex aggregation is noisy; seller-self-reported phase/indication;
"next wave" here is for BD/licensing buyers, not the public market; single
snapshot with no year-over-year growth (a BIO2025 diff would be stronger).
