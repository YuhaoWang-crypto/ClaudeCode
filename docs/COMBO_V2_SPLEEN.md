# Expanded library (v2) — spleen deep-screen

Enumerate: `analysis/enumerate_library_v2.py` · Ranked: `data/libraries/combo_v2/combo_v2_spleen_ranked.csv`
· Figure: `results/figures/combo_v2_spleen.png`

## What was added
On top of v1: **4 degradable/reduction-cleavable tails** (2 disulfides `T17/T18`,
2 internal-ester tails `T19/T20`) + the pre-existing thioether/ether hetero tails,
and **4 new heads** (polar bis-hydroxyethyl `H15`, longer-spacer dimethylaminobutyl
`H16`, rigid cyclohexane-diamine `H17`, large polyamine `H18`).

18 heads × {ester, amide} × 20 tails → **538 valid lipids** (170 degradable-tail,
50 disulfide). Deep-screened for **spleen** (spleen/IV/mouse context) with the
5-fold ensemble.

Spleen is the chosen target because it had the most distinct SAR in the organ
comparison (weakest transfer from in-vitro HeLa, ρ=0.60) and is Su 2026's headline
application (spleen-tropic mRNA vaccines).

## Do the degradable tails help spleen? — yes, modestly

Mean predicted spleen delivery by tail chemistry:

| tail chemistry | mean pred (z) |
|---|---|
| branched | **0.173** |
| thioether (–S–) | 0.083 |
| disulfide (–S–S–) | 0.032 |
| internal-ester | −0.003 |
| plain alkyl | −0.029 |
| ether (–O–) | −0.032 |
| unsaturated | **−0.165** |

Sulfur-containing degradable tails (thioether, disulfide) **beat plain alkyl** for
spleen, and branched is best; unsaturated is worst (opposite of the in-vitro HeLa
result). So adding reduction-cleavable tails is worthwhile here — you gain
biodegradability/tolerability *and* a small predicted-potency lift.

SAR (panel C): **amide ≫ ester**, and **2-tail ≈ 3-tail ≫ 4-tail** for spleen.

## Spleen shortlist (confidence-aware: mean ± cross-fold std)

| rank | head | linker | tail | n_tails | MW | pred | note |
|---|---|---|---|---|---|---|---|
| 1 | H7 (N-Me-DAP) | amide | C9 | 3 | 680 | **0.57 ± 0.21** | top, high conf |
| — | H7 | amide | **thioether-C8** | 3 | 818 | **0.49 ± 0.22** | best degradable |
| — | H2 (DEAE-amine) | amide | thioether-C8 | 2 | 603 | 0.44 ± 0.30 | small + degradable |
| — | H7 | amide | **disulfide-C6** | 3 | 999 | 0.32 ± 0.21 | best disulfide, high conf |
| — | H2 | amide | disulfide-C6 | 2 | 723 | 0.32 ± 0.31 | small + disulfide |

New head **H16 (dimethylaminobutyl)** performs well for spleen (mean 0.22, in the
top-15); the polar `H15` and rigid `H17` heads do not.

## Recommended spleen leads
- **Potency-first**: `H7 + amide + C9` (3-tail) — highest, tight ensemble.
- **Degradable/tolerability-first**: `H7 + amide + thioether-C8` (near-top,
  high conf) or `H7 + amide + disulfide-C6` (reduction-triggered release, high conf).
- **Small/simple**: `H2 + amide + C9 / thioether-C8` (2-tail, MW ~510–600).

## Caveat
Same as ORGAN_SCREEN.md: model conditioning on a context with sparser in-vivo
training data — treat as ranked hypotheses, favour the low-std picks. Validate the
disulfide leads' redox behaviour experimentally.
