# Data manifest — Dauphin-Ducharme et al. 2017, DNA20–MB SWV (frozen before any fitting)

**Source.** P. Dauphin-Ducharme, N. Arroyo-Currás, M. Kurnik, G. Ortega, H. Li, K. W. Plaxco,
"Simulation-Based Approach to Determining Electron Transfer Rates Using Square-Wave Voltammetry",
*Langmuir* 2017, 33, 4407–4413. doi:10.1021/acs.langmuir.7b00359 (PMC5660319, NIH author manuscript).
SI Figure S-7 (page S-8 of `NIHMS881272-supplement-supplemental_information.pdf`).

| item | value | source |
|---|---|---|
| probe | 5′-HS-(CH₂)₆-ATTATTTTTTATTTATTTTT-MB-3′ (MB on a 3′ thymine) | main text, Methods (rendered formula) |
| SAM | 6-mercapto-1-hexanol, overnight | Methods |
| electrode | 2 mm Au disc, ~0.03 cm² electroactive (average) | Methods |
| reference | Ag/AgCl; Pt counter; CHI 660D | Methods |
| scan | 0 → −0.45 V, step 1 mV (reduction direction) | Methods |
| conditions | f ∈ {20, 100, 500} Hz × E_sw ∈ {25, 50, 75} mV (9 panels) | Fig. S-7 |
| sign | cathodic (reduction) current plotted positive; forward red, backward black, net blue | Fig. S-7 |
| temperature | 298 K (stated for the numerical model) | main text |
| buffer | **not stated for the DNA measurements**; assumed aqueous phosphate-buffered saline (η = 0.89 mPa·s at 298 K). Assumption. | — |

**What these numbers are.** Figure readings from one published raster figure (1299×1086 px JPEG,
SHA-256 `577a0052…f1269e`, kept as `figS7_source.jpeg`). The main-text Fig. 5 (100 Hz, same data) is
kept as `fig5_main_source.jpg` for reference only and is **not** digitized separately.
Each panel shows one voltammogram (no replicate electrodes shown); the nine conditions are
presumably on one electrode and are **not** nine independent electrodes. Points along a curve are
strongly correlated; pixel samples are **not** experimental replicates.

**What was dropped.** The open circles (authors' COMSOL fit, k⁰ = 0.4 cm s⁻¹, α = 0.37) are removed by
morphological erosion and never used.

**Digitization** (`digitize_2017.py`):
- per-panel axis detection; linear least-squares map through the y tick marks (tick-position
  residuals 0.0004 / 0.0016 / 0.0036 µA for the 20 / 100 / 500 Hz rows — the source figure's ticks
  are not perfectly uniform) and the two x ticks (−0.2, −0.4 V);
- colour masks → 2×2 erosion → components spanning ≥15 px → per-column median;
- resampled on a 5 mV grid, −0.060 … −0.450 V (79 points per curve);
- `sigma_read_uA` = half the local line width ⊕ max tick residual (reading error only);
- `flag`: `visible` (1407 pts), `interp` (82; interior gaps, linear), `occluded_eq_net` (637; forward/backward
  hidden under the net trace at baseline, set equal to the net trace there).

**Internal consistency check** (after the frozen baseline step below): net − (forward − backward), as a
fraction of the net peak, RMS 1.7–4.4 %, worst single point 4.9–15.9 % (at occluded edges). This is the
**digitization floor**: NRMSE differences below ~3 % between routes are not interpretable.

**Frozen preprocessing (identical for every route).** For each of the 27 curves, subtract a straight
line fitted by least squares to the edge windows E > −0.12 V and E < −0.41 V. (The 100 Hz panels carry a
+0.01–0.02 µA offset in the source figure.) No smoothing, no rescaling.
