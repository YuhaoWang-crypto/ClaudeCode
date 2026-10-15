# Conformational dynamics in pulsed electron transfer — three-route comparison (MB–DNA, SWV)

**Question.** With the same data and the same electron-transfer assumption, does keeping conformational
exchange (R1, R2) predict square-wave current better than a static / fast-equilibrium average (R0)?

This round implements a real dynamic current backend: the joint conformation × redox master equation is
propagated through the exact SWV waveform, and current comes only from redox jumps. Contact probability
is never used as a current. Pre-registered protocol: [`PROTOCOL.md`](PROTOCOL.md), frozen in commit
`c05c2a4` before any fitting. Data provenance and digitization: [`data/dauphin2017/MANIFEST.md`](data/dauphin2017/MANIFEST.md).

<!-- RESULTS -->

## Pipeline

| step | file | what it does |
|---|---|---|
| data | `digitize_2017.py` | reads the experimental lines (not the fit circles) of SI Fig. S-7, per-panel axis calibration, reading errors, overlay audit |
| backend | `swv_backend.py` | exact half-pulse propagation of the deviation from Nernst equilibrium; spectrum from an SVD of the generator's square-root factor (stable when exchange is 10⁸–10¹³× faster than ET); R0, R1, generic reversible generators, two-population control |
| tests | `tests/test_backend.py` | vs full-system `expm`, vs stiff Radau integration over a whole scan, fast/slow exchange limits, generator properties, Nernst hold, CTMC generator recovery |
| molecule | `cg_md.py`, `md_campaign.py` | OpenMM CPU Brownian dynamics, 1 bead/nt, 64 independent replicas per context, 4 seeds, mixed extended/collapsed starts |
| time step | `md_analysis.py` | selection rule: fastest dt whose halving changes ⟨z⟩, ⟨R_g⟩, P(contact), ⟨e^{−β(z−z_min)}⟩, contacts and τ_z by ≤5 % |
| R2 | `msm.py`, `r2_build.py` | k-means states on (z, R_g, r_end, contacts), reversible generator by lag-τ likelihood, timescale / CK / sampling / block-PMF gates |
| fits | `fit_cv.py`, `run_analysis.py` | fit A (shared θ) and fit B (per-route + global s), frequency and amplitude folds, β and grid sensitivity, profiles |
| report | `summarize.py` | tables (`results/tables_dt*.md`) and figures (`figures/`) |

Reproduce (about 2–3 h on 4 cores):
```
python mb_dynamics/digitize_2017.py
python -m pytest mb_dynamics/tests -q
cd mb_dynamics && python md_campaign.py dt --dts 1.0,0.5,0.25 && python md_analysis.py
python r2_build.py <dt> && python run_analysis.py --dt <dt> && python summarize.py <dt>
```
