# Data access — fetching exact literature models & real data

The difference between "reproduced the topology" and "reproduced the paper to
the decimal" is real curated data. Here is what works behind the agent proxy
and what is blocked.

## Exact curated models: BioModels GitHub mirror + libRoadRunner

**This is the working recipe.** BioModels is mirrored one-repo-per-model on the
`biomodels` GitHub org, and `raw.githubusercontent.com` is reachable through the
proxy. Combined with libRoadRunner (an SBML simulator) this loads and integrates
any curated model with **zero transcription error**.

```python
MIRROR = "https://raw.githubusercontent.com/biomodels/{id}/master/{id}/{id}.xml"
# note the NESTED path: {id}/{id}.xml
# fetch with curl (see m20b_biomodels_exact.fetch_sbml), cache under figures/_sbml_cache/
import roadrunner
rr = roadrunner.RoadRunner(path)
rr["MAPKK"] = 50            # set a parameter/species
rr.reset(); rr.simulate(0, 200000, 400)
mpp = rr["Mpp"]            # read a species
lam = max(__import__("numpy").linalg.eigvals(rr.getFullJacobian()).real)  # λ_max
```

Worked examples in `m20b_biomodels_exact.py`:
- **Markevich2004 MAPK** = `BIOMD0000000027`. Confirms the hand-coded M15:
  official `Km5 = 78.0` (the value M15 reverse-engineered), Mpp OFF=49.42 /
  ON=481.95 at MAPKK=50 — matches to the decimal. This is how you *independently
  validate* a hand-coded model.
- **Legewie2006 apoptosis** = `BIOMD0000000102`. Bistable life/death caspase
  switch; the control knob is XIAP synthesis `k18prod` (window ≈ [0.08, 0.14]).

To find a model's BIOMD id: search PubMed / the BioModels curated list for the
paper, or look for the `biomodels/BIOMD…` repo. Set the death/ON basin
explicitly (initial conditions) to map both branches of a bistable model.

## What is BLOCKED (don't waste time)

- `biomodels.org` and the EBI download endpoints → 403 through the proxy.
- `api.github.com` → scoped/limited; use `raw.githubusercontent.com` (and
  jsDelivr) for raw file content instead.
- The session's GitHub MCP is scoped to specific repos only.
- PubMed open-access full text exists only for models with a PMCID — some
  classic papers (e.g. Yao 2008, Eissing 2004) have **no** PMCID, so their exact
  rate constants are not fetchable; reproduce topology + documented behavior and
  label it ⚠️ (or find the model on the biomodels mirror instead, as with M20b).

## Real single-cell / time-series data (for validation)

- Pertz-lab single-cell EKAR traces ship with `dmattek/shiny-timecourse-inspector`
  and are fetchable from `raw.githubusercontent.com` (see `m17_realdata`).
- **Honest lesson from M17:** to decisively confirm a critical-slowing biomarker
  you need a dataset where a **dose gradient crosses the switch threshold** (some
  cells sitting near the bifurcation). Public ERK datasets checked
  (Pertz EGF/FGF = all supra-threshold; Goglia `idr0064` = single 2.5 µM dose,
  single-cell trajectories not deposited; Albeck-Lab = analysis code only) did
  **not** provide this. The decisive experiment is a MEKi titration gradient with
  single-cell ERK-KTR/EKAR imaging crossing the OFF threshold. The M16→M18
  pipeline is ready to judge such data on arrival; the prediction is a ~5×
  variance/AR1 peak at near-threshold doses (M18 positive control).

## Drug / target / trial data (MCP tools, when connected)

ChEMBL (targets, IC50, mechanism), Boltz-2.1 (structure + binding — remember to
disable the reactive-group SMARTS filter for covalent warheads), Inductive Bio
(physchem), ClinicalTrials (endpoints), PubMed / bioRxiv (methods). All were
used live in M5–M10. Validate computational rankings against measured values
(M10): use Boltz `optimization_score` (tracks potency, ρ≈+0.6), **not**
`binding_confidence` (pose confidence, ρ≈−0.2).
