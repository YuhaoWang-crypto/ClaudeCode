---
name: pathway-biomarker-triangulation
description: Discover and CROSS-VALIDATE the core biomarkers of a pathway / target set with three orthogonal methods, then map them to clinical use (enrollment / progression-PD / efficacy). Use when asked to mine or screen biomarkers for a gene set or pathway, build a companion-diagnostic (CDx) enrichment panel, or sanity-check a biomarker list. Runs Method 1 = evidence×measurability attribution + greedy minimal-essential-unit cover (Open Targets + STRING + measurability priors); Method 2 = random-walk-with-restart network diffusion (pure STRING topology); and hands off to Method 3 = the user's `network-biomarker` skill (grn_pipeline: irreducibility/CRNT/DNB critical-slowing early-warning). Convergence across methods = the defensible core. Triggers - "find/screen/mine biomarkers", "pathway biomarker", "CDx panel", "minimal core biomarkers", "enrollment/PD/efficacy biomarkers", "cross-validate biomarkers".
---

# Pathway biomarker triangulation

Three methods answer three DIFFERENT questions about a pathway's biomarkers. Where they
converge is the most defensible core; where they diverge is informative (method signature).

| Method | Question | Nature | Adds |
|---|---|---|---|
| **1 Attribution** | which *measurable* node has *disease evidence*? | static · evidence | the clinical CDx panel |
| **2 RWR diffusion** | which node is *network-central* to the seeds? | static · topology | topology sanity-check |
| **3 network-biomarker skill** | which module is *dynamically irreducible* + which observable *warns of a tipping point*? | dynamical · systems | switch capacity (CRNT δ) + DNB early-warning |

## Run Methods 1 + 2 (this skill)
```bash
pip install networkx matplotlib
python scripts/triangulate.py \
  --seeds GPR75,INHBE,ACVR1C,ACVR2B,MSTN,GDF15,GFRAL \
  --disease MONDO_0011122 --out out/          # obesity; use any EFO/MONDO id
# outputs: out/triangulation.tsv, out/triangulation.png, out/triangulation.json
```
- **Method 1 (attribution)**: `attribution = measurability × (0.65·disease_assoc + 0.35·centrality)`, then a **greedy minimal cover** takes the smallest node subset reaching ≥80% of total attribution = the *minimal essential biomarker unit*. Disease association = Open Targets Platform GraphQL; centrality = STRING-subnetwork betweenness+degree; measurability = analyte-class prior (secreted/genotype 1.0 · phospho-node 0.8 · receptor/occupancy 0.7 · intracellular 0.3) applied as a **multiplicative gate** ("if you can't measure it, it isn't a biomarker").
- **Method 2 (RWR)**: personalized PageRank (restart r=0.3) from the seeds on the STRING network — pure topology, no disease/measurability. Take top-k (= minimal-set size). Report overlap + Jaccard vs Method 1.

Tune with `--assoc-w`, `--coverage`, `--rwr-r`, `--string-score`. Edit the `MEAS/SECRETED/RECEPTOR/PHOSPHO` sets in the script for your analytes.

## Method 3 hand-off (the user's `network-biomarker` skill)
For the *dynamical* layer — irreducible core (graph automorphism/fibration), switch capacity (CRNT deficiency δ), and a **DNB / critical-slowing early-warning biomarker** (variance↑ / AR1↑ / leading-eigenvalue→0 at a saddle-node) — run the separate **`network-biomarker`** skill (`grn_pipeline` package) on the same pathway: add it via `reference/adding-a-pathway.md`, reuse engines `m1_symmetry`/`m11_fibration` (core), `m2_crnt`/`m12_dualphos` (δ), `m4_dnb_lyapunov`/`m19` (DNB). It uniquely produces switch-capacity + a tipping-point observable that Methods 1–2 structurally cannot.

## Clinical mapping (of the convergent core)
- **Enrollment (入组)**: genotype/eQTL + baseline-measurable (circulating ligands, imaging).
- **Progression / PD (进程)**: central node — ligand knockdown ↓, phospho-node (pSMAD2/3), receptor occupancy; set re-dose cadence off suppression kinetics.
- **Efficacy (效果确认)**: downstream phenotype (imaging / body-composition / disease lab).

## Interpreting results (worked obesity example — see references/methodology.md)
On the obesity activin/SMAD + GDF15 set: Method 1 minimal core (13/24) and Method 2 RWR top-13 shared **GDF15, ACVR2B, GFRAL, MSTN, RET, FST, ACVR1** (Jaccard 0.37); Method 3 gave irreducible core `ligand*→ACVR2B→ACVR1C→SMAD2/3→SMAD4`, CRNT δ=2 (bistable-capable), and a saddle-node **pSMAD2/3 variance/AR1/DNB** early-warning. **All three converge on ACVR2B→SMAD2/3.**

## Honest caveats
- Method 1 disease-association is an algorithmic prior — a target missing from Open Targets (e.g. INHBE's recent pLOF signal not yet ingested) scores 0 and must be **rescued by other evidence** (e.g. AlphaGenome splice prediction). This is *why* triangulation matters: no single method is sufficient. `references/example_obesity_triangulation.tsv` is the worked output.
- GPR75-type orphan nodes with **no STRING interactions are invisible to Method 2** (topology) — only Method 1's evidence layer sees them. Never rely on topology alone.
