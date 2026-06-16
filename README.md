# Allosteric distal-pocket → His-pair/His-triplet metal-switch pipeline

Workflow for the project brief: take an enzyme whose **active site changes state**
(open↔closed, apo↔ligand), align the two conformers, find the **5–20 distal
regions that respond**, then design a small **His-pair / His-triplet** library at
the top hinge/interface loops to screen **Zn²⁺ / Ni²⁺ / Co²⁺** effects on
`kcat`, `Km`, `kcat/Km`.

```
align conformer pair  →  per-residue Cα displacement (rigid-core fit)
                      →  rank distal candidate regions + hinge residues
                      →  His-pair/His-triplet library (metal-geometry scored)
                      →  results/<TARGET>.md report
```

## Quick start

```bash
pip install -r requirements.txt
python src/run_pipeline.py            # all targets in targets.yaml
python src/run_pipeline.py AdK        # one target
```

Reports land in `results/<TARGET>.md`; `results/SUMMARY.md` is the overview.

## Targets (`targets.yaml`)

| target | enzyme | pair (ref ⇄ alt) | role |
|--------|--------|------------------|------|
| GCK    | human glucokinase | 1V4T ⇄ 1V4S | monomer positive control, big domain motion + activator pocket |
| PTP1B  | protein tyrosine phosphatase 1B | 2HNP ⇄ 1T49 | known distal allosteric inhibitor site (~20 Å) |
| AdK    | *E. coli* adenylate kinase | 4AKE ⇄ 1AKE | **de-novo Zn-switch object** (large open/closed, small) |
| PFK    | *E. coli* phosphofructokinase-1 | 2PFK ⇄ 1PFK | oligomeric allostery benchmark |
| ATCase | *E. coli* aspartate transcarbamoylase | 6AT1 ⇄ 1D09 | strong allostery, very large assembly |
| GP     | human liver glycogen phosphorylase | 1FA9 ⇄ 3CEH | drug-type AMP-site benchmark |
| cpTEM1 | circularly-permuted TEM-1 β-lactamase | 1BTL (seq-only) | source-paper reproduction |

Edit `active_site` / `known_allosteric` / cutoffs in `targets.yaml` freely.

## Method notes

- **Rigid-core superposition** (Kabsch + iterative outlier trimming): superposing
  on the whole protein smears a domain motion across every residue, so we
  auto-detect the rigid core and measure displacement against *that* frame.
- **Distal candidate region** = contiguous residues that are far from the
  active-center centroid (`distal_cutoff_A`, default 15 Å) **and** move more than
  the `displacement_pctl` (default P75) threshold.
- **Hinge residues** = distal residues with the steepest local displacement
  gradient — where a moving domain meets the rigid core; the best single-metal
  clamp points.
- **His library** = backbone-geometry heuristics (Cα–Cα window + converging side
  chains + solvent exposure). **Coarse pre-filter only** — validate every site
  with rotamer-level metal modelling (Rosetta / MIB / HADDOCK-metal) before
  synthesis. Regulation direction (activate vs inhibit) is geometry-dependent and
  must be measured per construct, not assumed.

## ⚠️ Network / data access in this environment

The structural-data hosts are **not in the network egress allowlist**, so direct
RCSB download returns `HTTP 403 "Host not in allowlist"`:

```
files.rcsb.org   www.ebi.ac.uk   rest.uniprot.org   files.pdbj.org   alphafold.ebi.ac.uk   → 403
pypi.org   github.com   raw.githubusercontent.com                                          → 200
```

Two ways to get the other six targets running:

1. **Allowlist a host** — add `files.rcsb.org` (and optionally `www.rcsb.org`,
   `data.rcsb.org`) to the environment's network egress settings, then re-run.
   See https://code.claude.com/docs/en/claude-code-on-the-web .
2. **Drop files in manually** — put `1V4T.pdb`, `1V4S.pdb`, … into `data/raw/`.

**AdK runs offline today**: its conformer pair (`4AKE` open, `1AKE` closed) is
sourced from the ProDy-bundled test data, so the demo works with no network and
proves the pipeline end-to-end.

## Layout

```
src/fetch_structures.py   resolve a PDB id → local file (cache → bundled → RCSB/PDBe/PDBj mirrors)
src/analyze_conformers.py alignment, rigid-core fit, displacement, distal + hinge ranking
src/design_his_sites.py   His-pair / His-triplet library with metal-geometry scoring
src/run_pipeline.py       driver → results/*.md
targets.yaml              curated target table (PDB ids, chains, active site, known allosteric)
```
