# Multi-model epitope / immunogenicity pipeline (DTU Health Tech)

Give it **one protein FASTA**, it pushes the sequence through **many epitope /
immunogenicity predictors in parallel** (MHC-I binding, MHC-II binding, CTL
epitopes, B-cell epitopes, …) and merges every model's output into **one
normalized table plus a consensus ranking** of the most immunogenic peptides.

```
protein.fasta ──┬─► NetMHCpan   (MHC-I binding, CD8)      ─┐
                ├─► NetCTLpan   (CTL: MHC+TAP+cleavage)    ─┤
                ├─► NetMHCIIpan (MHC-II binding, CD4)      ─┼─► normalized CSV
                ├─► BepiPred    (linear B-cell epitopes)   ─┤   + consensus
                └─► SignalP …   (aux filtering)            ─┘   epitope ranking
        (all models run concurrently; failures/uninstalled tools are skipped)
```

---

## How to actually run these tools programmatically

`services.healthtech.dtu.dk` (former CBS) exposes **no public REST/JSON API**
and the old SOAP/WSDL endpoints were retired. There are two real automation
routes, and this pipeline supports **both** via a per-model `backend`:

| Route (`backend`) | What it is | Needs | Verdict |
|---|---|---|---|
| **BioLib cloud** (`biolib`) | Some DTU tools are published on [biolib.com](https://biolib.com/DTU/) and run in the cloud via the `pybiolib` Python API — no download, no local license | `pip install pybiolib` + a free BioLib account | ✅ Easiest; **a real API** — use it wherever the tool exists |
| **Stand-alone CLI** (`cli`) | Per-tool tarball from the DTU download page, run locally | Academic licence download; binary on `PATH` | ✅ Required for licensed tools not on BioLib |
| ~~Scrape the web form~~ | Emulate the website | — | ⚠️ Against the terms of use, fragile — avoid |

**The catch: not every tool is on BioLib.** Verified against the live BioLib API:

| Tool | On BioLib? | app id | Backend to use |
|------|-----------|--------|----------------|
| BepiPred-3.0 (B-cell linear) | ✅ | `DTU/BepiPred-3` | `bepipred-cloud` |
| DiscoTope-3.0 (B-cell conformational) | ✅ (needs PDB) | `DTU/DiscoTope-3` | `discotope-cloud` |
| NetSurfP-3.0 (surface/2° structure) | ✅ | `DTU/NetSurfP-3` | `netsurfp-cloud` |
| DeepTMHMM (transmembrane) | ✅ | `DTU/DeepTMHMM` | `deeptmhmm-cloud` |
| **NetMHCpan / NetMHCIIpan / NetCTLpan** (MHC binding) | ❌ **not published** (separately licensed) | — | `cli` (local download) |

So the practical setup is **hybrid**: B-cell / structural predictors run on
BioLib cloud with zero install, while the licensed MHC binders run from the
local CLI packages. One orchestrator, one merged report.

---

## Setup

```
pip install pybiolib          # for the cloud backend
biolib login                  # or export BIOLIB_TOKEN=...  (free account)
python -m pipeline.cli --list-models      # show the registry + each backend
```

Install only what you need. Missing tools are **skipped, not fatal**, so you can
dry-run first. Registry:

| key | tool | backend | how to enable | category |
|-----|------|---------|---------------|----------|
| `bepipred-cloud`  | BepiPred-3.0    | cloud | `pip install pybiolib` | B-cell linear |
| `discotope-cloud` | DiscoTope-3.0   | cloud | pybiolib + `--structure` (PDB) | B-cell conformational |
| `netsurfp-cloud`  | NetSurfP-3.0    | cloud | `pip install pybiolib` | aux (surface/2°) |
| `deeptmhmm-cloud` | DeepTMHMM       | cloud | `pip install pybiolib` | aux (transmembrane) |
| `netmhcpan`   | NetMHCpan-4.1   | cli   | `netMHCpan` on PATH        | MHC-I (CD8) |
| `netctlpan`   | NetCTLpan-1.1   | cli   | `netCTLpan` on PATH        | MHC-I / CTL |
| `netmhciipan` | NetMHCIIpan-4.3 | cli   | `netMHCIIpan` on PATH      | MHC-II (CD4) |
| `bepipred`    | BepiPred-3.0    | cli   | `bepipred3_CLI.py` on PATH | B-cell linear (local) |
| `signalp`     | SignalP-6.0     | cli   | `signalp6` on PATH         | aux |

For a **licensed local tool**, unpack its tarball and either add its dir to
`PATH` or point the pipeline straight at the executable:

```
--binary netMHCpan=/opt/netMHCpan-4.1/netMHCpan
```

---

## Usage

```bash
# Cloud-only run — zero local install, just pybiolib + a BioLib account.
# Runs B-cell + structural predictors in the BioLib cloud, in parallel.
python -m pipeline.cli \
    --fasta my_protein.fasta \
    --models bepipred-cloud,netsurfp-cloud,deeptmhmm-cloud \
    --out results/

# Hybrid run — cloud B-cell/structure + local licensed MHC binders together.
python -m pipeline.cli \
    --fasta my_protein.fasta \
    --models netmhcpan,netctlpan,netmhciipan,bepipred-cloud \
    --alleles-i  HLA-A02:01,HLA-A01:01,HLA-B07:02 \
    --alleles-ii DRB1_0101,DRB1_0401,DRB1_1501 \
    --lengths-i  8,9,10,11 --lengths-ii 15 \
    --workers 8 --out results/

# Conformational B-cell epitopes need a structure:
python -m pipeline.cli --fasta my_protein.fasta --structure my_protein.pdb \
    --models discotope-cloud --out results/
```

### Outputs (`results/`)

| file | contents |
|------|----------|
| `all_predictions.csv` | every peptide × model × allele, normalized to one schema |
| `consensus_epitopes.csv` | peptides ranked by how many models/categories flag them |
| `run_summary.json` | per-model status, record counts, timings, exact command run |

The **consensus** step is the payoff: a peptide called a strong binder by
NetMHCpan *and* NetCTLpan, or one sitting in both an MHC-I hotspot and a
BepiPred B-cell region, ranks above a peptide only one model liked.

---

## How it's built (extending it)

- `pipeline/models.py` — the model **registry**. Each `ModelSpec` bundles a
  backend (`cli` or `biolib`), a command/arg builder, and an output parser.
  **Add a tool by adding one entry** — a `cli` model needs a `_cmd_*` builder +
  `_parse_*` (stdout); a `biolib` model needs `app_uri` + `_args_*` +
  `_pfiles_*` (output files). No other file changes.
- `pipeline/biolib_backend.py` — thin `pybiolib` wrapper: `biolib.load(app_uri)`
  → `app.cli(args=...)` → `job.save_files()`. Imported lazily so the CLI-only
  path never requires pybiolib.
- `pipeline/runner.py` — parallel executor (`ThreadPoolExecutor`); dispatches
  each model to its backend, runs it in its own workdir, and never lets one
  tool's failure sink the run.
- `pipeline/aggregate.py` — normalization, %Rank binder-calling (DTU thresholds:
  MHC-I SB ≤0.5 / WB ≤2.0; MHC-II SB ≤1.0 / WB ≤5.0), consensus ranking, writers.
- `pipeline/cli.py` — argument parsing / entry point.

## Test

```bash
python tests/test_pipeline.py        # or: python -m pytest tests/
```

Tests use `tests/mock_netMHCpan`, a stand-in that emits real NetMHCpan-4.1
column layout, so the full path (command build → parse → consensus) is verified
without any licensed binary.

## Optional: cross-check with hosted models

This environment also exposes hosted immunogenicity/binding models via MCP
(EDEN `predict_immunogenicity`, Boltz structure+binding). They can be run as an
independent cross-check against the DTU consensus — see `INTEGRATIONS.md`.
