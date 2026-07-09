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

## Important: there is no DTU REST API

`services.healthtech.dtu.dk` (former CBS) does **not** expose a public REST/JSON
API, and the old SOAP/WSDL endpoints were retired. The officially supported way
to automate these tools is to **download the free-for-academia stand-alone
command-line packages** and run them locally. This project is the orchestration
layer around those binaries:

| Route | What it is | Verdict |
|-------|-----------|---------|
| **Stand-alone CLI packages** (what this repo uses) | Per-tool tarballs, run locally, batch-friendly | ✅ Supported, robust, this is the right way |
| Scripting the web form (POST + scrape) | Emulate the website | ⚠️ Against the terms of use, fragile, rate-limited — avoid |

Each tool's download tab is linked from its service page (e.g.
<https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/>). Academic users get
a free licence; commercial users email `health-software@dtu.dk`.

---

## Install the underlying tools

The pipeline runs whatever is on your `PATH`. Install only the models you need;
missing ones are **skipped, not fatal**, so you can dry-run first.

```
python -m pipeline.cli --list-models      # show the registry
```

| key | tool | binary expected | category |
|-----|------|-----------------|----------|
| `netmhcpan`   | NetMHCpan-4.1   | `netMHCpan`        | MHC-I (CD8) |
| `netctlpan`   | NetCTLpan-1.1   | `netCTLpan`        | MHC-I / CTL |
| `netmhciipan` | NetMHCIIpan-4.3 | `netMHCIIpan`      | MHC-II (CD4) |
| `bepipred`    | BepiPred-3.0    | `bepipred3_CLI.py` | B-cell linear |
| `signalp`     | SignalP-6.0     | `signalp6`         | aux |

After unpacking a tarball, either add its directory to `PATH` or point the
pipeline straight at the executable:

```
--binary netMHCpan=/opt/netMHCpan-4.1/netMHCpan
```

> DiscoTope-3.0 (conformational B-cell epitopes) needs a **3D structure (PDB)**,
> not a bare sequence — it is registered but marked `needs_structure` until a
> `--structure` input is wired in.

---

## Usage

```bash
# Dry run — proves the orchestration; skips tools you haven't installed yet
python -m pipeline.cli --fasta examples/example.fasta --out results/

# Real run once tools are on PATH
python -m pipeline.cli \
    --fasta my_protein.fasta \
    --models netmhcpan,netctlpan,netmhciipan,bepipred \
    --alleles-i  HLA-A02:01,HLA-A01:01,HLA-B07:02 \
    --alleles-ii DRB1_0101,DRB1_0401,DRB1_1501 \
    --lengths-i  8,9,10,11 \
    --lengths-ii 15 \
    --workers 8 \
    --out results/
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
  command builder + an output parser. **Add a new DTU tool by adding one entry**
  (a `_cmd_*` builder and a `_parse_*` function) — no other file changes.
- `pipeline/runner.py` — parallel executor (`ThreadPoolExecutor`); resolves
  binaries, runs each model in its own workdir, captures stdout, never lets one
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
