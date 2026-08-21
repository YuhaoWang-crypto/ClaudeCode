#!/usr/bin/env python3
"""Build the distributable `neoantigen-selection` skill zip.

    python package_skill.py                 # -> dist/neoantigen-selection-skill-v1.0.0.zip

The result is self-contained: the `neoantigen_pipeline` package is vendored
under `scripts/`, the reference human proteome and the TESLA mirror ship in the
cache directory the code already looks in, and a finished demo run is included
so a reader can see the outputs without spending 40 minutes of cloud
predictions first.

What is deliberately NOT shipped: the prediction cache (~46 MB and specific to
one patient and one HLA set) and `candidates.csv` (4.3 MB of intermediate rows).
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import zipfile
from datetime import datetime, timezone

ROOT = os.path.dirname(os.path.abspath(__file__))
SKILL_SRC = os.path.join(ROOT, ".claude", "skills", "neoantigen-selection")
PKG_SRC = os.path.join(ROOT, "neoantigen_pipeline")
DEMO_SRC = os.path.join(ROOT, "demo_out")
NAME = "neoantigen-selection"

# Demo artefacts worth shipping. candidates.csv is excluded on size.
DEMO_KEEP = [
    "REPORT.md", "audit_manifest.json", "selection_qc.json",
    "gate_waterfall.csv", "ranked.csv", "selected.csv", "coverage.csv",
    "peptides_skipped.csv", "minigenes.csv", "junction_scan.csv",
    "construct.fasta",
    "benchmark_metrics_allele_matched.csv",
    "benchmark_metrics_presentation_controlled.csv",
    "benchmark_metrics_by_stratum.csv",
    "tesla_metrics.csv", "tesla_per_patient.csv",
    "summary.png", "junctions.png", "benchmark.png",
]

REQUIREMENTS = """\
# Required
pandas>=2.0
numpy>=1.24
requests>=2.28

# Optional: figures in the demo report
matplotlib>=3.7

# Optional: local MHC-I prediction instead of the IEDB cloud REST API.
# Removes the network round-trip for large runs. Not required.
# mhcflurry>=2.1
"""


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def build(outdir: str = "dist", version: str = None) -> str:
    sys.path.insert(0, ROOT)
    from neoantigen_pipeline import __version__
    version = version or __version__

    stage = os.path.join(outdir, f"{NAME}-skill")
    if os.path.exists(stage):
        shutil.rmtree(stage)
    os.makedirs(stage)

    # --- skill docs -------------------------------------------------------
    shutil.copy2(os.path.join(SKILL_SRC, "SKILL.md"), stage)
    shutil.copytree(os.path.join(SKILL_SRC, "reference"),
                    os.path.join(stage, "reference"))

    # --- code -------------------------------------------------------------
    scripts = os.path.join(stage, "scripts")
    os.makedirs(scripts)
    shutil.copy2(os.path.join(SKILL_SRC, "scripts", "neoantigen.py"), scripts)
    shutil.copytree(
        PKG_SRC, os.path.join(scripts, "neoantigen_pipeline"),
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "cache"))

    # --- bundled data -----------------------------------------------------
    # The proteome goes exactly where fetch.human_proteome() looks, so a fresh
    # install never has to page 20,400 UniProt entries over the network.
    cache = os.path.join(scripts, "neoantigen_pipeline", "data", "cache")
    os.makedirs(cache, exist_ok=True)
    proteome = os.path.join(PKG_SRC, "data", "cache", "uniprot_human_reviewed.fasta.gz")
    bundled = {}
    if os.path.exists(proteome):
        shutil.copy2(proteome, cache)
        bundled["uniprot_human_reviewed.fasta.gz"] = sha256(proteome)
    tesla = os.path.join(PKG_SRC, "data", "tesla_deepimmuno_public.csv")
    if os.path.exists(tesla):
        bundled["tesla_deepimmuno_public.csv"] = sha256(tesla)

    # --- finished demo ----------------------------------------------------
    demo = os.path.join(stage, "demo")
    os.makedirs(demo)
    shipped = []
    for name in DEMO_KEEP:
        src = os.path.join(DEMO_SRC, name)
        if os.path.exists(src):
            shutil.copy2(src, demo)
            shipped.append(name)

    # --- install notes + requirements ------------------------------------
    open(os.path.join(stage, "requirements.txt"), "w").write(REQUIREMENTS)
    open(os.path.join(stage, "INSTALL.md"), "w").write(install_md(version, shipped))
    manifest = {
        "skill": NAME,
        "version": version,
        "built_at": datetime.now(timezone.utc).isoformat(),
        "bundled_data": bundled,
        "demo_files": shipped,
        "excluded": {
            "prediction cache": "~46 MB, specific to one patient and HLA set",
            "candidates.csv": "4.3 MB of intermediate peptide x allele rows",
        },
        "research_use_only": True,
    }
    json.dump(manifest, open(os.path.join(stage, "package_manifest.json"), "w"), indent=2)

    # --- zip --------------------------------------------------------------
    zip_path = os.path.join(outdir, f"{NAME}-skill-v{version}.zip")
    if os.path.exists(zip_path):
        os.remove(zip_path)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for base, _dirs, files in os.walk(stage):
            for f in sorted(files):
                full = os.path.join(base, f)
                z.write(full, os.path.join(NAME, os.path.relpath(full, stage)))
    return zip_path


def install_md(version: str, demo_files) -> str:
    return f"""\
# neoantigen-selection — install and run

Version {version}. Research use only. Not a clinical tool, and not a
reproduction of any company's proprietary neoantigen score or construct rules.

## Install

Unzip into your skills directory:

```bash
# for one project
unzip {NAME}-skill-v{version}.zip -d .claude/skills/
# or for every project
unzip {NAME}-skill-v{version}.zip -d ~/.claude/skills/
```

You should end up with `.claude/skills/{NAME}/SKILL.md`. Claude Code picks the
skill up on the next session; ask for it by name with `/{NAME}` or just describe
a neoantigen-selection task.

Then install the Python dependencies:

```bash
pip install -r .claude/skills/{NAME}/requirements.txt
```

## Check it works — offline, about two seconds

```bash
cd .claude/skills/{NAME}
python scripts/neoantigen.py selftest
```

43 checks over peptide tiling, wild-type pairing, the self k-mer index, junction
enumeration, the codon optimizer's translation invariant, the ordering
optimizer, and the metrics. No network.

## Run the public demo

```bash
python scripts/neoantigen.py demo --out demo_out --benchmark --tesla
```

Real TCGA-SKCM melanoma variants (cBioPortal open API), real NetMHCpan-4.1 EL
predictions (IEDB cloud REST), real IEDB and TESLA ground truth. First run takes
roughly 40 minutes, almost all of it waiting on cloud predictions; everything is
cached on disk afterwards, so a re-run is close to instant.

A finished run of exactly this command is already in `demo/` — read
`demo/REPORT.md` first if you would rather see the output than wait for it.

## Run on a patient

```bash
python scripts/neoantigen.py run \\
    --maf patient.maf \\
    --patient PT-014 \\
    --hla HLA-A*02:01 HLA-A*24:02 HLA-B*07:02 HLA-B*44:02 HLA-C*05:01 HLA-C*07:02 \\
    --expression tumor_tpm.csv \\
    --purity 0.62 \\
    --out PT-014_out
```

- `--maf` — a somatic MAF from your own tumor/normal caller. This skill does not
  call variants; pair it with a somatic caller upstream.
- `--hla` — **four-digit class-I typing from the NORMAL sample**, from
  OptiType / xHLA / HLA-HD. Never type from the tumor: LOH at the HLA locus
  corrupts the call. A wrong allele invalidates the whole ranking silently, and
  it is the single largest error source in the pipeline.
- `--expression` — a CSV with `entrez,tpm` from the patient's tumor RNA-seq.
  Omit it and the expression gate cannot run; the tool warns loudly, because
  ranking an unexpressed gene highly is the failure mode this pipeline exists to
  prevent.
- `--purity` — from ABSOLUTE / FACETS / Sequenza / PureCN. Without it the CCF
  column is an estimate, not a measurement, and the report says so.

Outputs: `ranked.csv`, `selected.csv` (one row per payload slot with
`why_selected`), `coverage.csv`, `minigenes.csv`, `junction_scan.csv`,
`construct.fasta`, `audit_manifest.json` and a labelled `REPORT.md`.

## What ships in the box

- `SKILL.md` — the skill itself.
- `reference/scoring.md` — every feature, its formula, its citation, how it fails.
- `reference/workflow.md` — the public workflow this mirrors; what is open, what
  is proprietary, and what you must supply for a real patient.
- `reference/benchmark.md` — both benchmarks, the TESLA results, and why the
  IEDB-derived numbers are lower bounds.
- `reference/comparison.md` — head-to-head with another neoantigen package.
- `scripts/neoantigen_pipeline/` — the code.
- `scripts/neoantigen_pipeline/data/cache/uniprot_human_reviewed.fasta.gz` — the
  UniProt reviewed human proteome (20,338 canonical sequences), bundled so the
  first run does not have to page it over the network. UniProt data is
  CC BY 4.0.
- `scripts/neoantigen_pipeline/data/tesla_deepimmuno_public.csv` — the public
  TESLA mirror, 522 peptide-HLA pairs with real T-cell assay labels.
- `demo/` — a finished run ({len(demo_files)} files), including the figures.

## Read this before quoting a number

The composite score is **not validated**. On the TESLA mirror, with real assay
labels and real negatives, NetMHCpan-4.1 EL %rank alone scores AP 0.207 and
recovers 31 of 35 experimentally immunogenic peptides in a 34-slot budget; this
package's composite scores 0.163 and 26/35. The extra peptide-intrinsic
features do not beat the binding predictor, which is why the default weights are
presentation-dominant and `config.PRESENTATION_ONLY` ships as a preset.

Every immunogenicity number this skill produces is an unverified prediction. The
`REPORT.md` it writes labels each claim `[computed]`, `[assumed]` or
`[unverified]`, and `audit_manifest.json` records the input hashes, weights,
predictor and evidence level behind the run.

## Network

Presentation prediction uses the IEDB cloud REST API
(`tools-cluster-interface.iedb.org`); the demo also reaches `www.cbioportal.org`
and `query-api.iedb.org`. Install `mhcflurry` and the MHC-I backend runs
locally instead. All prediction results are cached per batch on disk.
"""


if __name__ == "__main__":
    out = build()
    print(f"{out}  ({os.path.getsize(out) / 1e6:.1f} MB)")
