# immunogenicity-multimodel

Multi-model protein immunogenicity screening + deimmunization.

## Contents
- `SKILL.md`         — skill definition (loads into Claude Science)
- `kernel.py`        — auto-loaded helper functions
- `reference_scripts/` — standalone backend modules used to build the skill:
  - `mhc_iedb_backend.py`      — IEDB cloud NetMHCpan (MHC-I), no install
  - `mhc_standalone_backend.py`— local NetMHCpan/IIpan via mhctools (license-gated)
  - `merge_layer.py`           — one-protein-in multi-backend merge
  - `dtu_orchestrator.py`      — parallel BioLib DTU app runner

## Install as a Claude Science skill
Place the `immunogenicity-multimodel/` folder (SKILL.md + kernel.py) in your
skills directory, or re-create via host.skills.edit()/publish().

## Requirements
pip: pybiolib, requests, pandas, mhctools
Network allowlist: biolib.com, tools-cluster-interface.iedb.org

## Quick use (inside Claude Science)
    load skill immunogenicity-multimodel
    merged = run_multimodel("protein.fasta", alleles_mhci=["HLA-A*02:01"])
    models, M = per_position_matrix(merged, seq_len=len(seq))   # -> heatmap
    outdir, _ = run_immunogenn("protein.fasta", mode="deimmunize")  # lower immunogenicity
