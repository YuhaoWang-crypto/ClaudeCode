---
name: immunogenicity-multimodel
description: Screen a protein for T-cell and B-cell immunogenicity across multiple DTU/IEDB models in one call, and optimize a sequence to lower immunogenicity. Use when the task is to assess or reduce the immunogenicity/antigenicity of a protein or peptide, predict MHC-II population immunogenicity risk (pIRS), predict MHC-I binding epitopes (NetMHCpan %rank), compare self vs foreign proteins, deimmunize a biologic by scanning single-point variants, or produce a combined per-residue immunogenicity landscape from a FASTA. Backends: BioLib DTU/ImmunoGeNN (anonymous, no token), IEDB NetMHCpan cloud REST, and optional local NetMHCpan/NetMHCIIpan via mhctools.
---

# Immunogenicity multi-model screening & deimmunization

One protein FASTA in -> several immunogenicity models run in parallel -> one
unified per-peptide table + a combined per-residue landscape. Plus sequence
optimization to lower immunogenicity via ImmunoGeNN's deimmunize mode.

## Environment
Needs `biolib`, `requests`, `pandas` (and `mhctools` only for the optional
local NetMHCpan backend). Install into the analysis env if missing:
`manage_packages(mode="install", environment="<env>", packages=["pybiolib","requests","mhctools"], use_pip=True)`.
Network: `biolib.com` and `tools-cluster-interface.iedb.org` must be
allowlisted (request_network_access if blocked).

## Backends (verified)
| model | layer | scale | access |
|---|---|---|---|
| DTU/ImmunoGeNN | MHC-II population immunogenicity (pIRS, per 15mer) | 0-100 rank; >=83 = immunogenic | BioLib cloud, **anonymous** (no token) |
| NetMHCpan (netmhcpan_el) | MHC-I binding (%rank, per 9mer x allele) | strong<=0.5, weak<=2.0 | IEDB cloud REST |
| NetMHCpan/IIpan (standalone) | MHC-I/II binding | MHC-I strong<=0.5, weak<=2.0; MHC-II strong<=2.0, weak<=10.0 | local binary via mhctools (license-gated) |

NetMHCpan-family tools are NOT on BioLib — MHC-I binding goes through IEDB
cloud or a locally-installed binary.

## Kernel helpers (auto-loaded from kernel.py)
- `read_fasta(path)` -> {name: seq}
- `run_multimodel(fasta_path, alleles_mhci=(...), mhci_lengths=(9,), use=("immunogenn","iedb"))`
  -> unified DataFrame (columns: id, source, model, allele, length, peptide,
  position, affinity_nM, percentile_rank, score, is_binder)
- `run_immunogenn(fasta_path, mode="screen"|"deimmunize")` -> (output_dir, df).
  In `deimmunize` mode reads nothing back; the output_dir holds `SAVs.csv`
  (9600 single-point variants scored), `scores.csv` (per-variant DRB1_pIRS_sum),
  and `deimmunized_variants.fasta` (top variants). Rank variants by lowest
  `DRB1_pIRS_sum` vs the wild-type row.
- `run_iedb_mhci(sequences, alleles, lengths=(9,), method="netmhcpan_el")` -> unified DataFrame
- `per_position_matrix(merged, seq_len)` -> (models, matrix[len(models) x seq_len]),
  intensity 0-100 (higher = stronger binding / immunogenicity), for a heatmap
- `classify_binder(rank, mhc2=False)` -> "strong"/"weak"/"none"

## Workflow
1. **Screen**: `merged = run_multimodel("prot.fasta", alleles_mhci=["HLA-A*02:01","HLA-A*01:01"])`
2. **Landscape**: `models, M = per_position_matrix(merged, seq_len=len(seq))` then `imshow(M)`.
3. **Self vs foreign sanity check**: a human self-protein should give ~0
   ImmunoGeNN-immunogenic peptides (cores filtered against the human proteome);
   a foreign antigen gives many. Read per-peptide `pIRS_rank`, NOT the
   `scores.csv` `DRB1_pIRS_sum` (which is a pre-filter aggregate and stays
   non-zero even when every per-peptide pIRS is 0).
4. **Deimmunize**: `outdir, _ = run_immunogenn("prot.fasta", mode="deimmunize")`;
   read `scores.csv`, subtract each variant's `DRB1_pIRS_sum` from the
   wild-type row to get the reduction; each variant carries an ESM2
   log-likelihood so you keep only "natural" mutations.

## Notes
- ImmunoGeNN needs sequences >=15 residues.
- ProteinMPNN is NOT immunogenicity-aware: to use it for deimmunization, treat
  it as a structure-conditioned variant generator and score each variant with
  `run_immunogenn` — it needs a structure (PDB/AlphaFold) and a GPU.
- Only `netmhcpan_el` is verified against the IEDB endpoint; confirm other
  method names against the live IEDB method list before use.
