# Optional cross-checks with hosted models (MCP)

The DTU stand-alone tools are the backbone of this pipeline. When you also have
access to hosted models — e.g. through the MCP servers wired into this Claude
Code environment — you can run them as an **independent second opinion** and
compare against the DTU consensus. These are *not* DTU services and are kept out
of the core pipeline so it stays self-contained and license-clean.

## Available in this environment

| MCP tool | Use | Maps to |
|----------|-----|---------|
| `EDEN_by_Basecamp_Research.predict_immunogenicity` | sequence → immunogenicity score | orthogonal check on DTU MHC/CTL binders |
| `Boltz_API.start_structure_and_binding` | pMHC / peptide-target structure + binding | structural confirmation of a top epitope |

## Suggested workflow

1. Run the DTU pipeline (`python -m pipeline.cli ...`) → `consensus_epitopes.csv`.
2. Take the top N consensus peptides.
3. For each, call `predict_immunogenicity` on the peptide sequence.
4. Flag peptides where DTU consensus **and** EDEN agree → highest-confidence set.
5. (Optional) Fold the strongest peptide–HLA pair with Boltz for a structural view.

Because MCP tool availability is session-scoped (servers connect/disconnect),
this step is intentionally driven from the orchestration layer (a notebook or
the assistant) rather than hard-coded into the batch pipeline. The join key is
the `peptide` column already present in `consensus_epitopes.csv`, so merging an
EDEN score column back in is a one-line `pandas.merge`.
