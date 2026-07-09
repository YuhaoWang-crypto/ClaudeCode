# De-immunization by inverse folding (ProteinMPNN + epitope re-scoring)

Goal: lower a protein's predicted immunogenicity **without breaking its fold**,
by redesigning only the surface residues inside predicted epitopes.

```
antigen (seq + PDB)
   │
   ├─▶ 1. predict epitopes         BepiPred / NetMHCpan  →  epitope positions
   │
   ├─▶ 2. ProteinMPNN redesign     fix everything EXCEPT epitope residues
   │        (structure-preserving) →  N candidate sequences (same backbone)
   │
   ├─▶ 3. re-score each candidate  BepiPred / NetMHCpan  →  new epitope load
   │
   └─▶ 4. accept variants with     lower epitope load  AND  good ProteinMPNN
            score (foldability) AND enough sequence identity to wild type
```

Implemented in `pipeline/deimmunize.py`; tested end-to-end with stubs
(`tests/test_pipeline.py::test_deimmunize_loop`).

## The key methodological point (read this)

**EDEN cannot be the re-scoring oracle here.** EDEN takes a *native nucleotide
CDS* and explicitly rejects amino-acid / back-translated / synthetic sequences
as out-of-distribution. ProteinMPNN emits amino-acid sequences with no natural
CDS, so scoring them with EDEN is invalid. Step 3 therefore uses **amino-acid
level** predictors:

| Reduce which immunogenicity | Oracle | Metric (lower = better) |
|---|---|---|
| Antibody / B-cell | BepiPred-3.0 | mean epitope-region probability |
| T-cell (CD8) | NetMHCpan | count of peptides with %Rank ≤ 2 |
| T-cell (CD4) | NetMHCIIpan | count of peptides with %Rank ≤ 5 |

EDEN stays useful as a **whole-antigen sanity check on natural inputs**, not as
the design loop's scorer.

## What ProteinMPNN does and does not guarantee

- ✅ Sequences compatible with the given backbone (fold-preserving by design).
- ✅ Free, CPU-friendly, MIT-licensed (github.com/dauparas/ProteinMPNN).
- ❌ Does **not** guarantee function, stability, expression, or *actual* reduced
  immunogenicity. Every accepted variant is a hypothesis for wet-lab testing.
- Fix catalytic / binding / disulfide residues via `keep_fixed` so the loop
  never mutates them, even if they fall inside an epitope.

## Run it

```python
from pipeline.deimmunize import deimmunize, run_proteinmpnn_cli, parse_mpnn_fasta
from pipeline.deimmunize import fixed_positions_dict, positions_to_redesign

wt   = "MK...."                 # antigen amino-acid sequence
epi  = [34,35,36,37,38,39]      # epitope positions from BepiPred/NetMHCpan (1-indexed)

def redesign(design_pos):
    fp = fixed_positions_dict("myprot", "A", len(wt), design_pos)
    fa = run_proteinmpnn_cli("myprot.pdb", fp, "mpnn_out", num_seq=32)
    return parse_mpnn_fasta(fa)

def epitope_score(seq):
    # wrap BepiPred/NetMHCpan on `seq`; return mean epitope prob or SB count
    ...

result = deimmunize(wt, epi, redesign, epitope_score, keep_fixed=[50, 120])
print(result["best"])           # lowest-immunogenicity fold-preserving variant
```

### Getting the pieces in this environment

- **Epitope prediction** — `bepipred-cloud` (BioLib) or `netmhcpan` (local CLI),
  already wired in `pipeline/models.py`. BioLib cloud needs a free `BIOLIB_TOKEN`;
  MHC binders need the local academic package.
- **ProteinMPNN** — the official CLI (`protein_mpnn_run.py`). It is **not** on
  BioLib under any obvious id (checked: `DTU/ProteinMPNN`, `dauparas/ProteinMPNN`,
  etc. all 404). Clone the GitHub repo and put the script on PATH; it runs on CPU.
- **A structure** — supply a PDB, or fold the sequence first (AlphaFold /
  ESMFold / Boltz) to get the backbone ProteinMPNN needs.
```
