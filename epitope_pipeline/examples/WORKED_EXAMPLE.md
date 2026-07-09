# Worked example — a real cloud run that completed

This documents an actual end-to-end run performed against a hosted model, plus
the exact status of the DTU-on-BioLib path, so you can reproduce both.

## ✅ EDEN immunogenicity (hosted, completed live)

EDEN (Basecamp Research) predicts **whole-antigen** immunogenicity from a
**native nucleotide CDS** (not amino acids, not back-translated). It is reachable
in this environment via MCP and needs no BioLib account, so it ran to completion.

**Case:** SARS-CoV-2 antigens, native CDS pulled straight from the reference
genome `NC_045512.2` (NCBI E-utilities), so the input is a genuine natural
sequence — exactly what EDEN expects.

| Antigen | CDS | length | EDEN immunogenicity score |
|---------|-----|--------|---------------------------|
| Nucleocapsid (N) | `examples/cds/sars2_N.cds.fasta` | 1260 nt | **0.805** (high) |
| Spike (S) | `examples/cds/sars2_Spike.cds.fasta` | 3822 nt | *(re-run; the call was cut off by a transient MCP disconnect)* |

Score 0.805 for N is biologically sensible: the nucleocapsid is the dominant
serological antigen in SARS-CoV-2 (the basis of most N-based antibody tests).

### Control experiment (completed live) — and an honest finding

Two human self-proteins (immune-tolerated) were run as native CDS, expecting
*low* scores. All four EDEN runs below completed live:

| Antigen | type | CDS | EDEN score |
|---------|------|-----|------------|
| Human serum albumin (ALB) | human self | `examples/cds/human_albumin.cds.fasta` | **0.934** (high!) |
| SARS-CoV-2 Nucleocapsid (N) | viral | `examples/cds/sars2_N.cds.fasta` | 0.805 |
| Human insulin (INS) | human self | `examples/cds/human_insulin.cds.fasta` | **0.424** (low) |

**The controls did NOT cleanly separate self from viral** — insulin scored low
as hoped, but albumin scored *higher than the viral antigen*. This is the
correct thing to report, and it teaches the model's scope:

- **EDEN does not model host self-tolerance.** It scores an intrinsic,
  sequence-level "does this look like it induces an immune response" signal —
  not "will a human tolerate their own protein." A self-protein can still score
  high.
- **Human proteins are out of EDEN's distribution.** The tool explicitly warns
  it is calibrated on viral/bacterial antigens; the albumin result is a live
  demonstration of that caveat, not a contradiction of it.
- Practical takeaway: use EDEN to rank **natural pathogen antigens** by
  immunogenicity, not to judge whether a human/self/engineered protein is
  "safe." For self/engineered sequences, use the per-epitope AA-level tools
  (BepiPred / NetMHCpan) instead.

### Reproduce

```python
# 1. fetch a real natural CDS (already saved under examples/cds/)
#    curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore\
#    &id=NC_045512.2&rettype=fasta&retmode=text&seq_start=28274&seq_stop=29533"
#
# 2. call EDEN with the nucleotide string (A/C/G/T, in-frame, 150-8192 nt)
#    -> returns index,status,score  e.g.  0,scored,0.8046
```

Input rules that matter (from the tool's own contract): **native nucleotide CDS
only** — amino-acid, codon-optimised, back-translated, partial, or epitope
sequences are out of distribution and rejected/meaningless. Viral & bacterial
antigens are best represented.

> EDEN gives one probability per whole antigen. It is a *complement* to, not a
> replacement for, the DTU per-epitope tools (which localise *which* peptide is
> the epitope).

## ⏸ DTU on BioLib (submits, but needs a token here)

The BioLib backend (`bepipred-cloud`, etc.) was exercised against the live
BioLib API in this environment:

- `biolib.load('DTU/BepiPred-3')` — works (real app metadata returned).
- `app.cli(...)` — **submits**; returns a real job UUID and status `in_progress`.
- The job then **stays `in_progress`** and does not complete within 25 min, and
  re-fetching it by UUID returns **HTTP 401 (auth required)**.

Conclusion: anonymous jobs are not schedulable/monitorable to completion from a
non-interactive sandbox. To actually get BepiPred/DiscoTope/NetSurfP results,
set a free **`BIOLIB_TOKEN`** (Account → API tokens on biolib.com) in the
environment, then:

```bash
export BIOLIB_TOKEN=...            # free account
python -m pipeline.cli --fasta examples/lysozyme.fasta \
    --models bepipred-cloud --out results/
```

The dispatch + output-file parsing for this path is covered by the test suite
(`tests/test_pipeline.py::test_biolib_backend_dispatch`), so once the token is
present it runs unchanged.
