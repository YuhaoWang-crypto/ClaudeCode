# MHC-I immunogenicity report — E. coli L-asparaginase II (ANSII, P00805)

A real target run end-to-end through the **integrated** stack (this package's
schema + the `immunogenicity-multimodel` skill's IEDB backend). All numbers
below were produced live via IEDB cloud NetMHCpan — **no token, no local
license.**

## Target & panel

- **Protein**: *E. coli* L-asparaginase II, mature chain (326 aa; 22-aa signal
  peptide removed). A first-line ALL chemotherapy biologic and a **textbook
  de-immunization target** — anti-asparaginase antibodies drive hypersensitivity
  and silent inactivation in the clinic.
- **HLA-I panel** (broad population coverage): A\*02:01, A\*01:01, A\*03:01,
  A\*24:02, B\*07:02, B\*08:01.
- **Model**: IEDB NetMHCpan_el, 9-mers. 1908 peptide×allele predictions.

## Result: epitope load

| call (%Rank) | count |
|---|---|
| strong binders (≤0.5) | **23** |
| weak binders (≤2.0) | 53 |

Top promiscuous / strongest MHC-I epitopes (`results_ansb/mhc1_consensus.csv`):

| peptide | alleles | best %Rank |
|---------|---------|-----------|
| SADGPFNLY | A\*01:01 | 0.01 |
| LYKSVFDTL | A\*24:02, B\*08:01 | 0.02 |
| NPQKARVLL | B\*07:02, B\*08:01 | 0.03 |
| NLVNAVPQL | A\*02:01, B\*08:01 | 0.05 |
| SVNYGPLGY | A\*01:01, A\*03:01 | 0.12 |

Per-residue landscape: `figures/asparaginase_mhc1_landscape.png` (bright bands =
epitope hotspots shared across alleles → priority de-immunization regions).

## De-immunization demo (sequence-level, live)

Using MHC-I %Rank as the oracle (valid for AA sequences, unlike EDEN), the
strong A\*02:01 epitope **NLVNAVPQL** (%Rank 0.05) was probed by single
substitutions at its anchor positions P2 and P9
(`results_ansb/deimmunize_NLVNAVPQL.csv`):

| variant | peptide | %Rank | call |
|---------|---------|-------|------|
| WT | NLVNAVPQL | 0.05 | **strong** |
| L9P | NLVNAVPQP | 3.6 | none |
| L2G | NGVNAVPQL | 7.4 | none |
| L9R | NLVNAVPQR | 7.6 | none |
| … | … | … | … |

**All 10** anchor substitutions abolish the strong binder (%Rank 0.05 → 3.6–22).
That is the T-cell-epitope removal step. To keep the protein **folded and
functional** while making those changes, feed these positions to the
structure-aware redesigner instead of substituting blindly:

- `pipeline/deimmunize.py` — ProteinMPNN loop: fix all residues except the
  epitope positions, redesign, re-score each variant with this same MHC-I
  oracle, keep variants that lose the epitope while preserving fold. (ProteinMPNN
  needs a PDB + GPU; wired as the official CLI.)

## How this was run (integrated entry point)

```python
from pipeline import integrate, aggregate
alleles = ["HLA-A*02:01","HLA-A*01:01","HLA-A*03:01",
           "HLA-A*24:02","HLA-B*07:02","HLA-B*08:01"]
res  = integrate.run_iedb_mhci("examples/asparaginase_mature.fasta", alleles)
cons = aggregate.consensus([res], top_n=15)     # same machinery as native adapters
```

`integrate.run_immunogenn(...)` adds the MHC-II / population-immunogenicity layer
(and ImmunoGeNN's own single-point de-immunization scan) — it needs a
`BIOLIB_TOKEN` in this sandbox, since anonymous BioLib jobs do not complete here.
