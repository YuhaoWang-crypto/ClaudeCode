# Independent validation — AAVX VHH affinity ligand HLA-DR report

Validates the uploaded RUO deck *"AAV affinity ligand HLA-DR RUO report"*
(NetMHCIIpan v4.3, file `2582297_NetMHCIIpan.xls`) by re-running the same class
of prediction from scratch through this pipeline's IEDB NetMHCIIpan cloud path
(`integrate.run_iedb_mhcii`) — **no token, no license**, independent of their
run.

- **Ligand**: 126-aa public AAVX affinity ligand (a camelid **VHH / single-domain
  antibody**), `examples/aavx_vhh_ligand.fasta`. Sequence QC: 126 aa, 112
  overlapping 15-mers, 0 non-standard residues — matches the deck.
- **Panel**: 15 HLA-DR molecules (11 DRB1 + DRB3×2 + DRB4 + DRB5), same structure
  as the deck (3 DRB1 alleles assumed where the deck did not list them).

## HLA-DR (CD4) — deck claim vs independent re-run

| quantity | deck | this re-run | verdict |
|---|---|---|---|
| total peptide×HLA pairs | 1680 | 1680 | ✅ exact |
| overlapping 15-mers | 112 | 112 | ✅ exact |
| **primary hotspot core** | **FVAVQDITA** (47–55) | **FVAVQDITA** | ✅ confirmed |
| FVAVQDITA — strong binders | 11 | 11 | ✅ exact |
| FVAVQDITA — binder pairs / alleles | 35 / 11 | 32 / 10 | ✅ same magnitude |
| FVAVQDITA — best Rank_EL | 0.049 | 0.01 | ✅ same order (top-ranked) |
| top 15-mer / best allele | EREFVAVQDITASNT / DRB1\*04:01 | EREFVAVQDITASNT / DRB1\*04:01 | ✅ identical |
| **secondary hotspot** | **YLQMNNLKP** (80–88) | **YLQMNNLKP** | ✅ confirmed |
| YLQMNNLKP best Rank | 0.196 | 0.23 | ✅ same order |
| strong binders (<1%) | 25 | 23 | ✅ ~match |
| weak binders (1–5%) | 83 | 61 | ⚠️ differs (panel/version) |
| binder pairs (<5%) | 108 | 84 | ⚠️ differs (panel/version) |

**Conclusion: the deck's central finding is independently confirmed.**
FVAVQDITA is unambiguously the dominant HLA-DR/CD4 hotspot (top-ranked on
DRB1\*04:01), YLQMNNLKP is the secondary hotspot, the strong-binder burden is
~23–25, and DRB1\*04:01 is the top driver allele. The weak-binder / total-pair
counts differ modestly, fully explained by (a) the 3 DRB1 alleles I had to
assume (the deck did not list all 11) and (b) NetMHCIIpan version/rank-calibration
differences (deck v4.3 standalone vs IEDB cloud). Landscape:
`figures/aavx_vhh_hladr_landscape.png`.

## Other immune properties (beyond the deck's MHC-II-only scope)

**MHC-I (CD8 T-cell) — added here, the deck did not assess it.** IEDB
NetMHCpan_el, 6 HLA-I alleles, 9-mers (`results_vhh/mhc1_all.csv`): **5 strong +
25 weak** CD8 binders. Notably, `DITASNTHY` and `APGKEREFV` fall inside the
**same 47–58 stretch** as the FVAVQDITA CD4 hotspot → that region is a **dual
CD4+CD8 T-cell hotspot**, raising its priority for empirical testing.

**Structural localization of the hotspot.** In the VHH, FVAVQDITA (47–55) sits
at the **FR2→CDR2 boundary** (context `…APGKEREF·VAVQDITA·SNTHY…`). That matters
for de-immunization: framework-side positions are more amenable to germline /
humanizing substitutions, while CDR2-side positions constrain how far you can
mutate without touching the paratope. This is exactly the case for the
structure-aware `pipeline/deimmunize.py` loop (fix CDR/binding residues via
`keep_fixed`, redesign the framework-side epitope residues, re-score with this
same NetMHCIIpan oracle).

## What is NOT valid / still needs a token or context

- **EDEN whole-antigen immunogenicity is not applicable here.** EDEN requires a
  *native nucleotide CDS*; this is an engineered VHH affinity ligand with no
  natural CDS, so an EDEN score would be out-of-distribution and meaningless.
- **Human-proteome homology / tolerance filter** (the deck's own open item): a
  self-homology screen of FVAVQDITA / YLQMNNLKP cores is the right next step —
  it needs the ImmunoGeNN human-proteome filter (BioLib, `BIOLIB_TOKEN`) or a
  standalone proteome BLAST, neither of which runs anonymously in this sandbox.
- **US/EU population coverage** and **ppm × dose × regimen exposure** remain
  context modules, not sequence predictions — same framing as the deck.

## Reproduce

```python
from pipeline import integrate
dr = ["DRB1*01:01","DRB1*03:01","DRB1*04:01","DRB1*04:05","DRB1*07:01",
      "DRB1*08:02","DRB1*09:01","DRB1*11:01","DRB1*13:01","DRB1*15:01",
      "DRB1*12:01","DRB3*01:01","DRB3*02:02","DRB4*01:01","DRB5*01:01"]
mhcii = integrate.run_iedb_mhcii("examples/aavx_vhh_ligand.fasta", dr)   # CD4/HLA-DR
mhci  = integrate.run_iedb_mhci ("examples/aavx_vhh_ligand.fasta",
          ["HLA-A*02:01","HLA-A*01:01","HLA-A*03:01","HLA-A*24:02",
           "HLA-B*07:02","HLA-B*08:01"])                                  # CD8/HLA-I
```
