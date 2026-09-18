# The measured limits, and the API mechanics behind them

Everything here was established by running it. Where a number came from a single
worked example rather than a general result, it says so.

## Class II length floor — measured

```
POST https://tools-cluster-interface.iedb.org/tools_api/mhcii/
sequence_text = ">x\nVHAAHAEIN\n",  allele = H2-IAd,  length = 9
-> "Input length should be between 11-30."
```

**11 residues is a hard floor on input.** An 8-mer or a bare 9-mer core cannot be
submitted at all. Separately, the requested scan `length` imposes its own floor:
asking for a 15-mer window refuses any input shorter than 15. These are two
different limits and it is easy to mistake one for the other.

**Score a panel with ONE window across every construct.** A `%Rank` is a
percentile against a background of peptides of that same length, so scoring a
17-mer at length 17 and a 30-mer at length 30 makes them incomparable and
manufactures a difference that is an artefact of the window. In the worked
example, scoring each construct at its own length made a 30-mer diblock look
8× worse than its 17-mer uniblock; scored at a common 15-mer window both return
`VHAAHAEIN` at `%Rank 0.17` — identical, which is the real result.

## Padding a short core is not a workaround — measured

Padding `VHAAHAEIN` with glycine to clear the floor, each scored at its own
length:

| padded length | sequence | core | %Rank |
|---|---|---|---|
| 11 | `GVHAAHAEING` | VHAAHAEIN | 0.08 |
| 13 | `GGVHAAHAEINGG` | VHAAHAEIN | 0.30 |
| 15 | `GGGVHAAHAEINGGG` | VHAAHAEIN | 0.24 |
| 17 | `GGGGVHAAHAEINGGGG` | VHAAHAEIN | 0.39 |

The core is stable but `%Rank` spans **4.9×** purely from how much padding was
added. Anyone padding a short peptide to get a number is importing an artefact of
that size. Report the parent construct's core instead, and say the bare core was
not scoreable.

## Non-natural residues — a refusal, not a caveat

D-amino acids, N-methylation, β-residues, cyclisation. Every MHC predictor in
service is trained on L-peptides over the 20 canonical residues; there is no
input channel for stereochemistry. Submitting the sequence returns the score of
the all-L twin.

Do not score these and quote the number with a footnote. Exclude them, and state
which comparisons that removes from reach — usually the study's central one,
since a D-substituted test peptide normally exists precisely because the sponsor
expects stereochemistry to matter.

## IEDB class I mechanics

```
POST https://tools-cluster-interface.iedb.org/tools_api/mhci/
method = netmhcpan_el | netmhcpan_ba
sequence_text = multi-record FASTA, one peptide per record
allele = "HLA-A*32:01,HLA-B*57:01"
length = "9,9"                     # one entry per allele, parallel lists
```

Submitting each peptide as its own FASTA record returns exactly one row per
(peptide, allele) with no offset bookkeeping — different from the class II case,
where concatenation into pseudo-proteins is worth the trouble.

Both heads return `percentile_rank`. `netmhcpan_ba` **also** returns `ic50`,
which is the quantity a cell-free binding assay measures and therefore the one
to compare a wet-lab result against.

**The IC50 < 500 nM standard rejects real epitopes.** `FSPEVIPMF`, a documented
HLA-B*57:01 epitope, scores `EL %Rank 0.47` (strong) but `BA IC50 3017 nM` —
six-fold outside the classic affinity standard. A screen gated on IC50 < 500 nM
would have discarded it. Report both, and do not let the affinity standard veto
an eluted-ligand call.

Batch ~900 peptides per request. Wall-clock tracks request count more than
payload, and long runs need an append-as-you-go cache so a timeout does not cost
the whole run.

## IEDB query API for labelled data

```
https://query-api.iedb.org/mhc_search
  ?select=linear_sequence,mhc_allele_name,qualitative_measure,quantitative_measure,...
  &mhc_class=eq.I
  &mhc_allele_name=eq.HLA-B*57:01        # '*' raw — a percent-encoded '*' is rejected
  &order=structure_id&limit=1000&offset=0 # offset without order is refused
```

Column names differ from `tcell_search`; `assay_type_name` does not exist here.
Fetch one row with no filters and read the keys rather than guessing.

## The class-imbalance trap — measured

For the two alleles in the worked example:

| allele | records | positive 9-mers | negative 9-mers | positive fraction |
|---|---|---|---|---|
| HLA-A*32:01 | 9,021 | 4,044 | 165 | **96 %** |
| HLA-B*57:01 | 46,369 | 7,262 | 1,636 | 82 % |

Eluted ligands are positive by construction, so an IEDB-derived class I benchmark
is positive-enriched to a degree that makes PPV and NPV computed at its own
prevalence meaningless — they describe IEDB's collection policy.

Handle it by reporting prevalence-independent quantities (AUC, sensitivity,
specificity) from the data, then recomputing PPV/NPV at a **stated** prior and
sweeping it:

```
PPV = sens·p / (sens·p + (1-spec)·(1-p))
NPV = spec·(1-p) / (spec·(1-p) + (1-sens)·p)
```

Note also that A*32:01's AUC rests on only 165 negatives — report the bootstrap
CI, which will be wide, and do not present a data-poor allele's number with the
same confidence as a data-rich one.

## Allele frequencies for donor arithmetic

The IEDB population-coverage package ships a pickle with class I **and** class II
frequencies per population:

```python
sys.path.insert(0, f"{pkg}/deps/population-coverage-pickle")
sys.path.insert(0, pkg)
from population_coverage_pickle import population_coverage
freqs = dict(population_coverage["I"]["Europe"]["HLA-A"])   # {allele: frequency}
```

Loci available under class I: `HLA-A`, `HLA-B`, `HLA-C`. Class II is where the
DRB1 work in the sibling skill draws from, and the renormalisation trap
documented there (populations whose frequencies sum above 1) applies to class I
tables too — check the sum before trusting a coverage figure.

Carrier frequency is single-locus Hardy-Weinberg, `1 - (1-f)²`. Combining two
loci as independent is an approximation that linkage disequilibrium violates;
state the direction of the error for the specific pair. In the worked example,
B*57:01 in Europeans travels mostly on the 57.1 ancestral haplotype
(A*01:01-B*57:01-C*06:02-DRB1*07:01), which carries A*01:01 — so pairing it with
A*32:01 makes the independent product an **over**-estimate of double carriers,
and the real study is harder than the arithmetic suggests.
