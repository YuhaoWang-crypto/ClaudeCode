# Target Gene Set: The Human Kinome

## Definition options

| Scope | Approx. genes | What it includes | When to use |
|-------|---------------|------------------|-------------|
| **Protein kinome (default)** | **518** | All eukaryotic protein kinases (ePKs) + atypical protein kinases (aPKs), per the Manning et al. 2002 classification | Standard kinase-function/dependency screens |
| Protein + lipid | ~538 | Adds PI3K family, PIP kinases, sphingosine/diacylglycerol kinases | If lipid-signaling nodes matter (e.g., PI3K/AKT axis) |
| Extended "kinome" | ~635 | Adds metabolic small-molecule kinases (sugar, nucleotide, adenylate kinases) | Metabolism-focused screens |
| Broad kinome + regulators | ~750–900 | Adds pseudokinases, kinase-associated subunits, phosphatases (if desired) | Comprehensive signaling-network screens |

**Recommendation:** target the **518 protein kinases** as the core, and add
the ~20 lipid kinases (total ~538). This is the most common and comparable
scope in the literature and keeps the library small and well-powered.

## Manning classification (protein kinome group counts)

| Group | Description | Approx. members |
|-------|-------------|-----------------|
| TK | Tyrosine kinases (RTKs + non-RTKs; EGFR, ABL1, SRC…) | ~90 |
| TKL | Tyrosine-kinase-like (RAF, IRAK, LRRK2…) | ~43 |
| STE | STE7/11/20 (MAP2K, MAP3K…) | ~47 |
| CK1 | Casein kinase 1 family | ~12 |
| AGC | PKA/PKG/PKC (AKT1-3, PKC, S6K, ROCK…) | ~63 |
| CAMK | Ca²⁺/calmodulin-dependent (CAMK, AMPK/PRKAA, CHK…) | ~74 |
| CMGC | CDK/MAPK/GSK3/CLK (CDK1-13, MAPK1/3, GSK3B…) | ~61 |
| RGC | Receptor guanylate cyclases | ~5 |
| Other | Aurora, PLK, NEK, IKK, WEE1, TTK… | ~83 |
| Atypical | ATM/ATR/mTOR/DNA-PK/PIKKs, ABC1, RIO, TAF1… | ~40 |

Totals are approximate; exact membership depends on annotation version.

## Authoritative sources for the exact gene list

Do **not** hand-curate the list — pull it from a maintained authority so gene
symbols map cleanly to Ensembl/RefSeq for guide design:

1. **HGNC gene group "Protein kinases"** (group ID 1258) — canonical, versioned
   HGNC symbols. REST: `https://rest.genenames.org/fetch/gene_group/1258`
2. **KinBase / Manning kinome** — `http://kinase.com/web/current/kinbase/`
   (the original 518 classification and phylogeny).
3. **KinMap / Coral** — interactive kinome tree, exportable gene lists
   (`http://www.kinhub.org/kinmap/`).
4. **UniProt keyword KW-0418 (Kinase)** + taxonomy 9606 for a superset
   including metabolic kinases.

### Fetch the HGNC protein-kinase set (reproducible)

```bash
# Requires curl + jq. Emits one HGNC symbol per line.
curl -s -H 'Accept: application/json' \
  https://rest.genenames.org/fetch/gene_group/1258 \
| jq -r '.response.docs[].symbol' | sort -u > data/kinome_symbols.txt
wc -l data/kinome_symbols.txt   # expect ~500+ symbols
```

Then map symbols → Ensembl transcript IDs (needed for constitutive-exon
targeting in guide design):

```bash
# Via Ensembl BioMart or pyensembl; example with pyensembl:
# pyensembl install --release 110 --species homo_sapiens
python3 - <<'PY'
from pyensembl import EnsemblRelease
data = EnsemblRelease(110)
for sym in open('data/kinome_symbols.txt').read().split():
    try:
        g = data.genes_by_name(sym)[0]
        print(sym, g.gene_id, g.contig, g.start, g.end, g.strand)
    except Exception:
        print(sym, "UNMAPPED")
PY
```

## Illustrative high-priority targets (sanity anchors)

Not the full list — a representative slice to confirm your mapping is sane and
to seed positive controls / known dependencies:

- **Cell cycle:** CDK1, CDK2, CDK4, CDK6, AURKA, AURKB, PLK1, WEE1, BUB1, TTK
- **Growth/survival:** EGFR, ERBB2, MET, ALK, KIT, AKT1, MTOR, PIK3CA, MAPK1, BRAF
- **DNA damage:** ATM, ATR, PRKDC (DNA-PK), CHEK1, CHEK2
- **Metabolism/stress:** PRKAA1 (AMPK), MAP2K1, MAPK14 (p38), GSK3B
- **Pan-essential positive controls (expected dropouts):** PLK1, CDK1, KIF11
  is not a kinase — use CDK1/PLK1/AURKB and pan-essential non-kinase controls
  from Hart CEGv2 for the essential-gene benchmark.

> Note on completeness: NGS files (FASTQ), full BAM/BED off-target enumerations,
> and any per-guide genome-wide off-target tables are large (>1 MB) and are
> **referenced by name**, not linked, in the top-level summary.
