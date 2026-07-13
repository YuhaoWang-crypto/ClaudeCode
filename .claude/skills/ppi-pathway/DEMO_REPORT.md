# ppi-pathway skill — demo run report

**Date:** 2026-07-13 · **Environment:** Claude Code cloud container (Linux,
Python 3.11) · **Command:** `python -m ppi_pathway.demo` (+ `--with-string`)

> ⚠️ **Scope of this report.** This is a **validation / smoke test of the skill
> pipeline**, not a scientific analysis of a real dataset. The input is a fixed,
> deliberately *well-understood* gene set (16 DNA-damage / cell-cycle genes)
> chosen so that a correct pipeline must reproduce known biology. "Conclusions"
> below are therefore about **whether the tool works**, plus the example biology
> it recovers. To get research conclusions, run the skill on **your own** omics
> hit list.

Input gene set (n=16): `TP53 BRCA1 BRCA2 EGFR MYC CDK2 CDK4 CDK6 RB1 MDM2 ATM
ATR CHEK1 CHEK2 CCNE1 E2F1`

---

## 1. Engineering conclusions — does the skill work? ✅

**All five data sources are reachable from the cloud and produce correct output.
Demo result: 6/6 checks pass (full run), 4/4 (online-only).**

| Source | Access | Result | Status |
|---|---|---|---|
| STRING REST API | online | enrichment + coherence returned | ✅ works |
| humanPPI (Cong lab) | download ~14 MB | 12,268-node graph built | ✅ works |
| BioGRID | download (738 MB tab3) | 957,872 physical edges parsed | ✅ works |
| STRING bulk files | download ~105 MB | 16,201-node / 236,930-edge graph (score ≥700) | ✅ works |
| InterPro REST API | online | 9 entries for TP53 | ✅ works |
| Reactome REST API | online | **Cloudflare-blocked** → auto-fallback to STRING RCTM | ⚠️ blocked, handled |

Notable engineering outcomes:

- **humanPPI TLS chain fixed the compliant way.** `conglab.swmed.edu` ships an
  incomplete certificate chain (omits the *InCommon RSA Server CA 2*
  intermediate). The skill fetches that one intermediate and adds it to the
  trust store — **verification stays on, nothing disabled** — and the download
  then succeeds.
- **Reactome Cloudflare block is detected, not silently failed.** The skill
  raises `ReactomeBlocked` and transparently falls back to STRING's `RCTM`
  (Reactome) enrichment category, which is not blocked — so Reactome pathways
  are still delivered.
- **Raw data is never committed.** ~1 GB of downloads live in a git-ignored
  cache and are re-fetchable on demand — correct for an ephemeral container.

---

## 2. Example biology recovered (on the demo gene set)

These are **sanity-check** conclusions: the pipeline recovers textbook biology of
this gene set, which is exactly what a working tool should do.

### 2a. The gene set is a coherent module ✅
STRING PPI-enrichment: **113 observed interactions vs ≈39 expected** at random →
**p ≈ 0** (below machine precision). The genes interact ~2.9× more than chance,
so downstream enrichment is trustworthy (not a random grab-bag).

### 2b. Pathway enrichment points at DNA-damage response & cell cycle ✅
Top STRING enrichment (FDR-corrected):

| Category | FDR | Pathway |
|---|---|---|
| WikiPathways | 1.6e-29 | miRNA regulation of DNA damage response |
| WikiPathways | 1.6e-29 | DNA damage response |
| WikiPathways | 3.7e-26 | Integrated cancer pathway |
| WikiPathways | 3.5e-24 | Cell cycle |
| KEGG | — | Cell cycle |

Reactome (via STRING/RCTM fallback): **Cell Cycle** (3.0e-18), **Cellular
Senescence** (1.2e-12), **Regulation of TP53 Activity through Phosphorylation**
(1.6e-12), **Defective binding of RB1 mutants to E2F1** (3.2e-12).

→ Correct: this *is* the DNA-damage/cell-cycle checkpoint machinery.

### 2c. Network hubs are the expected control points ✅
humanPPI module around the seeds (13/16 seeds present, 254 nodes / 442 edges after
1 neighbour shell). Top hubs by degree/betweenness:

| Gene | Degree | Betweenness | Role |
|---|---|---|---|
| EGFR | 98 | 0.591 | seed — dominant hub/bottleneck |
| MYC | 26 | 0.158 | seed |
| MDM2 | 24 | 0.145 | seed (TP53 regulator) |
| CDK2 | 24 | 0.119 | seed |
| BRCA1 | 21 | 0.098 | seed |
| CDK1 | 15 | 0.020 | **connector** (pulled in, not a seed) |

→ EGFR and the CDK/MDM2/BRCA1 axis are the module's control points — plausible
drug-target candidates. CDK1 surfaced as a connector the input list didn't
contain, illustrating the value of the `--expand` neighbour shell.

### 2d. Domain annotation explains the hubs ✅
InterPro on TP53 (P04637) returns 9 entries: *p53 tumour suppressor family*, *p53
DNA-binding domain*, *tetramerisation domain*, *transactivation domain*, with GO
terms `DNA binding`, `DNA-binding transcription factor activity`. → correct
mechanistic annotation.

### 2e. Multi-source consensus interactome ✅
Overlaying all three networks (`--with-string`):

| Support level | Edges |
|---|---|
| STRING (score ≥700) | 236,930 |
| BioGRID (physical) | 957,872 |
| humanPPI (precision 80%) | 29,191 |
| **Overlay total** | **1,142,496** |
| supported by ≥2 sources | 65,349 |
| **supported by all 3 sources** | **8,074** |

→ **8,074 interactions are backed by predicted (STRING) + curated (BioGRID) +
structural (humanPPI) evidence simultaneously** — a conservative,
high-confidence human interactome backbone. This tri-evidence intersection is
the main analytical payoff of combining the three resources.

---

## 3. Bottom line

- **Engineering:** the skill is production-ready. All requested platforms
  (STRING, BioGRID, humanPPI, InterPro, Reactome) download/run from the cloud;
  the one blocked service (Reactome) degrades gracefully; the demo is green.
- **Biology (demo set):** the pipeline reproduces the known DNA-damage/cell-cycle
  network — coherent module, correct pathways, sensible hubs, right domains,
  triple-supported consensus edges. That is *validation*, not discovery.
- **Next step:** replace the demo list with **your** differentially-expressed /
  mutated / screen-hit genes and re-run `enrich` → `subnet` → `annotate`. Any
  conclusions from *that* run are real findings; this report only certifies the
  instrument reads true.

*Reproduce:* `cd .claude/skills/ppi-pathway/scripts && python -m ppi_pathway.demo --with-string`
