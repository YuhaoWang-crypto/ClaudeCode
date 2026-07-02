# Stimulus-responsive, cell-type-selective synthetic promoters

Design assets + a proto-language (Evo2 / AlphaGenome) scoring pipeline for
synthetic promoters that switch **on** in response to a defined stimulus, in a
defined cell type.

## Tooling reality (read first)

Install status from this session:

- ✅ **`proto-language` (v0.1.0) installed from PyPI** (`pip install proto-language`).
  The pipeline is written against its **real, verified API** — `Segment` /
  `Construct` / `Constraint`, `Evo2Generator`, `MCMCOptimizer`, `Program`, and
  `alphagenome_interval_track_constraint` (signatures read from the installed
  source, incl. its doctests).
- ❌ **`proto-tools` (the GPU/model execution layer) NOT installable here.** It is
  GitHub-only (`evo-design/proto-tools`, not on PyPI), and this sandbox's proxy
  only allows your own repo, so the clone 403s. It also needs a GPU + HF-gated
  weights. It provides the Evo2 / AlphaGenome / Enformer / Borzoi wrappers.
- ❌ **Modal / Evo2 / AlphaGenome endpoints** — no tool, credential, or endpoint
  is exposed to this session. Check your own apps with `modal app list`.

So the split is:

| Step | Where it runs | Needs proto-tools + GPU? |
|------|---------------|--------------------------|
| RE modules + cassette assembly (`build_constructs.py`, `dual_and_designs.py`) | anywhere (done here) | no |
| Evo2 fill + AlphaGenome cell-type scoring (`design_pipeline.py`) | **your Modal env** | yes |

`design_pipeline.py` imports the execution layer lazily: without `proto-tools`
it prints the concrete **design plan** (generator + exact AlphaGenome
`ontology_terms` / `contrastive_ontology_terms`); with it, it runs the optimiser.
On your Modal box:
`pip install "git+https://github.com/evo-design/proto-tools.git"` then
`export HF_TOKEN=...`.

## The design logic

**A synthetic inducible promoter = [ response-element module ] x N + [ minimal promoter ].**

- The **module** decides *what signal* turns it on — it's the DNA element bound
  by the TF at the end of the pathway (the table below).
- **Cell-type selectivity** comes from three places, all scored in AlphaGenome
  per cell line: (1) is the pathway TF even expressed/active in that cell,
  (2) the minimal promoter choice, (3) optional cell-restricted enhancer.
  Some modules are *intrinsically* cell-selective (e.g. ERE fires only in
  ER+ cells like MCF-7).

### Stimulus -> pathway TF -> DNA module (in `elements.py`)

| Stimulus | Terminal TF | Module | Core motif |
|----------|-------------|--------|------------|
| Hypoxia | HIF1A-ARNT | HRE | `RCGTG` |
| TNF / LPS / inflammation | RELA / NF-kB | kB element | `GGGRNWYYCC` |
| cAMP / forskolin / beta-adrenergic | CREB1 | CRE | `TGACGTCA` |
| Type I IFN | STAT1/2/IRF9 | ISRE | `RNGAAANNGAAACT` |
| Type II IFN (gamma) | STAT1 | GAS | `TTCNNNGAA` |
| Estrogen | ESR1 | ERE | `GGTCANNNTGACC` |
| Glucocorticoid | NR3C1/GR | GRE | `GGTACANNNTGTTCT` |
| Oxidative stress | NFE2L2/NRF2 | ARE | `RTGACNNNGCR` |
| ER stress (IRE1) | XBP1s | UPRE | `TGACGTGG` |
| ER stress (ATF6) | ATF6 | ERSE-I | `CCAATN9CCACG` |
| ER stress (PERK) | ATF4 | AARE | `ATTGCATCA` |

The sequences in `elements.py` are literature-grounded **seeds**, not final
answers — the pipeline optimises copy number, spacers, flanks and minimal
promoter, then ranks by Evo2 likelihood + AlphaGenome inducibility/specificity.

## Composite designs (the two you asked for)

**1. ISRE + GAS pan-interferon sensor** (`dual_and_designs.py` -> `COMPOSITE`):
one enhancer interleaving ISRE and GAS, so it fires on **type I IFN (ISGF3->ISRE)
OR type II IFN (STAT1->GAS)**. Built at 2x/3x/4x on the low-baseline E1b core.

**2. Stimulus-AND-cell-type gates** (analog AND from enhancer synergy):
`[stimulus RE]xN + [lineage RE]xM + minimal promoter`. Strong output needs BOTH
the signal AND the lineage TF. Worked examples shipped:
IFN×myeloid (SPI1), hypoxia×hepatocyte (HNF4A), oxidative×hepatocyte (HNF1),
cAMP×neuronal (E-box). Cell selectivity is then *verified* by AlphaGenome
contrastive scoring — `ontology_terms`=target vs `contrastive_ontology_terms`=
off-targets. For a **digital** AND, use a two-component relay (stimulus-driven
split transactivator whose halves are each cell-restricted) rather than one
composite enhancer.

## Files

- `elements.py` — stimulus modules, **lineage modules**, **composite specs**,
  minimal promoters, per-cell-line context + AlphaGenome ontology map.
- `build_constructs.py` — assemble the 11 single-stimulus cassettes (no GPU).
- `dual_and_designs.py` — build pan-IFN sensor + AND-gate cassettes (no GPU).
- `design_pipeline.py` — real proto-language optimiser: Evo2 generator over a
  designable spacer + AlphaGenome contrastive cell-type scoring (Modal).
- `designs/*.fasta` — clone-ready cassettes (NheI/AgeI/KpnI/EcoRI handles).

## Run order on Modal

```bash
pip install proto-language
pip install "git+https://github.com/evo-design/proto-tools.git"
export HF_TOKEN=...            # gated Evo2 / AlphaGenome
# stimulus-responsive, cell-type-selective design:
python design_pipeline.py --stimulus interferon_typeII --target THP1
# stimulus-AND-cell-type gate:
python design_pipeline.py --stimulus hypoxia --target HepG2 --lineage hepatocyte
```

> Confirm the AlphaGenome ontology terms in `elements.CELL_CONTEXTS` /
> `LINEAGE_ELEMENTS` against your AlphaGenome build's `output_metadata` before
> running — tissue-level UBERON/CL terms are given; cell-line EFO terms sharpen
> the contrast.

## Caveats

- Seed element sequences should be confirmed against your references before
  ordering; the pipeline's model scores are only as good as the tracks you map.
- ARE overlaps AP-1/TRE — keep the NQO1 3' flank for NRF2 selectivity.
- ERSE-I needs the exact `CCAAT-N9-CCACG` geometry; don't collapse the spacer.
- These are research designs, not validated constructs — validate in a reporter
  assay (e.g. luciferase +/- stimulus) before use.
