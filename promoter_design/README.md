# Stimulus-responsive, cell-type-selective synthetic promoters

Design assets + a proto-language (Evo2 / AlphaGenome) scoring pipeline for
synthetic promoters that switch **on** in response to a defined stimulus, in a
defined cell type.

## Tooling reality (read first)

This repo was built in a session with **no access to Modal, Evo2, or
AlphaGenome** (no GPU, no credentials, no such endpoints). So the split is:

| Step | Where it runs | Needs GPU/model? |
|------|---------------|------------------|
| Response-element modules + cassette assembly (`build_constructs.py`) | anywhere | no |
| Evo2 naturalness + AlphaGenome expression scoring (`design_pipeline.py`) | **your Modal env** | yes |

Verify your Modal apps yourself with `modal app list`. The pipeline's two
`score_*` functions are documented **adapters** — wire them to your Evo2 /
AlphaGenome endpoints (names vary by proto-tools version; check
`help(proto_tools)`).

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

## Files

- `elements.py` — modules, minimal promoters, per-cell-line context map.
- `build_constructs.py` — assemble FASTA cassettes (no GPU). Run:
  `python build_constructs.py --copies 4 --min_promoter E1b_TATA`
- `design_pipeline.py` — Evo2 + AlphaGenome scoring/optimisation loop (Modal).
- `designs/*.fasta` — pre-built cassettes, clone-ready (with NheI/AgeI/KpnI/EcoRI handles).

## Run order on Modal

```bash
pip install git+https://github.com/evo-design/proto-language.git
export HF_TOKEN=...            # gated Evo2 / AlphaGenome
# implement score_evo2 / score_alphagenome against your endpoints, then:
python design_pipeline.py     # optimises one promoter per (stimulus, cell)
```

## Caveats

- Seed element sequences should be confirmed against your references before
  ordering; the pipeline's model scores are only as good as the tracks you map.
- ARE overlaps AP-1/TRE — keep the NQO1 3' flank for NRF2 selectivity.
- ERSE-I needs the exact `CCAAT-N9-CCACG` geometry; don't collapse the spacer.
- These are research designs, not validated constructs — validate in a reporter
  assay (e.g. luciferase +/- stimulus) before use.
