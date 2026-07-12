# Adding a new pathway

The pipeline is pathway-agnostic. Adding one is a single `Pathway` entry in
`psmos/pathways.py` (plus, if it introduces new organisms, `Species` entries).
Everything downstream — ortholog search, gate, CDS, Compara, AlphaGenome, Evo2,
scoring — reuses the same code.

## Recipe

1. **Pick the component families** and mark the **gate families** — the
   component(s) whose absence means "no canonical pathway output". For a
   receptor→TF cascade the gate is usually {ligand-receiver, DNA-binding TF}
   (Notch: receptor + CSL; Hippo: YAP/TAZ + TEAD).
2. **List human paralogues per family** (drives the redundancy prior; Compara
   makes it computed later).
3. **Write the ortholog seed** for the gate families: per-species gene symbols,
   as a fallback list (first hit wins). Use each species' *own* historical gene
   name (fly `N`/`Su(H)`, worm `lin-12`/`lag-1`). For true negative controls,
   query the human symbol so the search genuinely returns nothing.
4. **Add any new species** to `SPECIES` with curated `tract`/`thru` priors and
   the NCBI taxon id. If it's on Ensembl, also add it to
   `cds.ENSEMBL_SPECIES` (and it will get CDS + Compara + Evo2); if not, it's an
   honest CDS/Evo2 gap.
5. **Curated baseline**: for the PSMOS 6-layer path, add the species' G/D/N/R/E/X
   priors + arch-role fits (see `CURATED_HIPPO` in `scoring_psmos.py`). For the
   lighter Notch path, add a `CURATED_<PATHWAY>` row.
6. **Run** `python3 -m psmos.run_all <Pathway>`.

See `assets/pathway_template.py` for a fill-in-the-blanks `Pathway`.

## Worked example already in the package

- **Wnt/β-catenin (regeneration ↔ fibrosis)** — `WNT` in `pathways.py`,
  `CURATED_WNT` in `scoring_psmos.py`, `Wnt` in `build_psmos_dashboard.PATHWAY_UI`.
  Gate = β-catenin (CTNNB1) + TCF/LEF; panel = planaria, zebrafish, axolotl,
  mouse (liver), human (IPF), naked mole-rat (fibrosis-resistant), fly (arm/pan),
  **yeast (live negative control — zero hits → disqualified)**. All three
  computed layers run live. Copy this as the pattern for a new pathway.

## Good next candidates (from the framework's own case list)

- **cGAS–STING innate immunity** — drug-target weighting; the value is
  comparing ligand recognition / STING agonist pharmacology across species
  (human vs mouse STING pockets differ), flagging model–drug mismatch before
  animal work.
- **Insulin–mTOR–FOXO / longevity** — yeast, worm, fly, mouse, naked mole-rat,
  killifish. Separates high-throughput lifespan-screen organisms from
  mammalian-metabolism translation, and flags long-lived species with genuine
  network rewiring vs simple conservation.

## Gotchas when adding a pathway

- Keep the pathway's species panel self-consistent: the Notch scorer iterates
  `CURATED_NOTCH` keys (its own panel), not all of `SPECIES` — so adding
  Hippo-only species doesn't leak into Notch. Do the same for a new pathway.
- Gate genes with divergent names across clades are the usual failure: give a
  fallback list and verify the resolved accession is the canonical one
  (length ≈ human), not a fusion isoform.
- If AlphaGenome coverage matters, remember it is human/mouse only — non-model
  panel members keep curated R by design.
