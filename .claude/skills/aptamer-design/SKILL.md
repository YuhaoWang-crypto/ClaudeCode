---
name: aptamer-design
description: >-
  Design and in-silico rank DNA/RNA aptamer candidates against a protein target.
  Use when the user wants to design aptamers or nucleic-acid binders for a protein,
  build a structure-informed SELEX starting library, or computationally screen/rank
  aptamer candidates. Runs a rational, structure-informed pipeline: fetch target
  sequence + domain/epitope analysis, design folded scaffolds (G-quadruplex, hairpin,
  pseudoknot) and a SELEX-biased library, co-fold each candidate against the target
  with Boltz-2.1 (via the Boltz_API MCP server), then rank by interface confidence
  weighted for the stated use case (diagnostic probe vs therapeutic antagonist).
  Enforces honesty labeling: computational prioritization is NOT experimental validation.
---

# Aptamer design (DNA / RNA) against a protein target

## What this skill does

Turns a protein target into a **prioritized shortlist of DNA/RNA aptamer candidates**
plus a SELEX starting library, using structure-informed rational design and Boltz-2.1
co-folding as an in-silico filter. It does **not** replace SELEX — it narrows the
starting sequence space and ranks candidates by structural-confidence proxies.

## Two modes
- **Quick** (default): rational design + Boltz-2.1 co-folding + scramble decoy + paralog
  counter-screen. Fast, runs from the Boltz_API MCP alone. Steps 0–6 below.
- **Full in-silico SELEX** (`references/t-selex-integration.md`): an iterative loop that
  adds the **T-SELEX** toolchain (ViennaRNA → RNAComposer → HDOCK + IntaRNA) as a **second,
  orthogonal 3D scorer** and evolves winners with the LM over several rounds. Rank by
  **HDOCK ⊕ Boltz consensus** (`scripts/consensus_rank.py`). Use when you want a real
  SELEX-style search, not a one-shot shortlist. Needs a provisioned Linux env (see that doc).

## The golden rule (state this to the user, every time)

Real high-affinity aptamers come from **experimental SELEX**. In-silico co-folding
(Boltz-2.1 / any predictor) gives a *relative* prioritization signal — **ipTM is not
Kd**, and nucleic-acid interface confidence is less validated than protein–protein.
Never present ranked candidates as validated binders. Always attach the caveats in
`references/metric-interpretation.md`.

## Prerequisites

- **Boltz_API MCP server** connected (tools `boltz_estimate_structure_and_binding`,
  `boltz_start_structure_and_binding`, `boltz_get_structure_and_binding_prediction`).
  Estimates are free; runs are cheap (~$0.025 / complex) but **paid** — always run an
  estimate and get the user's OK before submitting paid jobs.
- Web access (UniProt REST, RCSB) for target sequence + structure/epitope facts.
- Optional: PubMed / bioRxiv MCP to check for prior aptamers and the ligand-binding footprint.

## Workflow

### 0. Clarify the use case (this changes the ranking)
Ask (or infer) whether the aptamer is for a **diagnostic/imaging probe** (specificity +
stability dominate; ligand-site competition is often a *liability* — endogenous ligand
can block the probe) or a **therapeutic antagonist** (overlap with the ligand-binding
footprint is desired). Pick the matching weight profile in `scripts/rank_candidates.py`.

### 1. Target intake
- Fetch the canonical sequence: `https://rest.uniprot.org/uniprotkb/<ACC>.fasta`.
- Fetch features: `.txt` — signal peptide, mature chain, domains, GPI/TM, disulfides,
  and any ligand-binding / interface residues.
- Identify the **folded functional module** to target and **crop** it (keep disulfide
  partners together). Cropping keeps Boltz fast/cheap and focuses design pressure.
- If a complex structure exists (RCSB), extract the ligand-binding footprint (residues
  within ~4.5 Å of the partner). Note the target species vs your construct species —
  **map residue numbers by alignment** if they differ (a common silent error).
- Note electrostatics: a basic (Arg/Lys-rich) or heparin/GAG groove is a natural docking
  site for polyanionic aptamers — but electrostatic grooves also risk cross-reactivity
  with other polyanion-binding proteins (specificity caveat for diagnostics).

### 2. Candidate design
See `references/design-principles.md`. Produce a diverse set across scaffolds:
G-quadruplex (compact, nuclease-resistant — best for diagnostics), hairpin/stem-loop,
dual-hairpin, pseudoknot, plus a random-fold family. Provide both DNA and RNA versions.
Generate a **SELEX-biased library** (fixed primers + partial structured core + randomized loops).
Optionally pre-filter folds with ViennaRNA (DNA: Mathews-2004 params; RNA: Turner-2004).

### 3. In-silico co-folding (Boltz-2.1)
See `references/boltz-workflow.md`. For each candidate, co-fold [target crop + aptamer]
as a 2-chain complex (protein auto-MSA, aptamer as `dna`/`rna` entity). **Always include
a scrambled/decoy control** to establish an ipTM baseline. Record ipTM, pTM,
complex_iplddt, complex_ipde. There is **no affinity head for nucleic acids** — do not
report a Kd.

### 4. Ranking
Run `scripts/rank_candidates.py --use-case diagnostic|therapeutic candidates.json`.
It computes a transparent composite (interface confidence + fold robustness + chemistry
suitability + specificity margin) and prints tiers. Cross-campaign ipTM differences
< ~0.03 are within noise — treat as ties.

**Full-SELEX mode:** when you also have HDOCK scores (T-SELEX half), rank by consensus
instead: `scripts/consensus_rank.py candidates.json --decoy-iptm <baseline>`. It ranks by
HDOCK⊕Boltz rank-agreement (Borda), flags single-scorer/disagree/NON-SPECIFIC/below-decoy,
and sinks candidates that fail the specificity or decoy gate. RNAComposer is RNA-only, so
DNA candidates are Boltz-only (single-scorer) — that's expected.

### 5. Specificity counter-screen (decisive for diagnostics)
Re-run the top candidates against **paralog targets** (e.g. related receptor family
members) AND an **unrelated protein** under the same construct/MSA. Prefer candidates whose
on-target ipTM clearly exceeds off-target ipTM. Feed off-target ipTM back into the ranker.

**Critical (learned the hard way):** absolute ipTM/HDOCK scores are unreliable for small
protein × ssRNA — a promiscuous RNA can score *higher* on an unrelated protein than on its
target (both Boltz and HDOCK). Always include a paralog **and** an unrelated negative, and
gate with `scripts/specificity_gap.py` (PASS only if on-target beats its decoy AND every
off-target on every scorer). For targets with close paralogs, use the full
**specificity-first SELEX** workflow in `references/specificity-first-selex.md`
(MSA-divergent-epitope targeting + counter-SELEX + this calibrated gate).

### 6. Deliverable
Use `templates/report_template.md`. Include: target rationale, methods, ranked tables,
predicted-complex structure files, the honesty caveats, and an experimental validation
plan (synthesis + modifications, BLI/SPR/MST binding, competition assay if antagonist,
doped-SELEX affinity maturation, cell functional readout).

## Chemistry / modification defaults
- Diagnostics → **DNA first** (cheap, stable, easy to label with biotin/fluorophore);
  G-quadruplex DNA is intrinsically nuclease-resistant. Add 3'-inverted-dT.
- Therapeutics / RNA → 2'-F or 2'-OMe pyrimidines, 3'-inverted-dT, PEGylation for half-life.

## Honesty labeling
Mark every quantitative claim: ✅ measured/《from structure》 vs ⚠️ predicted/hypothesis.
ipTM, pTM, pLDDT are all ⚠️ predicted confidence, not affinity or binding proof.
