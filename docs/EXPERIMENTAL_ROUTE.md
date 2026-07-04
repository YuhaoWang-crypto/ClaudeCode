# Experimental route — alignment → distal candidates → His-library → metal screen

最稳妥的实验路线：先用 alignment 找 5–20 个远端候选区域，再围绕最高优先级的
hinge/界面 loop 做 His-pair / His-triplet 小文库，筛选 Zn²⁺ / Ni²⁺ / Co²⁺ 对
kcat、Km、kcat/Km 的影响。

## Stage 0 — structures / sequences (this repo)
- `python3 scripts/download_structures.py` fetches all 13 PDBs + FASTA.
- Pairs are defined in `data/targets.json`. Monomeric positive controls (GCK, PTP1B, AdK)
  are wired into the candidate finder; PFK/ATCase/GP need interface-aware extensions (see below).

## Stage 1 — alignment → 5–20 distal candidate regions (computational)
`python3 scripts/find_allosteric_candidates.py [SYS ...]`

Per system this produces:
- `results/<SYS>_candidates.csv` — per-residue Cα displacement + distance to active site + core flag.
- `results/<SYS>_regions.json` — ranked contiguous distal regions with His-pair suggestions.

**Selection logic** (tunable at the top of the script):
- `DISPLACEMENT_MIN = 2.0 Å` — residue counts as "moved".
- `DISTAL_MIN = 15.0 Å` — residue counts as "distal from active site" (allosteric, not orthosteric).
- Regions ranked by `mean_displacement × length` so real mobile domains beat single-residue termini.

**Prioritisation for engineering** (pick 3–6 of the 5–20 to take into the wet lab):
1. Distal + high mean displacement (state-coupled).
2. Sits at a **hinge or domain/subunit interface** (mechanically leveraged, not a floppy terminus).
3. Surface-accessible loop (mutable without destabilising the core).
4. Contains a geometrically plausible **His-pair/triplet** site (Cα–Cα 4.5–12 Å, see `his_pair_candidates`).

> Recommended first engineering targets from the current runs:
> - **AdK** — LID 117–163 and the LID/CORE hinge; NMP 31–70. Small protein, big motion → cleanest Zn-switch.
> - **GCK** — hinge/connecting region around 65–77 and 446–461 (activator-site adjacent).
> - **PTP1B** — WPD-loop edge (185–189) coupled to the distal 239–241 / 279–282 patches (drug-like allosteric site).

## Stage 2 — His-pair / His-triplet small library (molecular biology)
For each prioritised region:
- Turn the flagged Cα–Cα pairs/triples into **His substitutions** (a metal-binding His₂ or His₃ site).
  Verify side-chain geometry in a viewer/Rosetta before ordering — Cα–Cα spacing is a coarse filter;
  the Cβ vectors must point toward a shared coordination point (target His Nδ/Nε–metal ≈ 2.0–2.2 Å,
  His–metal–His ≈ 90–120°).
- Build a **small combinatorial library**: single His, His-pairs, and His-triplets per hotspot
  (typically 6–20 variants per region). Use site-directed / Kunkel / Golden-Gate or a spiked oligo.
- Keep WT and a distal-but-inert control mutant as negatives.
- Constructs: express the same tag/vector as WT for clean kinetic comparison. For **cpTEM-1** replication,
  build the circular permutant from the `1BTL` sequence and graft His₂/His₃ per the literature geometry.

## Stage 3 — metal screen on kinetics (biochemistry)
- Purify variants (IMAC/SEC); **strip adventitious metal** with EDTA then buffer-exchange into
  metal-free buffer (Chelex-treated, no Zn/Ni contamination) before adding defined metal.
- Assay ± **Zn²⁺ / Ni²⁺ / Co²⁺** across a concentration series (e.g. 0, 1, 10, 100, 1000 µM),
  measuring **kcat, Km, kcat/Km** by full Michaelis–Menten curves (not single-point).
  - GCK: coupled G6PDH (NADPH at 340 nm).
  - PTP1B: pNPP (405 nm) or DiFMUP (fluorescence).
  - AdK: coupled PK/LDH (NADH at 340 nm) or luciferase.
  - PFK: coupled aldolase/TIM/GDH (NADH at 340 nm).
  - cpTEM-1: nitrocefin or CENTA (chromogenic), purified enzyme only.
- Read-outs of interest:
  - Direction: **activation vs inhibition** (His-site geometry can flip the sign — key cpTEM-1 lesson).
  - Which parameter moves: **kcat** (catalytic step / conformational gating) vs **Km** (binding) vs **kcat/Km**.
  - Metal selectivity (Zn vs Ni vs Co) as evidence of a defined engineered site vs nonspecific binding.
- Confirm the mechanism is metal-mediated: reversibility on EDTA, dose-dependence, and no effect on the
  WT / inert-control mutant.

## Extending to oligomeric systems (PFK, ATCase, GP)
The current finder compares chain A only. For these, add:
- Multi-chain superposition and **inter-subunit contact-map differencing** (effector site ↔ active site).
- Interface His-sites spanning two chains (a metal bridging a subunit interface can tune the T↔R equilibrium).
- Symmetry-aware clustering so equivalent sites across protomers are merged.

## Caveats / good practice
- Cα displacement is a first filter; couple with contact-network / normal-mode analysis for subtle systems (PTP1B).
- Cα–Cα spacing is a coarse His-site proxy — always validate side-chain rotamer geometry before synthesis.
- Log every threshold; a silent top-N truncation hides real candidates. `results/<SYS>_regions.json`
  keeps the ranked list so nothing is dropped without a record.
