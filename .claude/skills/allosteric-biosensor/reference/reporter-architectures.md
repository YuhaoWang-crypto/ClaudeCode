# Biosensor reporter architectures & luminescent-reporter technology

Technical capabilities distilled from five papers, mapped onto the skill's
topology repertoire. These broaden the *readout* and *architecture* options
beyond single-component domain insertion.

## The architecture taxonomy (5 ways to couple binding → reporter output)

| # | architecture | mechanism | reporter examples | in this skill |
|---|---|---|---|---|
| 1 | **domain insertion** | binder (circularly permuted) inserted into a reporter loop; binding modulates activity | TEM-1 β-lac (colorimetric), PQQ-GDH (electrochemical), **cpFluc** (luminescent) | ✅ `insertion` mode |
| 2 | **CP-reporter + terminal fusion** | reporter circularly permuted; binder fused at a new terminus | NanoLuc | ✅ `cp_reporter_terminal` mode |
| 3 | **covalent / proteolytic activation** | cleavable linker in a dark cp-reporter; protease cleavage restores activity | cpFluc (CP234Luc) | ⬜ addable (`protease` config) |
| 4 | **split-fragment complementation (PCA)** | reporter split into 2 weak fragments; proximity/interaction reconstitutes activity | **NanoBiT** (LgBiT + SmBiT) | ⬜ addable (`split_complementation` topology) |
| 5 | **ligand-gated modules → mechanics/logic** | ligand-gated DNA-binding domains as switch "legs"; clocked multi-input control | TrpR, DtxR, MetJ | conceptual (logic gates / motors) |

## Promega luminescent-reporter lineage (Wood group)

### Firefly luciferase biosensors — 3 configurations (Fan et al., ACS Chem. Biol. 2008)
- **Circularly permuted firefly luciferase `CP234Luc`**: residues 234–544 placed
  *before* 4–233, joined by a polypeptide linker (new termini at the 234 loop).
  The cp-enzyme is conformationally constrained and nearly **dark** (~10³-fold
  below parental).
- **Covalent/protease config:** put a protease site in the CP linker → cleavage
  gives **185–2610-fold** luminescence increase. **Optimal linker 11–13 aa**
  (longer linker = higher basal, smaller activation — the conformational-
  constraint model; the same linker/DR trade-off the main paper and our
  `tune_linker` show).
- **Noncovalent config:** complementation (e.g. FRB/FKBP + rapamycin).
- **Allosteric config:** insert a ligand-binding domain (cAMP) into cpFluc →
  ligand modulates luminescence. This is architecture #1 with a luminescent
  reporter.

### Dynamic-range engineering (Binkowski et al., ACS Chem. Biol. 2011)
- Optimizing the cAMP cpFluc sensor (construct **22F** vs **20F**): signal/
  background **100× → 3500×** at saturation, EC50 **0.3 → 9 µM**, linear range
  **+30×**, in-cell response window **20-fold → 800-fold**.
- **Lesson for DR (ON/OFF):** the insertion geometry + linker set a
  sensitivity↔saturation trade-off. A "better" sensor is not the tightest binder
  but the one whose linear range matches the analyte's physiological span. This
  is precisely what our `tune_linker` / site×linker grid explores in silico —
  and why in-silico *binding retention* ≠ *dynamic range* (a high-affinity graft
  can saturate).

### NanoBiT — engineered split NanoLuc (Dixon et al., ACS Chem. Biol. 2016)
- NanoLuc dissected between **residues 156/157** → **LgBiT** (156 aa, 18 kDa;
  stability-optimized "11S") + **SmBiT** (13-aa native peptide, optimized to an
  **11-aa** peptide).
- Fragments engineered to associate **weakly by design**: intrinsic
  **KD ≈ 190 µM**, kon ≈ 500 M⁻¹s⁻¹, koff ≈ 0.2 s⁻¹ — *below* typical
  protein-interaction affinities, so complementation is dictated by the
  **appended proteins' interaction**, not the tags. Validated with the
  **SME-1 β-lactamase / BLIP** pair (β-lactamase recurs as a benchmark).
- Enables architecture #4: fuse LgBiT and SmBiT to two halves of a
  ligand-induced dimer (or to a binder + its target) → **ligand-induced
  proximity = luminescence**. This is the two-component counterpart to the main
  paper's single-component switches and its auxiliary-binding-domain designs.

## Ligand-gated modules & clocked control (Nilsson et al., Nat. Nanotechnol. 2026 — "Tumbleweed")
- An artificial protein **walker** whose three "legs" are **ligand-gated
  DNA-binding domains**: **TrpR** (gated by tryptophan), **DtxR** (Co²⁺/metal),
  **MetJ** (S-adenosylmethionine). Ligand present → foot binds its cognate DNA
  site; ligand withdrawn → foot releases.
- **Clocked control:** cycling ligand *pairs* (Trp+Co²⁺ → Co²⁺+SAM → SAM+Trp)
  drives directional 16-nm stepping; reversing the order reverses direction.
- **Capability for this skill:** a repertoire of **natural ligand-gated
  DNA-binding modules** (orthogonal small-molecule inputs) and the idea of
  **multi-input clocked/logic control** — directly relevant to the main paper's
  YES/AND logic gates and to building multi-analyte or sequential-logic sensors.

## What to add to the pipeline (menu)

- **cpFluc reporter** (`FLUC`, luminescent, insertion @ 234/544 CP site,
  linker 11–13 aa) — a drop-in `Reporter` like GDH/NanoLuc.
- **`protease` config** — cleavable-linker variant of a dark cp-reporter
  (architecture #3): a design function + a `protease_site` field.
- **`split_complementation` topology** (NanoBiT) — build two constructs
  (LgBiT-fusion + SmBiT-fusion) whose luminescence reports ligand-induced
  proximity (architecture #4).
- **ligand-gated module registry** — TrpR/DtxR/MetJ as orthogonal switch inputs
  for multi-analyte / logic designs.

Honesty labels carry through: sequence construction for any of these is ✅
deterministic; whether a given fusion *switches* (dynamic range) is ⚠️ until a
wet-lab luminescence titration — the Promega papers show DR is won by empirical
site/linker optimization, exactly what our in-silico scan can only *prioritize*.
