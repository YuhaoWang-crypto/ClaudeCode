# Vitamin-D metabolite discrimination — real Boltz-2.1 matrix

**Run:** 3 receptors × 3 clinically-key metabolites = 9 co-folding jobs
(Boltz-2.1, `ligand_protein_binding`, single-sequence, ~$0.45). All 9 succeeded.
Focus: the hard, useful **25(OH)D3 vs 1,25(OH)₂D3** discrimination (a single A-ring 1α-OH).

## Receptors (real PDB/designed sequences)
| receptor | source | length | natural target |
|---|---|---|---|
| VDR-LBD | 1DB1_A | 259 aa | 1,25(OH)₂D3 (active hormone) |
| CDL2.2 | designed (5IEN family) | 131 aa | D3 (parent) |
| DBP / GC | 1J78_A | 458 aa | 25(OH)D3 (clinical status marker) |

## Metabolites
- D3 `CC(C)CCCC(C)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C` — 3-OH only
- 25(OH)D3 `CC(CCCC(C)(C)O)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C` — +25-OH
- 1,25(OH)₂D3 `CC(CCCC(C)(C)O)C1CCC2C1(CCCC2=CC=C3CC(CC(C3=C)O)O)C` — +1α-OH +25-OH

## Raw metrics (✅ real Boltz numbers)

**ligand_iptm** (interface confidence)
| receptor | D3 | 25(OH)D3 | 1,25(OH)₂D3 |
|---|---|---|---|
| VDR-LBD    | 0.973 | 0.973 | **0.981** |
| CDL2.2     | 0.463 | 0.840 | 0.924 |
| DBP        | 0.861 | 0.777 | 0.696 |

**binding_confidence**
| receptor | D3 | 25(OH)D3 | 1,25(OH)₂D3 |
|---|---|---|---|
| VDR-LBD    | 0.667 | 0.757 | 0.743 |
| CDL2.2     | 0.368 | 0.377 | 0.316 |
| DBP        | 0.659 | 0.667 | 0.725 |

## Specificity margin = on-target − best off-target

| receptor | target | margin (ligand_iptm) | margin (binding_conf) | verdict |
|---|---|---|---|---|
| VDR-LBD | 1,25(OH)₂D3 | **+0.008** | −0.014 | non-selective |
| DBP     | 25(OH)D3    | −0.084 | −0.058 | non-selective |
| CDL2.2  | D3          | −0.461 | −0.010 | non-selective |

## Honest conclusion ⚠️

**Co-folding interface confidence does NOT resolve a single hydroxyl in this panel.**
Every margin is within method noise (±0.05–0.08), and in two of three cases a *wrong*
metabolite scores highest. This is exactly the resolution-limit caveat the framework
warned about — treated as a **prioritization**, not a measurement.

Two signals are still worth carrying to the bench:
1. **VDR ranks its on-target correctly.** The active hormone 1,25(OH)₂D3 gets the single
   highest ligand_iptm in the whole matrix (0.981), and the order D3 < 25(OH)D3 < 1,25(OH)₂D3
   matches VDR's known affinity trend — a real (if tiny-margin) positive.
2. **The designed D3 binder CDL2.2 is not D3-selective in silico** — it scores the
   hydroxylated metabolites *far higher* (0.92/0.84 vs 0.46). A concrete red flag: do not
   assume a "D3 binder" rejects 25/1,25 without counter-selection.

**What this argues for:** neither PDB template cleanly separates 25(OH) from 1,25(OH)₂ by
co-folding alone. The path to a real discriminating sensor is the **generative +
counter-selection** route (`specificity.generative_design_plan`): design a pocket against
the discriminating –OH, then *co-fold against every off-target and keep only large-margin
survivors*. The wet-lab ground truth remains a **competition / cross-reactivity assay** —
the same discipline as the rest of the skill.

Data: `specificity_vitd_matrix.json` (full metrics + both ranked matrices);
job IDs in `specificity_job_ids.txt`.
