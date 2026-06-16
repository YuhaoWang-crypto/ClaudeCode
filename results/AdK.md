# AdK — E. coli adenylate kinase

**Conformer pair:** `4AKE:A` (apo / open)  ⇄  `1AKE:A` (Ap5A-bound / closed (transition-state mimic))

- matched residues: **214**
- rigid-core residues: **121**  (core RMSD 1.21 Å, all-CA RMSD 8.20 Å)
- max Cα displacement after core fit: **24.5 Å**
- active-center = centroid of 17 active-site CA; distal cutoff 15 Å; responsive ≥ 9.71 Å (P75)

> _Expected signal:_ Large LID (118-160) and NMP (30-67) closure over CORE; hinges at ~30/60/118/160; distal mobile 'counterweight'.

## Distal candidate regions (ranked)

| # | residues | len | mean Δ (Å) | max Δ (Å) | min dist→active (Å) |
|---|----------|-----|------------|-----------|---------------------|
| 1 | 123–156 | 34 | 16.2 | 23.2 | 17.0 |
| 2 | 39–58 | 20 | 12.6 | 15.6 | 18.6 |

## Top hinge residues (clamp points)

Distal residues with the steepest local displacement gradient — where a moving domain meets the rigid core. Best single-metal clamp points.

| residue | resname | hinge score | Δ (Å) | dist→active (Å) |
|---------|---------|-------------|-------|-----------------|
| 154 | THR | 2.25 | 14.9 | 22.9 |
| 153 | LEU | 2.15 | 17.1 | 25.1 |
| 155 | THR | 2.13 | 11.7 | 19.7 |
| 156 | ARG | 1.88 | 9.7 | 17.0 |
| 124 | ARG | 1.86 | 12.5 | 20.4 |
| 123 | ARG | 1.81 | 9.4 | 17.9 |
| 152 | GLU | 1.79 | 20.4 | 28.8 |
| 59 | VAL | 1.72 | 7.4 | 18.4 |
| 58 | LEU | 1.71 | 10.0 | 18.6 |
| 122 | GLY | 1.64 | 8.6 | 18.5 |
| 125 | VAL | 1.64 | 13.9 | 20.3 |
| 157 | LYS | 1.58 | 8.8 | 16.8 |

## Suggested His-pair / His-triplet library

_Backbone-geometry heuristics — validate with rotamer-level metal modelling before synthesis. Screen each under Zn²⁺/Ni²⁺/Co²⁺ for kcat, Km, kcat/Km._

**His-pairs (bidentate clamp):**

| mutations | Cα–Cα (Å) | score |
|-----------|-----------|-------|
| 56H / 58H | 6.66 | 0.51 |
| 42H / 44H | 6.6 | 0.44 |
| 128H / 130H | 5.98 | 0.397 |
| 150H / 152H | 6.7 | 0.394 |
| 52H / 56H | 6.27 | 0.387 |
| 127H / 154H | 6.78 | 0.35 |
| 148H / 150H | 6.01 | 0.331 |
| 142H / 150H | 8.03 | 0.322 |
| 141H / 147H | 6.16 | 0.318 |
| 51H / 56H | 8.43 | 0.235 |
| 147H / 150H | 5.63 | 0.206 |
| 42H / 47H | 5.56 | 0.198 |

**His-triplets (facial site):**

| mutations | Cα–Cα pairwise (Å) | score |
|-----------|--------------------|-------|
| 52H / 56H / 58H | [6.27, 6.66, 6.48] | 0.735 |
| 127H / 154H / 156H | [6.78, 6.59, 8.6] | 0.294 |
| 142H / 147H / 150H | [6.72, 5.63, 8.03] | 0.267 |
| 127H / 129H / 131H | [5.77, 5.54, 7.34] | 0.255 |
| 42H / 44H / 47H | [6.6, 5.43, 5.56] | 0.249 |
| 127H / 152H / 154H | [8.84, 6.82, 6.78] | 0.244 |
| 51H / 53H / 55H | [5.51, 5.41, 5.9] | 0.232 |
| 52H / 54H / 56H | [5.47, 5.28, 6.27] | 0.205 |
| 53H / 55H / 57H | [5.41, 5.31, 5.53] | 0.177 |
| 40H / 43H / 47H | [8.3, 6.06, 8.48] | 0.166 |
