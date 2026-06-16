# PFK — E. coli phosphofructokinase-1

**Conformer pair:** `2PFK:A` (unliganded)  ⇄  `1PFK:A` (products + allosteric activator ADP/Mg2+)

- matched residues: **301**
- rigid-core residues: **240**  (core RMSD 0.38 Å, all-CA RMSD 0.67 Å)
- max Cα displacement after core fit: **2.1 Å**
- active-center = centroid of 8 active-site CA; distal cutoff 15 Å; responsive ≥ 0.63 Å (P75)

> _Expected signal:_ Effector-site and subunit-interface contact-map changes; analyze per-chain then interfaces, not raw global RMSD.

## Distal candidate regions (ranked)

| # | residues | len | mean Δ (Å) | max Δ (Å) | min dist→active (Å) |
|---|----------|-----|------------|-----------|---------------------|
| 1 | 106–118 | 13 | 1.3 | 1.7 | 15.6 |
| 2 | 195–216 | 22 | 1.1 | 1.5 | 20.1 |
| 3 | 296–300 | 5 | 1.0 | 1.3 | 16.6 |
| 4 | 73–90 | 18 | 1.0 | 1.3 | 17.6 |
| 5 | 287–290 | 4 | 0.7 | 0.8 | 21.1 |
| 6 | 1–2 | 2 | 0.7 | 0.7 | 34.0 |
| 7 | 58–60 | 3 | 0.7 | 0.7 | 23.3 |
| 8 | 185–186 | 2 | 0.6 | 0.6 | 20.9 |

## Top hinge residues (clamp points)

Distal residues with the steepest local displacement gradient — where a moving domain meets the rigid core. Best single-metal clamp points.

| residue | resname | hinge score | Δ (Å) | dist→active (Å) |
|---------|---------|-------------|-------|-----------------|
| 118 | PRO | 0.25 | 0.9 | 27.1 |
| 119 | CYS | 0.24 | 0.3 | 23.4 |
| 117 | PHE | 0.23 | 1.3 | 28.1 |
| 215 | HIS | 0.22 | 1.0 | 28.5 |
| 216 | ALA | 0.21 | 0.4 | 25.2 |
| 116 | GLY | 0.19 | 1.7 | 30.5 |
| 217 | ILE | 0.17 | 0.2 | 21.6 |
| 120 | ILE | 0.17 | 0.2 | 21.5 |
| 214 | LYS | 0.15 | 1.9 | 32.2 |
| 296 | ILE | 0.15 | 0.3 | 20.1 |
| 295 | ASP | 0.15 | 0.4 | 23.0 |
| 195 | GLU | 0.15 | 0.8 | 20.5 |

## Suggested His-pair / His-triplet library

_Backbone-geometry heuristics — validate with rotamer-level metal modelling before synthesis. Screen each under Zn²⁺/Ni²⁺/Co²⁺ for kcat, Km, kcat/Km._

**His-pairs (bidentate clamp):**

| mutations | Cα–Cα (Å) | score |
|-----------|-----------|-------|
| 208H / 212H | 6.54 | 0.417 |
| 75H / 81H | 6.88 | 0.378 |
| 204H / 208H | 6.37 | 0.343 |
| 212H / 214H | 5.96 | 0.299 |
| 208H / 210H | 5.55 | 0.217 |
| 207H / 212H | 8.73 | 0.175 |
| 208H / 213H | 5.51 | 0.17 |
| 208H / 214H | 8.61 | 0.145 |
| 203H / 208H | 8.74 | 0.128 |

**His-triplets (facial site):**

| mutations | Cα–Cα pairwise (Å) | score |
|-----------|--------------------|-------|
| 196H / 198H / 200H | [6.68, 5.46, 7.78] | 0.269 |
| 199H / 201H / 203H | [5.57, 5.41, 6.41] | 0.252 |
| 77H / 79H / 81H | [6.71, 5.44, 7.89] | 0.245 |
| 75H / 78H / 81H | [5.71, 5.35, 6.88] | 0.244 |
| 198H / 200H / 202H | [5.46, 5.45, 6.15] | 0.237 |
| 197H / 199H / 201H | [5.33, 5.57, 6.25] | 0.234 |
| 208H / 210H / 213H | [5.55, 7.22, 5.51] | 0.224 |
| 208H / 211H / 213H | [5.33, 6.23, 5.51] | 0.222 |
| 200H / 202H / 204H | [5.45, 5.36, 6.29] | 0.219 |
| 207H / 209H / 211H | [5.39, 5.41, 6.43] | 0.217 |

## Ground-truth recovery (known allosteric site)

- known allosteric residues: 10
- recovered inside top regions: **4** (40%) → [59, 211, 213, 214]
