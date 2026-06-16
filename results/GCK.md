# GCK — Human glucokinase (hexokinase IV)

**Conformer pair:** `1V4T:A` (super-open / inactive (apo))  ⇄  `1V4S:A` (closed / active (glucose + small-molecule activator))

- matched residues: **424**
- rigid-core residues: **244**  (core RMSD 0.58 Å, all-CA RMSD 11.07 Å)
- max Cα displacement after core fit: **39.1 Å**
- active-center = centroid of 10 active-site CA; distal cutoff 15 Å; responsive ≥ 10.83 Å (P75)

> _Expected signal:_ Large open->closed domain (large+small lobe) motion; activator pocket at the hinge/cleft.

## Distal candidate regions (ranked)

| # | residues | len | mean Δ (Å) | max Δ (Å) | min dist→active (Å) |
|---|----------|-----|------------|-----------|---------------------|
| 1 | 150–156 | 7 | 24.4 | 32.8 | 22.3 |
| 2 | 77–145 | 69 | 21.2 | 31.8 | 19.5 |
| 3 | 186–199 | 14 | 19.6 | 26.9 | 24.4 |
| 4 | 455–461 | 7 | 15.3 | 19.4 | 21.4 |
| 5 | 180–181 | 2 | 14.6 | 18.0 | 22.1 |
| 6 | 69–75 | 7 | 13.0 | 14.3 | 19.8 |

## Top hinge residues (clamp points)

Distal residues with the steepest local displacement gradient — where a moving domain meets the rigid core. Best single-metal clamp points.

| residue | resname | hinge score | Δ (Å) | dist→active (Å) |
|---------|---------|-------------|-------|-----------------|
| 180 | ASN | 4.53 | 5.2 | 22.1 |
| 151 | SER | 4.27 | 19.1 | 25.1 |
| 156 | HIS | 4.01 | 39.1 | 38.6 |
| 181 | VAL | 3.72 | 5.7 | 22.9 |
| 150 | PHE | 3.59 | 12.4 | 22.3 |
| 143 | LYS | 3.43 | 22.7 | 25.6 |
| 144 | LEU | 3.41 | 17.9 | 23.0 |
| 198 | ASP | 3.39 | 14.3 | 28.1 |
| 199 | VAL | 3.29 | 10.5 | 24.4 |
| 152 | PHE | 3.29 | 24.0 | 27.8 |
| 197 | MET | 3.09 | 21.2 | 28.9 |
| 145 | PRO | 3.07 | 12.6 | 19.5 |

## Suggested His-pair / His-triplet library

_Backbone-geometry heuristics — validate with rotamer-level metal modelling before synthesis. Screen each under Zn²⁺/Ni²⁺/Co²⁺ for kcat, Km, kcat/Km._

**His-pairs (bidentate clamp):**

| mutations | Cα–Cα (Å) | score |
|-----------|-----------|-------|
| 193H / 195H | 7.06 | 0.548 |
| 80H / 153H | 7.09 | 0.479 |
| 113H / 117H | 6.45 | 0.462 |
| 97H / 99H | 6.39 | 0.451 |
| 141H / 195H | 7.38 | 0.427 |
| 81H / 151H | 7.24 | 0.426 |
| 117H / 119H | 6.95 | 0.42 |
| 95H / 97H | 5.91 | 0.389 |
| 80H / 82H | 6.61 | 0.373 |
| 83H / 107H | 7.11 | 0.366 |
| 190H / 193H | 6.12 | 0.355 |
| 94H / 99H | 6.01 | 0.341 |
| 81H / 152H | 6.04 | 0.338 |
| 81H / 153H | 6.04 | 0.337 |
| 92H / 102H | 7.34 | 0.32 |

**His-triplets (facial site):**

| mutations | Cα–Cα pairwise (Å) | score |
|-----------|--------------------|-------|
| 83H / 108H / 110H | [6.19, 6.63, 6.06] | 0.587 |
| 81H / 151H / 153H | [7.24, 6.77, 6.04] | 0.571 |
| 190H / 194H / 196H | [7.32, 6.14, 6.39] | 0.544 |
| 92H / 94H / 99H | [6.13, 6.01, 5.8] | 0.451 |
| 100H / 102H / 459H | [6.44, 5.94, 7.65] | 0.411 |
| 82H / 109H / 111H | [6.32, 6.38, 8.05] | 0.378 |
| 82H / 111H / 115H | [8.05, 6.14, 6.47] | 0.356 |
| 113H / 115H / 117H | [6.08, 5.48, 6.45] | 0.354 |
| 111H / 113H / 115H | [5.5, 6.08, 6.14] | 0.351 |
| 94H / 98H / 100H | [5.49, 6.85, 6.03] | 0.334 |

## Ground-truth recovery (known allosteric site)

- known allosteric residues: 12
- recovered inside top regions: **3** (25%) → [455, 456, 459]
