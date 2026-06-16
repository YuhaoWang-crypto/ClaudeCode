# PTP1B — Protein tyrosine phosphatase 1B

**Conformer pair:** `2HNP:A` (apo (WPD loop open))  ⇄  `1T49:A` (allosteric-inhibitor bound (BB site))

- matched residues: **278**
- rigid-core residues: **247**  (core RMSD 0.21 Å, all-CA RMSD 0.57 Å)
- max Cα displacement after core fit: **5.0 Å**
- active-center = centroid of 11 active-site CA; distal cutoff 15 Å; responsive ≥ 0.28 Å (P75)

> _Expected signal:_ Smaller motion than GCK; signal lives in WPD-loop mobility + alpha7 ordering, not raw max displacement.

## Distal candidate regions (ranked)

| # | residues | len | mean Δ (Å) | max Δ (Å) | min dist→active (Å) |
|---|----------|-----|------------|-----------|---------------------|
| 1 | 277–282 | 6 | 1.4 | 2.0 | 19.8 |
| 2 | 235–245 | 11 | 1.3 | 2.5 | 22.4 |
| 3 | 61–66 | 6 | 0.6 | 0.8 | 24.4 |
| 4 | 147–151 | 5 | 0.5 | 0.5 | 18.4 |
| 5 | 122–123 | 2 | 0.3 | 0.4 | 16.4 |
| 6 | 129–131 | 3 | 0.3 | 0.3 | 29.0 |
| 7 | 12–16 | 5 | 0.3 | 0.3 | 17.5 |
| 8 | 47–48 | 2 | 0.3 | 0.3 | 15.4 |
| 9 | 205–207 | 3 | 0.3 | 0.3 | 25.1 |

## Top hinge residues (clamp points)

Distal residues with the steepest local displacement gradient — where a moving domain meets the rigid core. Best single-metal clamp points.

| residue | resname | hinge score | Δ (Å) | dist→active (Å) |
|---------|---------|-------------|-------|-----------------|
| 244 | VAL | 0.39 | 0.4 | 24.2 |
| 243 | SER | 0.38 | 0.6 | 27.8 |
| 237 | LYS | 0.37 | 0.3 | 29.0 |
| 238 | ARG | 0.36 | 0.6 | 30.2 |
| 236 | ASP | 0.35 | 0.5 | 27.0 |
| 242 | SER | 0.34 | 1.0 | 28.2 |
| 277 | GLY | 0.29 | 0.2 | 19.8 |
| 278 | ALA | 0.28 | 0.3 | 22.7 |
| 235 | MET | 0.28 | 0.2 | 25.2 |
| 245 | ASP | 0.24 | 0.3 | 22.4 |
| 279 | LYS | 0.23 | 1.4 | 24.8 |
| 280 | PHE | 0.23 | 2.6 | 25.5 |

## Suggested His-pair / His-triplet library

_Backbone-geometry heuristics — validate with rotamer-level metal modelling before synthesis. Screen each under Zn²⁺/Ni²⁺/Co²⁺ for kcat, Km, kcat/Km._

**His-pairs (bidentate clamp):**

| mutations | Cα–Cα (Å) | score |
|-----------|-----------|-------|
| 240H / 243H | 6.0 | 0.35 |
| 238H / 243H | 7.97 | 0.347 |

**His-triplets (facial site):**

| mutations | Cα–Cα pairwise (Å) | score |
|-----------|--------------------|-------|
| 238H / 240H / 243H | [5.04, 6.0, 7.97] | 0.123 |
| 238H / 240H / 242H | [5.04, 5.16, 9.27] | 0.017 |

## Ground-truth recovery (known allosteric site)

- known allosteric residues: 11
- recovered inside top regions: **2** (18%) → [280, 282]
