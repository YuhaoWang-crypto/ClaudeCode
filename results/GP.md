# GP — Human liver glycogen phosphorylase

**Conformer pair:** `1FA9:A` (R-state (active-like))  ⇄  `3CEH:A` (T-state + allosteric inhibitor AVE5688 (AMP site))

- matched residues: **794**
- rigid-core residues: **727**  (core RMSD 0.99 Å, all-CA RMSD 1.76 Å)
- max Cα displacement after core fit: **14.9 Å**
- active-center = centroid of 5 active-site CA; distal cutoff 15 Å; responsive ≥ 1.35 Å (P75)

> _Expected signal:_ Inhibitor at AMP site stabilizes inactive tense conformation; recover the AMP/nucleotide allosteric pocket.

## Distal candidate regions (ranked)

| # | residues | len | mean Δ (Å) | max Δ (Å) | min dist→active (Å) |
|---|----------|-----|------------|-----------|---------------------|
| 1 | 279–290 | 12 | 7.8 | 13.3 | 18.3 |
| 2 | 249–252 | 4 | 4.2 | 5.7 | 41.2 |
| 3 | 261–270 | 10 | 3.4 | 6.1 | 35.3 |
| 4 | 379–386 | 8 | 3.1 | 4.2 | 15.6 |
| 5 | 324–326 | 3 | 2.5 | 2.9 | 42.7 |
| 6 | 314–316 | 3 | 2.4 | 2.9 | 42.4 |
| 7 | 205–214 | 10 | 2.1 | 2.7 | 44.0 |
| 8 | 192–197 | 6 | 2.1 | 2.4 | 34.2 |
| 9 | 713–742 | 30 | 2.0 | 2.8 | 17.7 |
| 10 | 757–770 | 14 | 1.9 | 2.3 | 17.2 |
| 11 | 829–831 | 3 | 1.9 | 2.0 | 31.2 |
| 12 | 590–597 | 8 | 1.9 | 2.3 | 28.9 |
| 13 | 552–557 | 6 | 1.8 | 2.0 | 23.0 |
| 14 | 406–436 | 31 | 1.8 | 2.5 | 31.2 |
| 15 | 772–780 | 9 | 1.6 | 1.8 | 16.0 |
| 16 | 62–66 | 5 | 1.6 | 1.8 | 33.2 |
| 17 | 41–45 | 5 | 1.6 | 1.7 | 44.2 |
| 18 | 131–132 | 2 | 1.5 | 1.5 | 15.2 |
| 19 | 745–751 | 7 | 1.5 | 1.5 | 24.4 |
| 20 | 165–166 | 2 | 1.4 | 1.4 | 17.5 |

## Top hinge residues (clamp points)

Distal residues with the steepest local displacement gradient — where a moving domain meets the rigid core. Best single-metal clamp points.

| residue | resname | hinge score | Δ (Å) | dist→active (Å) |
|---------|---------|-------------|-------|-----------------|
| 288 | GLY | 2.06 | 6.8 | 26.1 |
| 281 | PRO | 1.95 | 6.6 | 21.3 |
| 289 | LYS | 1.88 | 1.2 | 25.8 |
| 282 | ASN | 1.84 | 8.4 | 21.6 |
| 287 | GLU | 1.78 | 11.0 | 29.8 |
| 280 | TYR | 1.76 | 1.6 | 18.3 |
| 290 | GLU | 1.44 | 1.0 | 28.7 |
| 283 | ASP | 1.39 | 13.0 | 22.7 |
| 279 | LEU | 1.32 | 1.1 | 19.1 |
| 286 | PHE | 1.22 | 13.3 | 29.6 |
| 291 | LEU | 1.00 | 0.7 | 27.2 |
| 250 | ASN | 0.92 | 2.9 | 43.4 |

## Suggested His-pair / His-triplet library

_Backbone-geometry heuristics — validate with rotamer-level metal modelling before synthesis. Screen each under Zn²⁺/Ni²⁺/Co²⁺ for kcat, Km, kcat/Km._

**His-pairs (bidentate clamp):**

| mutations | Cα–Cα (Å) | score |
|-----------|-----------|-------|
| 593H / 596H | 7.16 | 0.544 |
| 433H / 435H | 7.22 | 0.54 |
| 723H / 725H | 6.24 | 0.48 |
| 209H / 214H | 6.97 | 0.464 |
| 284H / 288H | 6.64 | 0.445 |
| 722H / 725H | 6.23 | 0.433 |
| 41H / 45H | 6.49 | 0.425 |
| 261H / 265H | 5.95 | 0.386 |
| 212H / 214H | 6.38 | 0.368 |
| 208H / 212H | 6.01 | 0.353 |
| 721H / 725H | 6.06 | 0.353 |
| 418H / 424H | 7.82 | 0.291 |
| 210H / 212H | 5.63 | 0.286 |
| 207H / 212H | 8.35 | 0.238 |
| 426H / 433H | 8.36 | 0.228 |

**His-triplets (facial site):**

| mutations | Cα–Cα pairwise (Å) | score |
|-----------|--------------------|-------|
| 208H / 212H / 214H | [6.01, 6.38, 5.96] | 0.513 |
| 41H / 43H / 45H | [7.08, 5.68, 6.49] | 0.435 |
| 282H / 285H / 289H | [7.47, 5.83, 7.27] | 0.426 |
| 282H / 284H / 289H | [5.64, 6.76, 7.27] | 0.409 |
| 284H / 286H / 288H | [7.14, 5.56, 6.64] | 0.395 |
| 208H / 210H / 212H | [7.0, 5.63, 6.01] | 0.36 |
| 207H / 209H / 213H | [5.94, 5.78, 5.69] | 0.354 |
| 733H / 735H / 739H | [5.46, 6.29, 6.17] | 0.353 |
| 733H / 736H / 739H | [6.87, 5.37, 6.17] | 0.314 |
| 422H / 424H / 426H | [5.68, 5.76, 7.23] | 0.297 |

## Ground-truth recovery (known allosteric site)

- known allosteric residues: 11
- recovered inside top regions: **5** (45%) → [42, 43, 44, 45, 314]
