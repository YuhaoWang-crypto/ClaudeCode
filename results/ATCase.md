# ATCase — E. coli aspartate transcarbamoylase

**Conformer pair:** `6AT1:A` (T-state (apo/CTP))  ⇄  `1D09:A` (R-state + PALA (bisubstrate analogue))

- matched residues: **310**
- rigid-core residues: **257**  (core RMSD 0.80 Å, all-CA RMSD 2.45 Å)
- max Cα displacement after core fit: **13.3 Å**
- active-center = centroid of 9 active-site CA; distal cutoff 15 Å; responsive ≥ 1.22 Å (P75)

> _Expected signal:_ Global T<->R quaternary transition; regulatory-chain nucleotide site is the allosteric ground truth (separate chain).

## Distal candidate regions (ranked)

| # | residues | len | mean Δ (Å) | max Δ (Å) | min dist→active (Å) |
|---|----------|-----|------------|-----------|---------------------|
| 1 | 229–248 | 20 | 7.4 | 10.7 | 15.8 |
| 2 | 78–92 | 15 | 2.7 | 4.2 | 15.7 |
| 3 | 187–200 | 14 | 1.9 | 2.2 | 21.7 |
| 4 | 307–310 | 4 | 1.8 | 2.1 | 28.1 |
| 5 | 162–164 | 3 | 1.6 | 1.7 | 17.9 |
| 6 | 169–170 | 2 | 1.3 | 1.4 | 15.6 |
| 7 | 115–116 | 2 | 1.3 | 1.3 | 15.3 |

## Top hinge residues (clamp points)

Distal residues with the steepest local displacement gradient — where a moving domain meets the rigid core. Best single-metal clamp points.

| residue | resname | hinge score | Δ (Å) | dist→active (Å) |
|---------|---------|-------------|-------|-----------------|
| 246 | GLN | 1.64 | 4.2 | 28.2 |
| 247 | PHE | 1.52 | 1.4 | 26.1 |
| 245 | ALA | 1.43 | 10.6 | 32.0 |
| 248 | VAL | 1.25 | 0.5 | 25.0 |
| 233 | GLU | 1.18 | 7.4 | 19.0 |
| 232 | LYS | 1.16 | 6.0 | 19.6 |
| 231 | GLN | 1.13 | 4.8 | 16.5 |
| 234 | ARG | 1.03 | 7.4 | 18.9 |
| 230 | VAL | 1.02 | 2.8 | 17.9 |
| 244 | LYS | 1.01 | 12.4 | 31.8 |
| 229 | ARG | 0.90 | 1.4 | 15.8 |
| 249 | LEU | 0.89 | 0.3 | 26.0 |

## Suggested His-pair / His-triplet library

_Backbone-geometry heuristics — validate with rotamer-level metal modelling before synthesis. Screen each under Zn²⁺/Ni²⁺/Co²⁺ for kcat, Km, kcat/Km._

**His-pairs (bidentate clamp):**

| mutations | Cα–Cα (Å) | score |
|-----------|-----------|-------|
| 81H / 85H | 6.98 | 0.531 |
| 83H / 85H | 6.9 | 0.529 |
| 78H / 82H | 6.9 | 0.433 |
| 235H / 239H | 6.44 | 0.432 |
| 82H / 85H | 6.46 | 0.42 |
| 82H / 87H | 7.03 | 0.42 |
| 81H / 86H | 6.37 | 0.418 |
| 79H / 82H | 6.66 | 0.417 |
| 78H / 83H | 6.06 | 0.366 |
| 308H / 310H | 6.19 | 0.358 |
| 231H / 234H | 6.03 | 0.357 |
| 85H / 87H | 5.82 | 0.307 |
| 80H / 82H | 5.65 | 0.247 |
| 85H / 90H | 8.74 | 0.163 |

**His-triplets (facial site):**

| mutations | Cα–Cα pairwise (Å) | score |
|-----------|--------------------|-------|
| 82H / 85H / 87H | [6.46, 5.82, 7.03] | 0.496 |
| 232H / 235H / 237H | [5.83, 6.24, 6.91] | 0.483 |
| 81H / 83H / 86H | [5.64, 6.75, 6.37] | 0.44 |
| 81H / 83H / 85H | [5.64, 6.9, 6.98] | 0.438 |
| 232H / 234H / 236H | [5.92, 7.02, 7.7] | 0.425 |
| 81H / 84H / 86H | [6.05, 5.62, 6.37] | 0.398 |
| 78H / 80H / 83H | [5.72, 7.27, 6.06] | 0.361 |
| 239H / 241H / 243H | [5.43, 6.02, 6.74] | 0.322 |
| 78H / 80H / 82H | [5.72, 5.65, 6.9] | 0.315 |
| 235H / 237H / 239H | [6.24, 5.31, 6.44] | 0.313 |
