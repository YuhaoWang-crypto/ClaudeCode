# GNoME methodology → this skill

Source: Merchant, Batzner, Schoenholz, Aykol, Cheon, Cubuk, **"Scaling deep learning
for materials discovery"**, *Nature* 623, 80–85 (2023). doi:10.1038/s41586-023-06735-9.

## What GNoME does
A deep-learning **active-learning engine** for inorganic crystal discovery:

1. **Two candidate pipelines.**
   - *Structural*: substitute ions into known crystal prototypes (dominant source of hits).
   - *Compositional*: enumerate/─randomize chemical formulas (template-free arm).
2. **GNN stability predictor.** Graph networks predict formation energy / **energy above
   the convex hull** (E_hull) for millions of candidates cheaply.
3. **DFT verification + active learning.** The most promising candidates are labeled with
   DFT; labels retrain the GNN. Discovery hit-rate rose from <10 % to **>80 %**; stability
   MAE fell to ~11 meV/atom.
4. **Convex-hull filter.** Keep candidates with E_hull ≈ 0 — they don't decompose into any
   mixture of known phases. Result: **2.2 M crystals, ~380 k stable**, incl. 528 Li-ion
   conductors and 52 k graphene-like layered materials.

## How each piece is realized here
| GNoME piece | This skill |
|---|---|
| Structural pipeline | `generate.structural()` — ion substitution into `prototypes.PROTOTYPES` with charge-neutrality + Goldschmidt tolerance-factor filters |
| Compositional pipeline | `generate.compositional()` — charge-neutral integer stoichiometries over an element palette |
| GNN stability predictor | `score.py` scorer — **surrogate** (offline demo) or **UMA** (`uma_energy`, real ML) |
| DFT verification label | `score.qe_modal_energy` — Quantum ESPRESSO SCF on Modal |
| E_hull / convex hull | `hull.PhaseDiagram` — exact LP; `e_above_hull`, `decomposition` |
| Active-learning loop | `discover.discover()` / `screen()` |

## The convex-hull decision (exact)
Energies are **formation energy per atom** (eV/atom), elements = zero reference (MP
convention). The hull energy at composition *x* is the minimum energy reachable by any
phase mixture matching *x*, found by linear programming over atom-fraction weights
`a_i ≥ 0, Σa_i = 1, Σa_i f_i = x`. Then `E_hull(x) = ε(x) − hull(x)`. A candidate with
`E_hull ≤ tol` (excluding itself) **defines a new hull point** — a predicted-stable
discovery. This math is unit-tested against synthetic systems with known answers.

## Honesty labels
- ✅ **Exact / tested**: hull LP, decomposition, substitution enumeration, charge & tolerance filters, formula algebra.
- ⚠️ **Surrogate**: the offline ionicity energy (Pauling Δχ heuristic) is monotonic and reproducible but **not DFT** — it lets the loop run without weights/cloud and must not be read as a real stability prediction. In the demo, perovskites land near E_hull ≈ 0, which is precisely why real GNoME needs a trained GNN + DFT.
- ✅ **Real path**: UMA for the ML prescreen, DFT (QE/Modal) for the verification label; the deciding number for any stability claim is a DFT E_hull, and ultimately experimental synthesis.
