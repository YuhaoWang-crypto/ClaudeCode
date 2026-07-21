# GNoME-style materials discovery — reproduction summary

Learned from **Merchant et al., "Scaling deep learning for materials discovery"**,
*Nature* 623, 80–85 (2023), doi:10.1038/s41586-023-06735-9, and built a runnable skill
(`materials-discovery`) that reproduces the **discovery workflow**:

```
generate (structural substitution + compositional) → cheap prescreen (E_hull)
   → DFT-verify top-k → keep convex-hull points (E_hull ≈ 0) → retrain → repeat
```

## What's exact vs. what's a placeholder
- ✅ **Exact, unit-tested (8 tests):** convex-hull energy-above-hull LP, phase
  decomposition, ion-substitution generator, charge-neutrality + Goldschmidt
  tolerance-factor filters, compositional enumeration, formula algebra.
- ⚠️ **Offline surrogate energy** (Pauling-ionicity heuristic) — runs the whole loop
  without weights/cloud, but is **not DFT**. Real runs swap in UMA (prescreen) + Quantum
  ESPRESSO on Modal (DFT verification) — same scorer signature.

## Demo result (offline surrogate)
ABO₃ perovskite substitutions over {Ca,Sr,Ba}×{Ti,Zr}×O, ranked by E_hull vs the
binary oxides (CaO/SrO/BaO/TiO₂/ZrO₂):

| formula | E_form/atom (⚠️surrogate) | E_hull | verdict | decomposes to |
|---|---|---|---|---|
| BaZrO₃ | −3.583 | −0.211 | STABLE | ZrO₂ + BaO |
| BaTiO₃ | −3.267 | −0.159 | STABLE | TiO₂ + BaO |
| SrZrO₃ | −3.491 | −0.039 | STABLE | ZrO₂ + SrO |
| SrTiO₃ | −3.180 | +0.008 | STABLE | TiO₂ + SrO |
| CaZrO₃ | −3.416 | +0.140 | — | ZrO₂ + CaO |
| CaTiO₃ | −3.108 | +0.184 | — | TiO₂ + CaO |

The four flagged (BaZrO₃, BaTiO₃, SrZrO₃, SrTiO₃) are all real, well-known stable
perovskites — a plausible ranking even from the crude surrogate. ⚠️ The Ca perovskites
also exist experimentally but land above the surrogate hull — a concrete reminder that
resolving true ternary stability needs DFT/UMA, exactly why GNoME uses a trained GNN +
DFT active learning rather than a hand heuristic.

## Reproduce
```bash
python3 -m materials_discovery.test_materials   # exact tests
python3 -m materials_discovery.discover         # this demo
```

## Real discovery path
`reference/going-real.md`: UMA (fairchem, `uma-crystal-mof` skill) for the ML prescreen;
Quantum ESPRESSO on Modal (`qe-modal-bader-density` / `mlp-modal` skills) for the DFT
verification label; Materials-Project or same-settings DFT energies for the hull references.
