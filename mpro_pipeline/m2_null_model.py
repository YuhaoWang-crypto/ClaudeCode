"""
Module 2 — the null-model gate.

Rule enforced: NOTHING in this pipeline may report a rho or an R^2 before the
plain-descriptor baseline on the SAME compounds has been computed and printed
next to it. A structure-based or ML score only earns its keep by the margin it
puts on top of descriptors that need no structure, no model and no GPU.

Why this is not optional here. On the 11-compound calibration set the
strongest single descriptor is cLogP at |rho| = 0.410 (p = 0.21), and Boltz's
own predicted lipophilicity reaches |rho| = 0.656 (p = 0.028). Ranking Mpro
inhibitors by "less greasy is better" therefore scores ~0.41-0.66 while
knowing nothing about the protein. The reason is a confound in the compound
set rather than chemistry: the potent compounds are polar peptidomimetics
(nirmatrelvir, PF-00835231, GC373) and the weak ones are greasier repurposing
candidates (masitinib, telaprevir, disulfiram). Any headline "our score
correlates at rho = 0.6" is indistinguishable from that confound until this
baseline is printed.

gate() returns the margin. A negative or ~0 margin means the expensive score
carries no structural information over descriptors, whatever its absolute rho.
"""
import json
import os

import numpy as np
from scipy import stats

try:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
    HAVE_RDKIT = True
except ImportError:                                    # pragma: no cover
    HAVE_RDKIT = False

MARGIN_REAL = 0.05          # margin below this = indistinguishable from null

# Same label aliasing as M3, so the bar is computed on the SAME compounds the
# model scores. Without this, M2 silently dropped GC373 (its measured label is
# filed under the GC376 prodrug) and reported a bar from n=10 while M3 scored
# n=11 -- two different bars for one rule, which defeats the point of the gate.
ALIAS = {"GC373_aldehyde": "GC376"}
CONTROLS = {"nirmatrelvir_epimer_CTRL"}     # method control, no measured label

DESCRIPTORS = {
    "MW": lambda m: Descriptors.MolWt(m),
    "cLogP": lambda m: Crippen.MolLogP(m),
    "TPSA": lambda m: rdMolDescriptors.CalcTPSA(m),
    "HBD": lambda m: rdMolDescriptors.CalcNumHBD(m),
    "HBA": lambda m: rdMolDescriptors.CalcNumHBA(m),
    "heavy_atoms": lambda m: m.GetNumHeavyAtoms(),
    "rot_bonds": lambda m: rdMolDescriptors.CalcNumRotatableBonds(m),
    "fsp3": lambda m: rdMolDescriptors.CalcFractionCSP3(m),
}


def baseline(smiles, y):
    """Best |rho| over plain descriptors. Returns (best_abs_rho, per_descriptor)."""
    if not HAVE_RDKIT:
        raise RuntimeError("rdkit required for the null-model gate: pip install rdkit")
    mols = [Chem.MolFromSmiles(s) if s else None for s in smiles]
    ok = [i for i, m in enumerate(mols) if m is not None]
    yv = np.asarray(y, dtype=float)[ok]

    per = {}
    for name, fn in DESCRIPTORS.items():
        try:
            x = np.array([fn(mols[i]) for i in ok], dtype=float)
        except Exception:
            continue
        if np.all(x == x[0]):
            continue
        rho, p = stats.spearmanr(x, yv)
        if not np.isnan(rho):
            per[name] = (float(rho), float(p))
    best = max((abs(r) for r, _ in per.values()), default=0.0)
    return best, per


def gate(score_rho, best_null_rho):
    """Margin of a model score over the null model, plus a verdict."""
    margin = abs(score_rho) - best_null_rho
    if margin > MARGIN_REAL:
        verdict = "beats null model"
    elif margin > -MARGIN_REAL:
        verdict = "INDISTINGUISHABLE from null model"
    else:
        verdict = "WORSE than null model"
    return margin, verdict


def report():
    here = os.path.dirname(__file__)
    with open(os.path.join(here, "data", "ground_truth.json")) as fh:
        gt = json.load(fh)
    with open(os.path.join(here, "data", "predictions.json")) as fh:
        pred = json.load(fh)

    # the calibration set's SMILES are the ones actually submitted to Boltz
    with open(os.path.join(here, "data", "calib_smiles.json")) as fh:
        smi = json.load(fh)

    names = [n for n in smi
             if n not in CONTROLS and ALIAS.get(n, n) in gt]
    y = [gt[ALIAS.get(n, n)]["pIC50"] for n in names]

    print("=" * 70)
    print("MODULE 2 — null-model gate (no rho is reportable without this)")
    print("=" * 70)
    print(f"compounds with measured pIC50: n = {len(names)}")

    best, per = baseline([smi[n] for n in names], y)
    print("\nplain descriptors, no structure used:")
    for name, (rho, p) in sorted(per.items(), key=lambda kv: -abs(kv[1][0])):
        print(f"  {name:<12} rho = {rho:+.3f}   |rho| = {abs(rho):.3f}   p = {p:.4f}")
    print(f"\n  strongest null model |rho| = {best:.3f}"
          f"   <- the bar every model score must clear")

    print("\nBoltz metrics measured against that bar:")
    for cfg, per_mol in sorted(pred.items()):
        for metric in ("optimization_score", "binding_confidence", "iptm",
                       "complex_plddt"):
            rows = [(per_mol[n][metric], gt[ALIAS.get(n, n)]["pIC50"])
                    for n in names
                    if n in per_mol and metric in per_mol[n]]
            if len(rows) < 5:
                continue
            rho, p = stats.spearmanr([r[0] for r in rows], [r[1] for r in rows])
            margin, verdict = gate(rho, best)
            print(f"  {cfg:<12} {metric:<20} rho = {rho:+.3f}  p = {p:.4f}  "
                  f"margin {margin:+.3f}  {verdict}")

    return dict(best_null=best, per_descriptor=per)


if __name__ == "__main__":
    report()
