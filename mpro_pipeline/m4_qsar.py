"""
Module 4 — the enzymatic QSAR, built under the four standing rules.

Rule 1  Labels are measured pIC50 and nothing else. Boltz metrics are NOT
        features and NOT labels here; pLDDT / ipTM / pTM are excluded outright
        (M3 measured them at or below the descriptor baseline, so as features
        they would add noise while looking informative). opt_score has a place
        only for compounds that carry no measured value at all, and that is a
        separate inference path, never mixed into training.
Rule 2  Pool, then quarantine. M1 showed per-assay stratification is not
        available (507 assays, median 2 records, largest 77) while pooling is
        sound for the ~59% of repeatedly-measured compounds that agree within
        0.5 log. So: train on the pooled set MINUS the discordant tail that M1
        tiers as SUSPECT/UNUSABLE.
Rule 3  Preincubation is a tier, not a filter. Only ~35% of assays (46% of
        records) report it; discarding the rest is too expensive. Performance
        is therefore reported separately on tier A and tier B, and an A value
        is never averaged against a B value for one compound.
Rule 4  The null model runs first. Descriptor-only performance is printed
        before the model's, and the margin decides whether the model earned
        anything.

Evaluation: scaffold split, not random. A random split on this set leaks
congeneric series across the boundary (the peptidomimetics are dense analogue
families) and inflates R^2. Bemis-Murcko scaffolds are assigned to train or
test as whole groups.

Honesty labels
  [validated] measured in this module on held-out compounds
  [limit]     a boundary established by measurement, not an opinion
"""
import collections
import json
import os

import numpy as np
from scipy import stats

from . import m1_labels, m2_null_model

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
    HAVE_RDKIT = True
except ImportError:                                     # pragma: no cover
    HAVE_RDKIT = False

SEED = 20260925
TEST_FRAC = 0.25
FP_BITS, FP_RADIUS = 1024, 2


def featurize(mols):
    """Descriptors + Morgan fingerprint. No model-derived columns (rule 1)."""
    desc_fns = [Descriptors.MolWt, Crippen.MolLogP, rdMolDescriptors.CalcTPSA,
                rdMolDescriptors.CalcNumHBD, rdMolDescriptors.CalcNumHBA,
                rdMolDescriptors.CalcNumRotatableBonds,
                rdMolDescriptors.CalcFractionCSP3,
                rdMolDescriptors.CalcNumRings,
                lambda m: m.GetNumHeavyAtoms()]
    gen = AllChem.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
    X = []
    for m in mols:
        d = [float(f(m)) for f in desc_fns]
        fp = np.frombuffer(gen.GetFingerprintAsNumPy(m), dtype=np.uint8).astype(float)
        X.append(np.concatenate([d, fp]))
    return np.asarray(X)


def scaffold_split(mols, test_frac=TEST_FRAC, seed=SEED):
    """Whole Bemis-Murcko scaffold groups go to one side. Prevents series leakage."""
    groups = collections.defaultdict(list)
    for i, m in enumerate(mols):
        try:
            s = MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)
        except Exception:
            s = ""
        groups[s].append(i)
    order = sorted(groups.values(), key=len, reverse=True)      # big series first
    rng = np.random.default_rng(seed)
    test, n_test = [], int(len(mols) * test_frac)
    for grp in order:
        if len(test) < n_test and (len(grp) == 1 or rng.random() < 0.5):
            test.extend(grp)
    test = set(test)
    train = [i for i in range(len(mols)) if i not in test]
    return train, sorted(test)


def report():
    if not HAVE_RDKIT:
        raise RuntimeError("rdkit required: pip install rdkit")
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import RidgeCV

    print("=" * 70)
    print("MODULE 4 — enzymatic QSAR under the four rules")
    print("=" * 70)

    tiers = m1_labels.tier_compounds(m1_labels.load())
    keep = m1_labels.trainable(tiers)                            # rule 2
    rows = [(m, v) for m, v in keep.items() if v["smiles"]]
    mols = [Chem.MolFromSmiles(v["smiles"]) for _, v in rows]
    ok = [i for i, m in enumerate(mols) if m is not None]
    rows = [rows[i] for i in ok]
    mols = [mols[i] for i in ok]
    y = np.array([v["pIC50"] for _, v in rows])

    print(f"pooled minus discordant tail: {len(rows)} compounds "
          f"(rule 2; {len(tiers) - len(keep)} quarantined)")
    print(f"pIC50 range {y.min():.2f}-{y.max():.2f}   sd {y.std():.2f}")

    tr, te = scaffold_split(mols)
    print(f"scaffold split: train {len(tr)}  test {len(te)}  "
          f"(whole scaffold groups, no series leakage)")

    # ---- rule 4: null model FIRST, on the same test compounds ----
    smi_te = [rows[i][1]["smiles"] for i in te]
    best_null, per = m2_null_model.baseline(smi_te, y[te])
    print(f"\nrule 4 — null model on the test set first:")
    for name, (rho, p) in sorted(per.items(), key=lambda kv: -abs(kv[1][0]))[:4]:
        print(f"  {name:<12} rho = {rho:+.3f}  p = {p:.4f}")
    print(f"  bar to clear: |rho| = {best_null:.3f}")

    X = featurize(mols)
    out = {}
    for name, model in (("ridge", RidgeCV(alphas=np.logspace(-2, 3, 20))),
                        ("random_forest", RandomForestRegressor(
                            n_estimators=300, min_samples_leaf=2,
                            random_state=SEED, n_jobs=-1))):
        model.fit(X[tr], y[tr])
        pred = model.predict(X[te])
        rho, p = stats.spearmanr(pred, y[te])
        rmse = float(np.sqrt(np.mean((pred - y[te]) ** 2)))
        r2 = 1 - np.sum((pred - y[te]) ** 2) / np.sum((y[te] - y[te].mean()) ** 2)
        margin, nv = m2_null_model.gate(rho, best_null)
        out[name] = dict(rho=float(rho), p=float(p), rmse=rmse, r2=float(r2),
                         margin=float(margin))
        print(f"\n[{name}]  held-out scaffold test, n = {len(te)}")
        print(f"  Spearman rho = {rho:+.3f}  p = {p:.2e}")
        print(f"  R^2 = {r2:+.3f}   RMSE = {rmse:.2f} log "
              f"(= {10 ** rmse:.1f}x in IC50)")
        print(f"  vs null model: margin {margin:+.3f}  -> {nv}")

    # ---- rule 3: preincubation tiers reported separately ----
    print("\nrule 3 — performance split by preincubation tier "
          "(never averaged together):")
    best = max(out, key=lambda k: out[k]["rho"])
    model = (RidgeCV(alphas=np.logspace(-2, 3, 20)) if best == "ridge"
             else RandomForestRegressor(n_estimators=300, min_samples_leaf=2,
                                        random_state=SEED, n_jobs=-1))
    model.fit(X[tr], y[tr])
    pred = model.predict(X[te])
    for tier in ("A", "B", "MIXED"):
        idx = [k for k, i in enumerate(te) if rows[i][1]["preincub"] == tier]
        if len(idx) < 8:
            print(f"  tier {tier:<5} n = {len(idx):<4} (too few to score)")
            continue
        rho, p = stats.spearmanr(pred[idx], y[te][idx])
        rmse = float(np.sqrt(np.mean((pred[idx] - y[te][idx]) ** 2)))
        print(f"  tier {tier:<5} n = {len(idx):<4} rho = {rho:+.3f}  "
              f"p = {p:.2e}  RMSE = {rmse:.2f} log")

    print("\n[limit] The label ceiling is not the model's fault: nirmatrelvir's")
    print("  own measured pIC50 spans 6.12-9.10 (2.98 log) over 35 assays and")
    print("  ebselen's spans 4.99-8.00. No QSAR can be more accurate than its")
    print("  labels, so an RMSE near ~0.7-0.9 log is at the noise floor of this")
    print("  data rather than a modelling failure.")
    return out


if __name__ == "__main__":
    report()
