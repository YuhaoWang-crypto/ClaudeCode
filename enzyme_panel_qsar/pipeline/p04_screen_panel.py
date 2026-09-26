"""Stage 4 - score the clinical library against every modellable target.

This is the high-throughput part: one library, every target, one pass. The
library is featurised once and reused, because every target model was built on
the identical fingerprint definition.

Two things are reported per prediction, and they are deliberately not merged:

**Reliability** (`ad_tier`) - how far the compound sits from the training set,
by maximum Tanimoto to any training compound. Closer means the model is
interpolating and the prediction is more trustworthy.

**Novelty** (`novel_scaffold`) - whether the compound's Murcko scaffold is absent
from the training scaffolds. Interesting for discovery, and orthogonal to
reliability.

The 3CLpro work collapsed these into one tier, which made near-analogues of
training compounds come out as `out_of_domain`; its own limitations section then
had to explain that those predictions are in fact the *most* reliable ones and
that the tier should be read backwards for repurposing. Keeping the two axes
separate removes the need for that caveat: a high-similarity, non-novel compound
is exactly what a repurposing screen wants, and it reads that way here.

`pred_sd` is the random-forest per-tree spread. It is an uncertainty estimate,
not a second prediction - the point value always comes from the target's
selected model.
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p02_benchmark import DESCRIPTORS, FP_BITS, FP_RADIUS, featurise  # noqa: E402

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
MODELS = ROOT / "models"

# Reliability thresholds on max-Tanimoto to the training set. 0.4 on ECFP4 is a
# conventional "same chemical series" boundary; below ~0.25 a fingerprint model
# is extrapolating and the number should not be used to make a decision.
AD_HIGH = 0.40
AD_BORDERLINE = 0.25

_gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
_pains = FilterCatalog(
    FilterCatalogParams(FilterCatalogParams.FilterCatalogs.PAINS)
)


def bitvects(smiles: list[str]):
    out = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        out.append(_gen.GetFingerprint(mol) if mol is not None else None)
    return out


def pains_flags(smiles: list[str]) -> list[str]:
    flags = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            flags.append("")
            continue
        hits = _pains.GetMatches(mol)
        flags.append(";".join(sorted({h.GetDescription() for h in hits}))[:120])
    return flags


def main() -> None:
    lib = pd.read_csv(DATA / "library_clinical.csv")
    bench = {b["name"]: b for b in json.loads((RESULTS / "benchmark.json").read_text())}
    panel = json.loads((ROOT / "configs" / "panel.json").read_text())
    controls = {
        t["name"]: [(c["name"].upper(), c["expect"]) for c in t["controls"]]
        for t in panel["targets"]
    }

    prov = json.loads((RESULTS / "panel_provenance.json").read_text())
    noise_floors = {
        t["name"]: t["noise_floor"]["noise_floor_log"] for t in prov["targets"]
    }
    modellable = [n for n, b in bench.items() if b["modellable"]]
    skipped = [n for n, b in bench.items() if not b["modellable"]]
    print(f"screening {len(lib)} library compounds against {len(modellable)} targets")
    if skipped:
        print(
            f"  excluded (did not clear the null, or rho < 0.3): {', '.join(skipped)}"
        )

    print("featurising library ...")
    fp, desc = featurise(lib["std_smiles"].tolist())
    Xlib = np.hstack([fp, desc]).astype(np.float32)
    lib_bv = bitvects(lib["std_smiles"].tolist())
    lib["pains_alert"] = pains_flags(lib["std_smiles"].tolist())

    long_rows = []
    control_recall = {}
    for name in modellable:
        with open(MODELS / f"{name}.pkl", "rb") as fh:
            art = pickle.load(fh)
        pred = art["model"].predict(Xlib)
        # Per-tree spread from the forest kept for exactly this purpose.
        trees = np.stack([t.predict(Xlib) for t in art["uncertainty_rf"].estimators_])
        sd = trees.std(axis=0)

        train_bv = bitvects(art["train_smiles"])
        train_bv = [b for b in train_bv if b is not None]
        nn = np.array(
            [
                max(DataStructs.BulkTanimotoSimilarity(b, train_bv)) if b is not None else 0.0
                for b in lib_bv
            ]
        )
        train_scafs = set(art["train_scaffolds"])
        train_smi = set(art["train_smiles"])
        # Measured labels, for scoring controls against their own value rather
        # than against a percentile band.
        tr = pd.read_csv(DATA / f"train_{name}.csv")
        measured = dict(zip(tr["molecule_chembl_id"], tr["pAffinity"]))
        noise = noise_floors.get(name) or 0.5

        tier = np.where(
            nn >= AD_HIGH, "high", np.where(nn >= AD_BORDERLINE, "borderline", "out_of_domain")
        )
        sub = pd.DataFrame(
            {
                "target": name,
                "chembl_id": lib["chembl_id"],
                "pref_name": lib["pref_name"],
                "max_phase": lib["max_phase"],
                "pred_pAffinity": np.round(pred, 3),
                "pred_sd": np.round(sd, 3),
                "nn_tanimoto": np.round(nn, 3),
                "ad_tier": tier,
                "novel_scaffold": ~lib["scaffold"].isin(train_scafs),
                "in_training": lib["std_smiles"].isin(train_smi),
                "pains_alert": lib["pains_alert"],
                "std_smiles": lib["std_smiles"],
                "mw": lib["mw"],
            }
        )
        long_rows.append(sub)

        # Did the virtual assay rank this target's real reference inhibitors as
        # potent? Controls were never held out, so an in-training control is a
        # sanity check on the fit, not evidence of generalisation - the flag says
        # which is which.
        names = lib["pref_name"].fillna("").str.upper()
        rec = []
        for ctrl, expect in controls.get(name, []):
            hit = sub[names.eq(ctrl)]
            if hit.empty:
                rec.append({"control": ctrl, "expect": expect, "in_library": False})
                continue
            r = hit.iloc[0]
            pct = float((sub["pred_pAffinity"] < r["pred_pAffinity"]).mean() * 100)
            rec.append(
                {
                    "control": ctrl,
                    "expect": expect,
                    "in_library": True,
                    "pred_pAffinity": float(r["pred_pAffinity"]),
                    "percentile_in_library": round(pct, 1),
                    "ad_tier": r["ad_tier"],
                    "in_training": bool(r["in_training"]),
                    "measured_pAffinity": (
                        round(float(measured[r["chembl_id"]]), 3)
                        if r["chembl_id"] in measured
                        else None
                    ),
                }
            )
            # Scoring: where the control has a measured label, the honest test is
            # whether the prediction reproduces it to within the assay noise
            # floor. A percentile band was tried first and produced boundary
            # artifacts - cilomilast landed at exactly the 75th percentile, and
            # rolipram's correct ~1 uM prediction failed a "weak" band purely
            # because a clinical library skews weaker than that. Percentile is
            # kept as context, not as the criterion.
            rec[-1]["residual_vs_measured"] = (
                round(float(r["pred_pAffinity"]) - rec[-1]["measured_pAffinity"], 3)
                if rec[-1]["measured_pAffinity"] is not None
                else None
            )
            rec[-1]["as_expected"] = bool(
                abs(rec[-1]["residual_vs_measured"]) <= 1.5 * noise
                if rec[-1]["residual_vs_measured"] is not None
                else (pct >= 75 if expect == "potent" else pct <= 50)
            )
            _ = (
            )
        control_recall[name] = rec
        found = [r for r in rec if r.get("in_library")]
        ok = [r for r in found if r["as_expected"]]
        res = [r["residual_vs_measured"] for r in found if r.get("residual_vs_measured") is not None]
        print(
            f"  {name:<10s} pred {pred.min():.2f}-{pred.max():.2f}  "
            f"AD high {int((tier == 'high').sum()):>5d}  "
            f"controls {len(found)}/{len(controls.get(name, []))} found, "
            f"{len(ok)}/{len(found)} reproduced"
            + (f"  median |resid| {np.median(np.abs(res)):.2f} log" if res else "")
        )

    long = pd.concat(long_rows, ignore_index=True)
    long.to_csv(RESULTS / "panel_predictions_long.csv.gz", index=False)

    wide = long.pivot_table(index="chembl_id", columns="target", values="pred_pAffinity")
    tiers = long.pivot_table(
        index="chembl_id", columns="target", values="ad_tier", aggfunc="first"
    )
    meta = (
        long.groupby("chembl_id")
        .agg(pref_name=("pref_name", "first"), max_phase=("max_phase", "first"),
             pains_alert=("pains_alert", "first"), std_smiles=("std_smiles", "first"))
    )
    wide = meta.join(wide)
    wide.to_csv(RESULTS / "panel_matrix_pAffinity.csv")
    tiers.to_csv(RESULTS / "panel_matrix_ad_tier.csv")

    (RESULTS / "control_recall.json").write_text(json.dumps(control_recall, indent=2))
    print(
        f"\nwrote panel_predictions_long.csv.gz ({long.shape}), "
        f"panel_matrix_pAffinity.csv ({wide.shape}), control_recall.json"
    )


if __name__ == "__main__":
    main()
