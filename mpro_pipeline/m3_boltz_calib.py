"""
Module 3 — calibrate Boltz-2 metrics against measured enzymatic pIC50.

Validate the tool before trusting its ranking. Same discipline as
grn_pipeline/m10_validate, which found on KRAS G12C that opt_score tracked
measured potency (rho = +0.6) while binding_confidence anti-correlated
(rho = -0.2): metrics from ONE run disagree with each other, so "Boltz says
so" is not a result until you know which number to read.

PRE-REGISTERED criteria, fixed before any result was seen:
    primary     Spearman rho(metric, measured pIC50)
    PASS        rho >= +0.5 and p < 0.05      -> usable to rank and pick
    GREY        +0.3 <= rho < +0.5            -> coarse triage only
    FAIL        rho < +0.3                    -> metric unusable
Every reported rho is additionally passed through the m2 null-model gate; a
rho that does not clear the descriptor baseline is not evidence of anything.

TARGET SETUP
  3CLpro is a homodimer and only the dimer is catalytically competent: the
  N-finger (Ser1-Phe3) of one protomer completes the S1 pocket of the other.
  Docking into a monomer leaves S1 unformed, so the dimer (chains A+B, one
  sequence) is the correct construct and the monomer is kept as a control.
  Pocket residues (0-indexed, chain A): His41 Met49 Phe140 Leu141 Asn142
  Gly143 Ser144 Cys145 His163 His164 Met165 Glu166 Leu167 Pro168 His172
  Asp187 Arg188 Gln189 Thr190 Gln192.

CONFIGURATIONS PRESENT IN data/predictions.json
  dimer_AB_v2  authoritative: dimer, correct nirmatrelvir stereochemistry
  dimer_AB     first pass; carried the wrong gamma-lactam epimer (see below)
  monomer_A    single-chain control, same pocket, same ligands

TWO CORRECTIONS RECORDED HERE
  1. ChEMBL molecule/search returns an unnamed duplicate of nirmatrelvir
     (CHEMBL5201264, pref_name None, 2 records) ahead of the curated entry
     (CHEMBL4802135, NIRMATRELVIR, 35 records). They differ at one
     gamma-lactam stereocentre. Taking result [0] silently submitted the wrong
     epimer. Never trust search ordering: check pref_name and record count.
  2. Measured against the epimer control in dimer_AB_v2, Boltz barely
     separates them: opt_score 0.624 vs 0.616 (delta +0.008), and
     binding_confidence is 0.999 for both. So the error cost almost nothing
     numerically -- but it also means BOLTZ-2 DOES NOT RANK STEREOISOMERS
     here. Do not ask it to do stereo-SAR.

  ! SMARTS filtering must be disabled for a calibration set: at the default
    'recommended' level ebselen (Se) and disulfiram (thiuram disulfide) are
    dropped before prediction and the weak end of the potency range vanishes.
"""
import json
import os

import numpy as np
from scipy import stats

from . import m2_null_model

PASS_RHO, GREY_RHO, PASS_P = 0.5, 0.3, 0.05
BOOT = 20000
SEED = 20260925

# GC376 is the bisulfite prodrug of the GC373 aldehyde and dissociates to it in
# assay buffer, so the measured GC376 IC50 is the aldehyde's label.
ALIAS = {"GC373_aldehyde": "GC376"}
# thiol-reactive promiscuous compounds: labels unreliable at any model quality
ARTIFACTS = {"ebselen", "disulfiram", "carmofur"}
# the epimer is a method control, not a compound with its own measured value
CONTROLS = {"nirmatrelvir_epimer_CTRL"}


def verdict(rho, p):
    if rho >= PASS_RHO and p < PASS_P:
        return "PASS  (usable to rank)"
    if rho >= GREY_RHO:
        return "GREY  (coarse triage only)"
    return "FAIL  (unusable)"


def _paired(per_mol, gt, metric):
    names, x, y = [], [], []
    for mol, sc in per_mol.items():
        if mol in CONTROLS or metric not in sc or sc[metric] is None:
            continue
        key = ALIAS.get(mol, mol)
        if key not in gt:
            continue
        names.append(mol)
        x.append(float(sc[metric]))
        y.append(gt[key]["pIC50"])
    return names, np.asarray(x), np.asarray(y)


def bootstrap_rho(x, y, boot=BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(boot):
        i = rng.integers(0, len(y), len(y))
        if len(set(y[i])) < 3:
            continue
        r = stats.spearmanr(x[i], y[i])[0]
        if not np.isnan(r):
            vals.append(r)
    vals = np.asarray(vals)
    return np.percentile(vals, [2.5, 97.5]), float((vals <= GREY_RHO).mean())


def report(config="dimer_AB_v2"):
    here = os.path.dirname(__file__)
    with open(os.path.join(here, "data", "ground_truth.json")) as fh:
        gt = json.load(fh)
    with open(os.path.join(here, "data", "predictions.json")) as fh:
        pred = json.load(fh)
    with open(os.path.join(here, "data", "calib_smiles.json")) as fh:
        smi = json.load(fh)

    print("=" * 70)
    print(f"MODULE 3 — Boltz-2 metric calibration  [{config}]")
    print("=" * 70)

    per_mol = pred[config]
    metrics = sorted({m for v in per_mol.values() for m in v})

    # null-model bar on exactly the compounds that will be scored
    names0, _, y0 = _paired(per_mol, gt, "optimization_score")
    best_null, _ = m2_null_model.baseline([smi[n] for n in names0], y0)
    print(f"n = {len(names0)}   measured pIC50 "
          f"{min(y0):.2f}-{max(y0):.2f} ({max(y0) - min(y0):.1f} log span)")
    print(f"null-model bar |rho| = {best_null:.3f}\n")

    rows = []
    for metric in metrics:
        names, x, y = _paired(per_mol, gt, metric)
        if len(names) < 5:
            continue
        rho, p = stats.spearmanr(x, y)
        margin, null_verdict = m2_null_model.gate(rho, best_null)
        rows.append((metric, rho, p, margin, null_verdict))

    rows.sort(key=lambda r: -r[1])
    print(f"{'metric':<22}{'rho':>8}{'p':>9}{'vs null':>10}   verdict")
    for metric, rho, p, margin, nv in rows:
        flag = "" if margin > m2_null_model.MARGIN_REAL else "   <- " + nv
        print(f"{metric:<22}{rho:>+8.3f}{p:>9.4f}{margin:>+10.3f}   "
              f"{verdict(rho, p)}{flag}")

    # primary metric in detail
    names, x, y = _paired(per_mol, gt, "optimization_score")
    rho, p = stats.spearmanr(x, y)
    (lo, hi), p_fail = bootstrap_rho(x, y)
    print(f"\nprimary metric optimization_score: rho = {rho:+.3f}  p = {p:.4f}")
    print(f"  bootstrap 95% CI [{lo:+.3f}, {hi:+.3f}]   "
          f"P(rho <= {GREY_RHO}) = {p_fail:.4f}")

    keep = [i for i, n in enumerate(names) if n not in ARTIFACTS]
    r2, p2 = stats.spearmanr(x[keep], y[keep])
    print(f"  excluding {len(names) - len(keep)} thiol-reactive artifacts: "
          f"rho = {r2:+.3f}  p = {p2:.4f}  ({r2 - rho:+.3f})")

    order = np.argsort(-x)
    print("\n  predicted rank | measured pIC50:")
    for i in order:
        tag = "  [artifact label]" if names[i] in ARTIFACTS else ""
        print(f"    {names[i]:<18}{x[i]:>7.3f} | {y[i]:>5.2f}{tag}")

    # dimer vs monomer, on the metrics that actually discriminate
    if "monomer_A" in pred:
        print("\ndimer vs monomer (discriminative metrics only):")
        for metric in ("optimization_score", "binding_confidence"):
            _, xd, yd = _paired(pred[config], gt, metric)
            _, xm, ym = _paired(pred["monomer_A"], gt, metric)
            if len(xd) < 5 or len(xm) < 5:
                continue
            rd = stats.spearmanr(xd, yd)[0]
            rm = stats.spearmanr(xm, ym)[0]
            print(f"  {metric:<22} dimer {rd:+.3f}   monomer {rm:+.3f}   "
                  f"delta {rd - rm:+.3f}")
        print("  ! direction favours the dimer on both, but n is far too small"
              "\n    to establish it: the bootstrap CI on delta spans 0"
              "\n    (P(delta<=0) ~ 0.17). Structural biology requires the dimer;"
              "\n    this experiment does NOT demonstrate it.")

    # stereochemistry control
    ctrl = next((c for c in CONTROLS if c in per_mol), None)
    if ctrl and "nirmatrelvir" in per_mol:
        a = per_mol["nirmatrelvir"]["optimization_score"]
        b = per_mol[ctrl]["optimization_score"]
        print(f"\nstereochemistry control: nirmatrelvir {a:.3f} vs "
              f"gamma-lactam epimer {b:.3f}  (delta {a - b:+.3f})")
        print("  -> the metric is blind to this centre; not usable for stereo-SAR")

    return dict(rho=float(rho), p=float(p), ci=(float(lo), float(hi)),
                best_null=float(best_null))


if __name__ == "__main__":
    report()
