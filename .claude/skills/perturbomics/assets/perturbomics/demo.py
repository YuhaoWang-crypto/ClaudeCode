"""End-to-end demo: one representation, three modalities, one combination call.

Runs on synthetic-but-structured data (numpy/pandas/scipy only) so the whole
pipeline is reproducible without downloading L1000/CRISPR data. The biology is
*planted* so the answers are checkable:

  * disease pushes genes G0..G19 UP and G20..G39 DOWN
  * drug "full_corrector" reverses BOTH halves      (a coherent single reverser)
  * drug "agonist_like"   MIMICS the disease         (should score positive)
  * drug "hdac_like"      reverses only the UP half  (a partial agent)
  * drug "kinase_like"    reverses only the DOWN half (the complementary half)
  * CRISPRi_TARGETX / CRISPRa_TARGETY  partial genetic reversers of each half
  * plus random drugs / CRISPR guides as distractors

Two things this demonstrates, honestly:
  [1] single-agent WTCS cleanly ranks the coherent reverser (negative) and the
      mimic (positive); it deliberately ZEROS the one-tailed partial agents —
      that is the CMap "coherent connection" rule, not a bug.
  [3] combination analysis then RECOVERS those partial agents: hdac_like +
      kinase_like each look weak alone but together cover ~90% of the disease
      signature. That gap is exactly why combined analysis exists.

Run:  python3 -m perturbomics.demo
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .signature import Signature
from .connectivity import connectivity_matrix
from .combine import rank_reversers, best_combinations
from .pseudobulk import pseudobulk, signature_from_pseudobulk

RNG = np.random.default_rng(7)
GENES = [f"G{i}" for i in range(200)]


def _noisy(base: dict[str, float], scale: float = 0.3) -> pd.Series:
    s = pd.Series(0.0, index=GENES)
    for g, v in base.items():
        s[g] = v
    s += RNG.normal(0, scale, len(GENES))
    return s


def build_library() -> tuple[Signature, list[Signature]]:
    up = {g: 3.0 for g in GENES[0:20]}
    down = {g: -3.0 for g in GENES[20:40]}
    disease = Signature(_noisy({**up, **down}), name="IPF_fibrosis",
                        modality="disease")

    lib: list[Signature] = []
    # a coherent full reverser (fixes both halves) -> strong negative WTCS
    lib.append(Signature(_noisy({**{g: -3.0 for g in GENES[0:20]},
                                 **{g: 3.0 for g in GENES[20:40]}}),
                         name="full_corrector", modality="drug",
                         meta={"moa": "dual pathway", "target": "SMAD3"}))
    # complementary partial reversers (each fixes one half)
    lib.append(Signature(_noisy({g: -3.0 for g in GENES[0:20]}),
                         name="hdac_like", modality="drug",
                         meta={"moa": "HDAC inhibitor", "target": "HDAC1"}))
    lib.append(Signature(_noisy({g: 3.0 for g in GENES[20:40]}),
                         name="kinase_like", modality="drug",
                         meta={"moa": "kinase inhibitor", "target": "JAK2"}))
    # a mimic (makes disease worse) -> positive connectivity
    lib.append(Signature(_noisy({**up, **down}), name="agonist_like",
                         modality="drug", meta={"moa": "TGFb agonist"}))
    # CRISPR perturbations
    lib.append(Signature(_noisy({g: -3.0 for g in GENES[0:10]}),
                         name="CRISPRi_TARGETX", modality="crispr_i",
                         meta={"target": "TARGETX"}))
    lib.append(Signature(_noisy({g: 3.0 for g in GENES[30:40]}),
                         name="CRISPRa_TARGETY", modality="crispr_a",
                         meta={"target": "TARGETY"}))
    # distractors
    for i in range(6):
        hits = {g: RNG.choice([-3.0, 3.0]) for g in RNG.choice(GENES, 15,
                replace=False)}
        mod = "drug" if i % 2 else "crispr_ko"
        lib.append(Signature(_noisy(hits), name=f"decoy_{mod}_{i}",
                             modality=mod))
    return disease, lib


def synthetic_singlecell() -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    """A tiny single-cell perturbation dataset -> pseudobulk DGE demo."""
    genes = [f"G{i}" for i in range(60)]
    rows, samp, cond = [], [], []
    # 3 control + 3 perturbed biological replicates, 40 cells each
    up_genes = set(range(0, 8))
    for rep in range(6):
        treated = rep >= 3
        base = RNG.gamma(2.0, 1.0, len(genes))
        if treated:
            for gi in up_genes:
                base[gi] *= 8.0
        for c in range(40):
            lam = base * RNG.gamma(3.0, 1 / 3, len(genes))  # cell noise
            rows.append(RNG.poisson(lam))
            samp.append(f"rep{rep}")
            cond.append("perturbed" if treated else "control")
    counts = pd.DataFrame(rows, columns=genes)
    counts.index = [f"cell{i}" for i in range(len(counts))]
    sample = pd.Series(samp, index=counts.index)
    condition = pd.Series(cond, index=counts.index)
    return counts, sample, condition


def main() -> dict:
    print("=" * 68)
    print("perturbomics demo — drug + CRISPR + in-silico on one axis")
    print("=" * 68)

    disease, lib = build_library()
    print(f"\nDisease signature : {disease}")
    print(f"Library           : {len(lib)} perturbations "
          f"({sorted(set(s.modality for s in lib))})")

    # 1) single-agent reversers -------------------------------------------
    rev = rank_reversers(disease, lib, k=30)
    print("\n[1] Top single reversers (most negative connectivity = corrects "
          "disease):")
    print(rev.head(6).to_string(index=False,
          formatters={"connectivity": "{:+.3f}".format}))

    # 2) connectivity map among the top perturbations ---------------------
    top_names = list(rev.head(5)["name"])
    subset = [disease] + [s for s in lib if s.name in top_names]
    cmat = connectivity_matrix(subset, k=30)
    print("\n[2] Pairwise connectivity (WTCS; +mimic / -reverse):")
    print(cmat.round(2).to_string())

    # 3) best cross-modality combinations ---------------------------------
    combos = best_combinations(disease, lib, k_genes=40,
                               require_cross_modality=False, top_n=5)
    print("\n[3] Best combinations (coverage of disease signature reversed):")
    show = combos[["pert_a", "modality_a", "pert_b", "modality_b",
                   "coverage", "complementarity", "combo_score"]]
    print(show.to_string(index=False, float_format="{:.3f}".format))

    xmod = best_combinations(disease, lib, k_genes=40,
                             require_cross_modality=True, top_n=3)
    print("\n[3b] Best DRUG + CRISPR (cross-modality only) combination:")
    print(xmod[["pert_a", "modality_a", "pert_b", "modality_b",
                "coverage", "combo_score"]].head(1).to_string(
                index=False, float_format="{:.3f}".format))

    # 4) rigorous DGE: raw single cells -> pseudobulk -> signature ---------
    counts, sample, condition = synthetic_singlecell()
    pb = pseudobulk(counts, sample, min_cells=10)
    design = condition.groupby(sample).first()
    sig = signature_from_pseudobulk(pb, design, "perturbed", "control",
                                    name="scDGE_demo")
    top_up = sig.top(8, "up")
    planted = [f"G{i}" for i in range(8)]
    recovered = len(set(top_up) & set(planted))
    print(f"\n[4] Pseudobulk DGE ({sig.meta['engine']}): recovered "
          f"{recovered}/8 planted up-genes in the top-8. Signature: {sig}")

    print("\nRIGOR: scores/coverage are deterministic computations. Calling a "
          "combination 'synergistic' or a drug 'therapeutic' is a HYPOTHESIS "
          "to validate experimentally.")
    return {"reversers": rev, "connectivity": cmat, "combinations": combos,
            "dge_recovered": recovered}


if __name__ == "__main__":
    main()
