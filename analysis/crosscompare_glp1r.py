"""Cross-compare the evidence streams for GLP1R binders and quantify DrugCLIP
enrichment against measured ChEMBL actives.

Streams (all mined earlier into data/targets/GLP1R/):
  - ChEMBL measured actives  (GLP1R_ligands.csv)          -> ground truth
  - DrugCLIP predicted hits   (GLP1R_drugclip_predicted.csv)-> virtual screen
  - (humanPPI partners are proteins, a separate modality — reported, not fused here)

Method
------
Restrict ChEMBL to the small-molecule subset (MW < 900) so the comparison is
like-for-like with DrugCLIP's fragment-sized hits. For every DrugCLIP molecule
compute its nearest-neighbour Tanimoto (Morgan r2, 2048 bit) to any measured
active. Enrichment is tested against a DECOY set: DrugCLIP's predicted hits for
an unrelated target (EGFR, a kinase). If DrugCLIP is target-specific, GLP1R hits
should sit closer to GLP1R actives than EGFR hits do.

Outputs
-------
  data/targets/GLP1R/GLP1R_consensus.csv     ranked consensus molecules
  results/figures/glp1r_crosscompare.png     3-panel enrichment + consensus figure
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator, DataStructs
from scipy.stats import mannwhitneyu, spearmanr

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[1]
GLP1R = ROOT / "data" / "targets" / "GLP1R"
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# Okabe-Ito colourblind-safe palette
BLUE, ORANGE, GREY = "#0072B2", "#E69F00", "#8A8A8A"
_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fp(smiles: str):
    m = Chem.MolFromSmiles(str(smiles))
    return _GEN.GetFingerprint(m) if m else None


def mw(smiles: str):
    m = Chem.MolFromSmiles(str(smiles))
    return Descriptors.MolWt(m) if m else np.nan


def nn_tanimoto(query_fps, ref_fps) -> np.ndarray:
    """Max Tanimoto of each query fp to the set of reference fps."""
    out = np.zeros(len(query_fps))
    for i, q in enumerate(query_fps):
        out[i] = max(DataStructs.BulkTanimotoSimilarity(q, ref_fps)) if q and ref_fps else np.nan
    return out


def load_fps(df: pd.DataFrame):
    fps = [fp(s) for s in df["smiles"]]
    keep = [i for i, f in enumerate(fps) if f is not None]
    return df.iloc[keep].reset_index(drop=True), [fps[i] for i in keep]


def main() -> None:
    # --- load streams -----------------------------------------------------
    chembl = pd.read_csv(GLP1R / "GLP1R_ligands.csv")
    chembl_sm = chembl[chembl["class"] == "small_molecule"].copy()
    drug = pd.read_csv(GLP1R / "GLP1R_drugclip_predicted.csv").drop_duplicates("ligand_id")
    decoy = pd.read_csv(ROOT / "data" / "targets" / "EGFR_decoy" /
                        "EGFR_decoy_drugclip_predicted.csv").drop_duplicates("ligand_id")

    chembl_sm, chembl_fps = load_fps(chembl_sm)
    drug, drug_fps = load_fps(drug)
    decoy, decoy_fps = load_fps(decoy)
    print(f"ChEMBL small-molecule actives: {len(chembl_sm)}")
    print(f"DrugCLIP GLP1R hits: {len(drug)} | EGFR decoy hits: {len(decoy)}")

    # --- exact structural overlap (InChIKey) ------------------------------
    def ikey(s):
        m = Chem.MolFromSmiles(str(s))
        try:
            return Chem.MolToInchiKey(m) if m else None
        except Exception:
            return None
    ck = {ikey(s) for s in chembl_sm["smiles"]} - {None}
    dk = {ikey(s) for s in drug["smiles"]} - {None}
    exact = ck & dk
    conn = {k.split("-")[0] for k in ck} & {k.split("-")[0] for k in dk}
    print(f"exact InChIKey overlap: {len(exact)} | connectivity-layer overlap: {len(conn)}")

    # --- nearest measured-active similarity + MW --------------------------
    drug["nn_sim"] = nn_tanimoto(drug_fps, chembl_fps)
    decoy["nn_sim"] = nn_tanimoto(decoy_fps, chembl_fps)
    drug["mw"] = [mw(s) for s in drug["smiles"]]
    chembl_sm["mw"] = [mw(s) for s in chembl_sm["smiles"]]

    # --- enrichment stats -------------------------------------------------
    u, p = mannwhitneyu(drug["nn_sim"], decoy["nn_sim"], alternative="greater")
    thr = 0.35
    ef = (float((drug["nn_sim"] >= thr).mean()) /
          max(float((decoy["nn_sim"] >= thr).mean()), 1e-9))
    rho, prho = spearmanr(drug["drugclip_score"], drug["nn_sim"])
    print(f"\nENRICHMENT (GLP1R vs EGFR-decoy nearest-active similarity):")
    print(f"  median nn_sim  GLP1R={drug['nn_sim'].median():.3f}  decoy={decoy['nn_sim'].median():.3f}")
    print(f"  Mann-Whitney U p(one-sided) = {p:.2e}")
    print(f"  enrichment factor @Tanimoto>={thr}: {ef:.2f}x")
    print(f"  Spearman(drugclip_score, nn_sim) = {rho:.3f} (p={prho:.2g})")

    # --- consensus ranking ------------------------------------------------
    def z(x):
        x = np.asarray(x, float)
        s = x.std()
        return (x - x.mean()) / s if s > 0 else x * 0
    drug["combined_score"] = (z(drug["drugclip_score"]) + z(-drug["docking_score"])
                              + z(drug["nn_sim"]))
    drug["tier"] = np.where(drug["nn_sim"] >= 0.5, "corroborated_close",
                    np.where(drug["nn_sim"] >= thr, "corroborated_scaffold", "novel_chemotype"))
    drug = drug.sort_values("combined_score", ascending=False).reset_index(drop=True)
    cols = ["ligand_id", "smiles", "source", "drugclip_score", "docking_score",
            "nn_sim", "combined_score", "tier", "nearby_residues"]
    drug[cols].to_csv(GLP1R / "GLP1R_consensus.csv", index=False)
    print(f"\ntiers: {drug['tier'].value_counts().to_dict()}")
    print(f"wrote consensus -> {GLP1R/'GLP1R_consensus.csv'}")

    # --- figure -----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 130})
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.6))

    # Panel A: enrichment — nearest-active similarity distributions
    bins = np.linspace(0, 1, 26)
    ax[0].hist(decoy["nn_sim"], bins=bins, density=True, color=GREY, alpha=0.75,
               label=f"EGFR decoy (n={len(decoy)})")
    ax[0].hist(drug["nn_sim"], bins=bins, density=True, color=BLUE, alpha=0.75,
               label=f"GLP1R DrugCLIP (n={len(drug)})")
    ax[0].axvline(decoy["nn_sim"].median(), color=GREY, ls="--", lw=2)
    ax[0].axvline(drug["nn_sim"].median(), color=BLUE, ls="--", lw=2)
    ax[0].set_xlabel("Nearest Tanimoto to a measured GLP1R active")
    ax[0].set_ylabel("Density")
    ax[0].set_title(f"A  DrugCLIP enrichment\nEF@0.35={ef:.1f}×, MWU p={p:.1e}")
    ax[0].legend(frameon=False, fontsize=8.5)

    # Panel B: MW distributions — explains the chemical-space gap
    mbins = np.linspace(100, 950, 30)
    ax[1].hist(chembl_sm["mw"].dropna(), bins=mbins, density=True, color=ORANGE,
               alpha=0.75, label=f"ChEMBL actives (med {chembl_sm['mw'].median():.0f})")
    ax[1].hist(drug["mw"].dropna(), bins=mbins, density=True, color=BLUE,
               alpha=0.75, label=f"DrugCLIP hits (med {drug['mw'].median():.0f})")
    ax[1].set_xlabel("Molecular weight (Da)")
    ax[1].set_ylabel("Density")
    ax[1].set_title("B  Why: fragment vs drug-sized\nDrugCLIP hits are fragments")
    ax[1].legend(frameon=False, fontsize=8.5)

    # Panel C: DrugCLIP-internal ranking (no measured-active corroboration exists)
    for src, col in [("REAL", BLUE), ("ZINC", ORANGE)]:
        s = drug[drug["source"] == src]
        ax[2].scatter(s["drugclip_score"], -s["docking_score"], s=28, color=col,
                      edgecolor="white", linewidth=0.5, label=f"{src} (n={len(s)})")
    top = drug.head(5)
    for _, r in top.iterrows():
        ax[2].annotate(r["ligand_id"], (r["drugclip_score"], -r["docking_score"]),
                       fontsize=6.5, color="#333", xytext=(3, 3),
                       textcoords="offset points")
    ax[2].set_xlabel("DrugCLIP retrieval score")
    ax[2].set_ylabel("Docking score (−kcal/mol, higher=better)")
    ax[2].set_title("C  Best DrugCLIP-internal picks\n(fragment starting points only)")
    ax[2].legend(frameon=False, fontsize=8.5)

    fig.tight_layout()
    out = FIG / "glp1r_crosscompare.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote figure -> {out}")


if __name__ == "__main__":
    main()
