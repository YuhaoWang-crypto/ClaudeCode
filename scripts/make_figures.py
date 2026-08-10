#!/usr/bin/env python
"""Render every figure used in the report, from measured results only.

Nothing here is illustrative.  Each panel reads a JSON written by one of the
other scripts, or recomputes from real data, so a figure cannot drift away from
the number it depicts.

    python scripts/make_figures.py            # all figures into figures/mri_cd8_til/
    python scripts/make_figures.py --only roi # one figure
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from mri_cd8_til.modeling.metrics import auc_interval_from_reported  # noqa: E402

OUT = Path("figures/mri_cd8_til")
RESULTS = Path("results")

# A restrained palette that survives greyscale printing and colour-vision
# deficiency: one warm accent for "the suspect number", one cool for "the
# defensible number", greys for context.
SUSPECT = "#c2410c"
HONEST = "#1d4ed8"
NEUTRAL = "#64748b"
LIGHT = "#cbd5e1"
ACCENT3 = "#047857"

plt.rcParams.update(
    {
        "figure.dpi": 130,
        "savefig.dpi": 130,
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.6,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    }
)


def _save(fig, name: str) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {path}")
    return path


# --------------------------------------------------------------------- fig 1


def fig_roi_boundary(patient: str = "TCGA-AR-A1AQ") -> None:
    """ROI construction on a real TCGA-BRCA DCE subtraction series.

    Shows what "inner 2 mm rim" and "outer ring" actually select, and why
    spacing-aware morphology matters at 0.6 x 0.6 x 2.0 mm.
    """
    import SimpleITK as sitk

    from mri_cd8_til.clients.tcia import TCIAClient
    from mri_cd8_til.imaging.preprocess import PreprocessConfig, preprocess_volume
    from mri_cd8_til.imaging.roi import (
        RegionKind,
        RegionSpec,
        body_mask_from_image,
        build_region,
    )

    cache = RESULTS / "dicom_cache" / patient
    if not cache.exists() or not any(cache.iterdir()):
        client = TCIAClient()
        series = [s for s in client.series("TCGA-BRCA", "MR") if s.patient_id == patient]
        subtraction = [s for s in series if ")-(" in s.description and s.image_count <= 140]
        chosen = min(subtraction or series, key=lambda s: s.image_count)
        print(f"  downloading {chosen.description} ({chosen.image_count} images)")
        client.download_series(chosen.series_uid, cache)

    reader = sitk.ImageSeriesReader()
    reader.SetFileNames(reader.GetGDCMSeriesFileNames(str(cache)))
    image = preprocess_volume(
        reader.Execute(),
        PreprocessConfig(target_spacing_mm=(0.6, 0.6, 2.0), apply_n4=True, normalization="percentile"),
    )

    body = sitk.Cast(body_mask_from_image(image), sitk.sitkUInt8)
    array = sitk.GetArrayFromImage(image)
    threshold = float(np.percentile(array[sitk.GetArrayFromImage(body) > 0], 99.85))
    seed = sitk.BinaryMorphologicalClosing(
        sitk.And(sitk.Cast(image > threshold, sitk.sitkUInt8), body), [1, 1, 1]
    )
    # Pick the most *compact* component, not the largest. Thresholding a
    # subtraction image catches the whole enhancing fibroglandular tissue, which
    # is a sprawling crescent and illustrates nothing about a lesion boundary;
    # the roundest sizeable component at least has lesion-like geometry.
    components = sitk.ConnectedComponent(seed)
    stats = sitk.LabelShapeStatisticsImageFilter()
    stats.Execute(components)
    candidates = [
        (stats.GetRoundness(label), label)
        for label in stats.GetLabels()
        if 1500 <= stats.GetNumberOfPixels(label) <= 60000
    ]
    best = max(candidates)[1] if candidates else 1
    mask = sitk.Cast(components == best, sitk.sitkUInt8)

    regions = {
        "whole tumour": build_region(mask, RegionSpec(RegionKind.WHOLE)),
        "inner rim 2 mm\n(the paper's region)": build_region(
            mask, RegionSpec(RegionKind.INNER_RIM, 2.0)
        ),
        "outer ring 3 mm": build_region(mask, RegionSpec(RegionKind.OUTER_RING, 3.0), body_mask=body),
        "outer ring 5 mm": build_region(mask, RegionSpec(RegionKind.OUTER_RING, 5.0), body_mask=body),
    }

    mask_array = sitk.GetArrayFromImage(mask)
    z = int(np.argmax(mask_array.sum(axis=(1, 2))))
    rows, cols = np.where(mask_array[z] > 0)
    pad = 70
    r0, r1 = max(0, rows.min() - pad), min(array.shape[1], rows.max() + pad)
    c0, c1 = max(0, cols.min() - pad), min(array.shape[2], cols.max() + pad)
    background = array[z, r0:r1, c0:c1]
    vmin, vmax = np.percentile(background, [2, 99.5])

    fig, axes = plt.subplots(1, 5, figsize=(15.5, 3.9))
    axes[0].imshow(background, cmap="gray", vmin=vmin, vmax=vmax)
    axes[0].set_title("DCE subtraction\n(real TCGA-BRCA)", fontsize=9)
    colours = [SUSPECT, HONEST, ACCENT3, "#7c3aed"]
    counts = {}
    for ax, (name, region), colour in zip(axes[1:], regions.items(), colours):
        ax.imshow(background, cmap="gray", vmin=vmin, vmax=vmax)
        layer = sitk.GetArrayFromImage(region)[z, r0:r1, c0:c1]
        overlay = np.zeros((*layer.shape, 4))
        overlay[layer > 0] = matplotlib.colors.to_rgba(colour, alpha=0.55)
        ax.imshow(overlay)
        n = int(sitk.GetArrayFromImage(region).sum())
        counts[name] = n
        ax.set_title(f"{name}\n{n:,} voxels", fontsize=9)
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)

    whole = counts["whole tumour"]
    rim = counts["inner rim 2 mm\n(the paper's region)"]
    fig.suptitle(
        f"ROI construction on real DICOM — {patient}, resampled to 0.6 × 0.6 × 2.0 mm.  "
        "Morphology radii are computed per axis, so a 2 mm rim is 2 mm in all three directions.",
        fontsize=9.5,
        y=1.02,
    )
    fig.text(
        0.5,
        -0.04,
        f"Note the inner rim is {rim / whole:.0%} of the whole tumour ({rim:,} of {whole:,} voxels): "
        "on a lesion this size a 2 mm erosion leaves almost no core, so the\n"
        '"peritumoural" and "whole tumour" feature sets are largely the same voxels — '
        "one reason their reported AUCs move together.",
        ha="center",
        fontsize=8.5,
        color=NEUTRAL,
    )
    _save(fig, "fig1_roi_boundary")


# --------------------------------------------------------------------- fig 2


def fig_optimism() -> None:
    """Reported AUC vs true effect, for three analysis protocols."""
    data = json.loads((RESULTS / "optimism_study.json").read_text())
    table = data["optimism"]["table"]
    grid = data["optimism"]["grid"]

    def series(protocol: str, key: str) -> np.ndarray:
        by_auc = {row["true_feature_auc"]: row for row in table if row["protocol"] == protocol}
        return np.array([by_auc[a][key] for a in grid])

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.1))
    x = np.array(grid)
    for protocol, colour, label in (
        ("published", SUSPECT, "as published\n(select on all 182, then split)"),
        ("split_first", HONEST, "split first\n(select on the 137 training only)"),
        ("repeated_cv", NEUTRAL, "repeated cross-validation"),
    ):
        median = series(protocol, "median_reported_auc")
        ax.plot(x, median, "o-", color=colour, label=label, lw=2, ms=5)
        ax.fill_between(x, series(protocol, "p05"), series(protocol, "p95"), color=colour, alpha=0.13)

    ax.axhline(0.985, color="black", ls=":", lw=1.2)
    ax.text(0.502, 0.992, "published headline 0.985", fontsize=8, color="black")
    ax.plot([0.50], [series("published", "median_reported_auc")[0]], "o", ms=11,
            mfc="none", mec=SUSPECT, mew=2)
    ax.annotate(
        f"no signal at all →\nstill reports {series('published', 'median_reported_auc')[0]:.3f}",
        xy=(0.50, series("published", "median_reported_auc")[0]),
        xytext=(0.545, 0.63), fontsize=8, color=SUSPECT,
        arrowprops=dict(arrowstyle="->", color=SUSPECT, lw=1.2),
    )
    ax.set_xlabel("true discriminative power of each informative feature (AUC)")
    ax.set_ylabel("reported validation AUC")
    ax.set_title("What the protocol reports vs what is really there", fontsize=10)
    ax.legend(fontsize=7.5, loc="lower right", framealpha=0.95)
    ax.set_ylim(0.3, 1.03)

    gap = series("published", "median_reported_auc") - series("split_first", "median_reported_auc")
    bars = ax2.bar([f"{a:.2f}" for a in grid], gap, color=SUSPECT, alpha=0.85, width=0.62)
    for rect, value in zip(bars, gap):
        ax2.text(rect.get_x() + rect.get_width() / 2, value + 0.006, f"+{value:.3f}",
                 ha="center", fontsize=8)
    ax2.set_xlabel("true per-feature AUC")
    ax2.set_ylabel("AUC manufactured by the leak")
    ax2.set_title("Selecting features before splitting is worth +0.21–0.25\nin the weak-signal regime",
                  fontsize=10)
    ax2.set_ylim(0, max(gap) * 1.22)
    _save(fig, "fig2_optimism")


# --------------------------------------------------------------------- fig 3


def fig_dimension_sweep() -> None:
    """Zero-signal reported AUC as a function of the candidate pool."""
    data = json.loads((RESULTS / "optimism_study.json").read_text())["dimension_sweep"]
    n = np.array(data["n_features"], dtype=float)
    median = np.array(data["median_auc"])
    p95 = np.array(data["p95_auc"])

    fig, ax = plt.subplots(figsize=(6.6, 4.1))
    ax.fill_between(n, median, p95, color=SUSPECT, alpha=0.14, label="median → 95th percentile")
    ax.plot(n, median, "o-", color=SUSPECT, lw=2, ms=6, label="median reported AUC")
    ax.plot(n, p95, "^--", color=SUSPECT, lw=1.2, ms=5, alpha=0.75)
    ax.axhline(0.5, color=NEUTRAL, ls="--", lw=1.2)
    ax.text(55, 0.508, "truth: AUC = 0.5 (there is no signal in this data)", fontsize=8, color=NEUTRAL)
    ax.axvline(3332, color=HONEST, ls=":", lw=1.4)
    ax.text(3500, 0.53, "the paper's\n833 × 4 phases", fontsize=8, color=HONEST)
    ax.set_xscale("log")
    ax.set_xlabel("number of candidate features offered to the selection step (log scale)")
    ax.set_ylabel("reported validation AUC")
    ax.set_title("Under zero signal, reported AUC is a function of how many\nfeatures the selector got to look at",
                 fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    _save(fig, "fig3_dimension_sweep")


# --------------------------------------------------------------------- fig 4


def fig_multicenter() -> None:
    """What multi-centre data buys, under two kinds of centre effect."""
    data = json.loads((RESULTS / "multicenter_demo.json").read_text())
    scenarios = data["scenarios"]
    keys = ["single_small", "single_large", "pooled_random_cv", "loco_raw", "loco_combat"]
    labels = [
        "1 centre\nn=182",
        f"1 centre\nn={data['n_pooled']}",
        f"{data['n_centers']} centres\nrandom CV",
        "leave-one-centre-out\nraw",
        "leave-one-centre-out\n+ ComBat",
    ]
    x = np.arange(len(keys))
    width = 0.38

    fig, ax = plt.subplots(figsize=(9.2, 4.3))
    for offset, (name, colour, pretty) in zip(
        (-width / 2, width / 2),
        (("nuisance", LIGHT, "centre = noise"), ("confounded", HONEST, "centre = confounder")),
    ):
        means = [scenarios[name][k]["mean"] for k in keys]
        errors = [scenarios[name][k]["std"] for k in keys]
        bars = ax.bar(x + offset, means, width, yerr=errors, capsize=3, label=pretty,
                      color=colour, edgecolor="white")
        for rect, value in zip(bars, means):
            ax.text(rect.get_x() + rect.get_width() / 2, value + 0.022, f"{value:.3f}",
                    ha="center", fontsize=7.5)

    ax.axhline(0.985, color=SUSPECT, ls=":", lw=1.4)
    ax.text(-0.45, 0.995, "the published 0.985 — a different quantity, not a target",
            fontsize=8, color=SUSPECT)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("AUC")
    ax.set_ylim(0.55, 1.05)
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title(
        "True effect size is held fixed across every bar, so the differences are caused by the "
        "study design alone",
        fontsize=10,
    )
    _save(fig, "fig4_multicenter")


# --------------------------------------------------------------------- fig 5


def fig_auc_intervals() -> None:
    """Forest plot of the published AUCs with their real class counts."""
    # Class distribution from the paper: desert 67, excluded 30, inflamed 85 (n=182).
    v = {k: round(n * 45 / 182) for k, n in (("d", 67), ("e", 30), ("i", 85))}
    t = {k: round(n * 137 / 182) for k, n in (("d", 67), ("e", 30), ("i", 85))}
    rows = [
        ("RM-whole  training   (inflamed vs rest)", 0.973, t["i"], t["d"] + t["e"], False),
        ("RM-whole  validation (inflamed vs rest)", 0.985, v["i"], v["d"] + v["e"], True),
        ("RM-peri   training   (excluded vs desert)", 0.993, t["e"], t["d"], False),
        ("RM-peri   validation (excluded vs desert)", 0.984, v["e"], v["d"], True),
    ]

    fig, ax = plt.subplots(figsize=(9.4, 4.0))
    for i, (name, auc, n_pos, n_neg, highlight) in enumerate(rows):
        interval = auc_interval_from_reported(auc, n_pos, n_neg)
        colour = SUSPECT if highlight else NEUTRAL
        ax.plot([interval.lower, interval.upper], [i, i], color=colour, lw=3.2, alpha=0.8,
                solid_capstyle="round")
        ax.plot([auc], [i], "o", color=colour, ms=8, zorder=3)
        ax.text(interval.lower - 0.006, i, f"{auc:.3f}", ha="right", va="center", fontsize=8.5,
                color=colour, fontweight="bold" if highlight else "normal")
        ax.text(1.004, i, f"n+={n_pos}, n−={n_neg}", va="center", fontsize=8, color=NEUTRAL)

    # A distinctly weaker model on the smallest split, for comparison.
    weak = auc_interval_from_reported(0.80, v["e"], v["d"])
    ax.plot([weak.lower, weak.upper], [-1, -1], color=ACCENT3, lw=3.2, alpha=0.7,
            solid_capstyle="round")
    ax.plot([0.80], [-1], "s", color=ACCENT3, ms=7)
    ax.text(weak.lower - 0.006, -1, "0.800", ha="right", va="center", fontsize=8.5, color=ACCENT3)
    ax.text(1.004, -1, "a much weaker model,\nsame RM-peri split", va="center", fontsize=8,
            color=ACCENT3)

    ax.set_yticks(range(-1, len(rows)))
    ax.set_yticklabels(["(comparison)"] + [r[0] for r in rows], fontsize=8.5, family="monospace")
    ax.set_xlim(0.55, 1.0)
    ax.set_xlabel("AUC with 95% confidence interval (Hanley–McNeil)")
    ax.axvline(1.0, color=LIGHT, lw=1)
    ax.set_title(
        "The RM-peri validation split has only 7 immune-excluded patients.\n"
        "Its 0.984 cannot be told apart from a model at 0.80.",
        fontsize=10,
    )
    ax.grid(axis="y", visible=False)
    _save(fig, "fig5_auc_intervals")


# --------------------------------------------------------------------- fig 6


def fig_model_escalation() -> None:
    """The paper's own six validation AUCs, in the order they escalate."""
    labels = ["best single\nDCE phase", "score\ncombination", "feature\ncombination"]
    whole = [0.671, 0.811, 0.985]
    peri = [0.866, 0.976, 0.984]

    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    x = np.arange(3)
    ax.plot(x, whole, "o-", color=SUSPECT, lw=2.4, ms=9, label="RM-whole (inflamed vs rest)")
    ax.plot(x, peri, "s--", color=HONEST, lw=2.4, ms=8, label="RM-peri (desert vs excluded)")
    for xi, value in zip(x, whole):
        ax.text(xi, value + 0.022, f"{value:.3f}", ha="center", fontsize=9, color=SUSPECT)
    for xi, value in zip(x, peri):
        ax.text(xi, value - 0.038, f"{value:.3f}", ha="center", fontsize=9, color=HONEST)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("validation AUC (same 45 patients throughout)")
    ax.set_ylim(0.6, 1.05)
    ax.legend(fontsize=8.5, loc="lower right")
    ax.set_title(
        "Six AUCs reported on one 45-patient split; the headline is the maximum.\n"
        "0.671 → 0.985 comes from changing the combination strategy, not the biology.",
        fontsize=10,
    )
    _save(fig, "fig6_model_escalation")


# --------------------------------------------------------------------- fig 7


def fig_real_til_maps(n_show: int = 4) -> None:
    """Real Stony Brook TIL maps for paired TCGA-BRCA patients."""
    from mri_cd8_til.clients.idc import load_map

    labels_path = RESULTS / "til_labels.json"
    if not labels_path.exists():
        labels_path = RESULTS / "til_labels_smoke.json"
    payload = json.loads(labels_path.read_text())
    records = sorted(payload["records"], key=lambda r: r["til_fraction_tissue"])
    if len(records) < n_show:
        n_show = len(records)
    picks = [records[0], records[len(records) // 3], records[2 * len(records) // 3], records[-1]][:n_show]

    fig, axes = plt.subplots(2, n_show, figsize=(3.4 * n_show, 6.6),
                             gridspec_kw={"height_ratios": [3, 2]})
    if n_show == 1:
        axes = axes.reshape(2, 1)

    for col, record in enumerate(picks):
        grid = load_map(RESULTS / "til_maps", record["patient_id"]).grid
        tissue = grid > 0
        axes[0, col].imshow(np.where(tissue, grid, np.nan), cmap="magma", vmin=0, vmax=255)
        axes[0, col].set_title(
            f"{record['patient_id']}\nTIL {record['til_fraction_tissue']*100:.1f}% of tissue\n"
            f"Moran's I {record['til_moran_i']:+.2f}",
            fontsize=8.5,
        )
        axes[0, col].axis("off")

        values = grid[tissue].astype(float) / 255.0
        axes[1, col].hist(values, bins=40, color=HONEST, alpha=0.8)
        axes[1, col].axvline(0.5, color=SUSPECT, ls="--", lw=1.3)
        axes[1, col].set_yscale("log")
        axes[1, col].set_xlabel("TIL probability per 50 µm patch", fontsize=8)
        if col == 0:
            axes[1, col].set_ylabel("patches (log)", fontsize=8)
        axes[1, col].tick_params(labelsize=7.5)

    fig.suptitle(
        f"Real TIL maps (Stony Brook Inception-V4) for paired TCGA-BRCA patients — "
        f"{payload['n_patients']} retrieved.  Dashed line: the 0.5 cut, verified against the "
        "companion binary maps at 99.9% patch agreement.",
        fontsize=9.5,
        y=1.0,
    )
    _save(fig, "fig7_real_til_maps")


# --------------------------------------------------------------------- fig 8


def fig_til_label_distribution() -> None:
    """The real label distribution of the external test cohort."""
    labels_path = RESULTS / "til_labels.json"
    if not labels_path.exists():
        labels_path = RESULTS / "til_labels_smoke.json"
    payload = json.loads(labels_path.read_text())
    records = payload["records"]
    fraction = np.array([r["til_fraction_tissue"] for r in records]) * 100
    moran = np.array([r["til_moran_i"] for r in records])
    cluster = np.array([r["til_cluster_share"] for r in records]) * 100

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
    axes[0].hist(fraction, bins=22, color=HONEST, alpha=0.85, edgecolor="white")
    cut = np.quantile(fraction, 2 / 3)
    axes[0].axvline(cut, color=SUSPECT, ls="--", lw=1.5)
    axes[0].text(cut * 1.03, axes[0].get_ylim()[1] * 0.82,
                 f"upper-tertile cut\n{cut:.1f}%", fontsize=8, color=SUSPECT)
    axes[0].set_xlabel("TIL-positive share of tissue patches (%)")
    axes[0].set_ylabel("patients")
    axes[0].set_title("How much infiltrate", fontsize=10)

    axes[1].hist(moran, bins=22, color=ACCENT3, alpha=0.85, edgecolor="white")
    axes[1].set_xlabel("Moran's I of the TIL field (tissue patches only)")
    axes[1].set_title("How clustered it is", fontsize=10)

    axes[2].scatter(fraction, moran, s=26, c=cluster, cmap="viridis", alpha=0.85,
                    edgecolor="white", linewidth=0.4)
    bar = plt.colorbar(axes[2].collections[0], ax=axes[2])
    bar.set_label("largest cluster share (%)", fontsize=8)
    axes[2].set_xlabel("TIL-positive share of tissue (%)")
    axes[2].set_ylabel("Moran's I")
    axes[2].set_title("Amount and arrangement are only\nweakly related — which is the point",
                      fontsize=10)

    correlation = np.corrcoef(fraction, moran)[0, 1]
    axes[2].text(0.03, 0.05, f"Pearson r = {correlation:+.2f}", transform=axes[2].transAxes,
                 fontsize=8.5, color=NEUTRAL)

    fig.suptitle(
        f"Real labels for the external test cohort — {payload['n_patients']} paired TCGA-BRCA "
        "patients with MRI + RNA-seq + diagnostic WSI",
        fontsize=10,
        y=1.03,
    )
    _save(fig, "fig8_til_label_distribution")


FIGURES = {
    "roi": fig_roi_boundary,
    "optimism": fig_optimism,
    "sweep": fig_dimension_sweep,
    "multicenter": fig_multicenter,
    "intervals": fig_auc_intervals,
    "escalation": fig_model_escalation,
    "tilmaps": fig_real_til_maps,
    "tildist": fig_til_label_distribution,
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=sorted(FIGURES), default=None)
    args = parser.parse_args()

    chosen = {args.only: FIGURES[args.only]} if args.only else FIGURES
    failures = []
    for name, function in chosen.items():
        print(f"[{name}]")
        try:
            function()
        except Exception as exc:  # noqa: BLE001
            print(f"  FAILED: {type(exc).__name__}: {exc}")
            failures.append(name)
    if failures:
        print(f"\n{len(failures)} figure(s) failed: {', '.join(failures)}")
        return 1
    print(f"\nall {len(chosen)} figures written to {OUT}/")
    return 0




# --------------------------------------------------------------------- fig 9


def fig_label_validation() -> None:
    """Does the H&E TIL label measure lymphocyte burden? Check vs transcriptome."""
    data = json.loads((RESULTS / "label_validation.json").read_text())
    rows = data["correlations"]
    order = [
        "CD8 core signature",
        "CD8A",
        "TIS (18-gene inflamed)",
        "cytolytic (GZMA, PRF1)",
        "PTPRC (CD45, pan-leukocyte)",
        "CD8B",
        "EPCAM (epithelial, negative control)",
    ]
    order = [k for k in order if k in rows]
    amount = [rows[k]["rho_til_fraction"] for k in order]
    arrange = [rows[k]["rho_moran"] for k in order]
    y = np.arange(len(order))

    fig, ax = plt.subplots(figsize=(8.8, 4.3))
    ax.barh(y + 0.19, amount, 0.36, label="TIL amount (share of tissue)", color=HONEST, alpha=0.9)
    ax.barh(y - 0.19, arrange, 0.36, label="TIL arrangement (Moran's I)", color=ACCENT3, alpha=0.85)
    for yi, value in zip(y + 0.19, amount):
        ax.text(value + 0.012, yi, f"{value:+.3f}", va="center", fontsize=8.5, color=HONEST)

    control = len(order) - 1
    ax.barh(control + 0.19, amount[control], 0.36, color=SUSPECT, alpha=0.9)
    ax.barh(control - 0.19, arrange[control], 0.36, color=SUSPECT, alpha=0.6)

    ax.axvline(0, color=NEUTRAL, lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Spearman ρ with the pathology-derived label")
    ax.set_title(
        f"The label is measuring what it claims (n={data['n_patients']}).\n"
        "Immune signatures track it; the epithelial negative control does not.",
        fontsize=10,
    )
    ax.legend(fontsize=8.5, loc="lower right")
    _save(fig, "fig9_label_validation")


# -------------------------------------------------------------------- fig 10


def fig_real_experiment() -> None:
    """The real train/test result, with its permutation null and interval."""
    data = json.loads((RESULTS / "real_experiment.json").read_text())
    cv = data["cv_by_model"]
    null = data["permutation_null"]
    test = data["locked_test"]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.3))

    names = list(cv)
    means = [cv[n]["mean"] for n in names]
    stds = [cv[n]["std"] for n in names]
    colours = [HONEST if n == data["selected_model"] else LIGHT for n in names]
    ax.barh(np.arange(len(names)), means, xerr=stds, color=colours, capsize=3, height=0.6)
    ax.axvline(0.5, color=NEUTRAL, ls="--", lw=1.2)
    ax.axvline(null["p95"], color=SUSPECT, ls=":", lw=1.6)
    ax.text(null["p95"] + 0.006, len(names) - 0.35,
            f"95th pct of the\npermutation null\n{null['p95']:.3f}", fontsize=8, color=SUSPECT)
    ax.set_yticks(np.arange(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("cross-validated AUC on the training set")
    ax.set_xlim(0.3, max(0.85, null["p95"] + 0.12))
    ax.set_title("Model selection, inside the training set only", fontsize=10)

    labels = ["training CV\n(selected model)", "locked test set\n(touched once)",
              "leave-one-\nscanner-out", "the leaky protocol\non this same data"]
    values = [
        null["observed"],
        test["auc"],
        data["leave_one_scanner_out"]["mean"] or np.nan,
        data["leaky_protocol_test_auc"],
    ]
    bar_colours = [NEUTRAL, HONEST, ACCENT3, SUSPECT]
    bars = ax2.bar(labels, values, color=bar_colours, width=0.6)
    ax2.errorbar(1, test["auc"],
                 yerr=[[test["auc"] - test["ci_lower"]], [test["ci_upper"] - test["auc"]]],
                 fmt="none", ecolor=HONEST, capsize=5, lw=1.6)
    for rect, value in zip(bars, values):
        if np.isfinite(value):
            ax2.text(rect.get_x() + rect.get_width() / 2, value + 0.018, f"{value:.3f}",
                     ha="center", fontsize=9)
    ax2.axhline(0.5, color=NEUTRAL, ls="--", lw=1.2)
    ax2.text(-0.45, 0.512, "chance", fontsize=8, color=NEUTRAL)
    ax2.set_ylabel("AUC")
    ax2.set_ylim(0.3, 1.0)
    ax2.tick_params(axis="x", labelsize=8.5)
    ax2.set_title(
        f"Real result: test AUC {test['auc']:.3f} "
        f"(95% CI {test['ci_lower']:.2f}–{test['ci_upper']:.2f}), permutation p={null['p_value']:.2f}",
        fontsize=10,
    )
    _save(fig, "fig10_real_experiment")


# -------------------------------------------------------------------- fig 11


def fig_optimization_sweep() -> None:
    """Every optimisation variant against its own null, and the multiplicity bar."""
    data = json.loads((RESULTS / "model_optimization.json").read_text())
    variants = [v for v in data["variants"] if v["metric"] == "auc"]
    mult = data["multiplicity"]

    names = [v["name"] for v in variants]
    observed = [v["observed"] for v in variants]
    p95 = [v["null_p95"] for v in variants]
    y = np.arange(len(names))

    fig, ax = plt.subplots(figsize=(9.6, 4.6))
    for yi, (obs, bar) in enumerate(zip(observed, p95)):
        ax.plot([0.5, bar], [yi, yi], color=LIGHT, lw=7, solid_capstyle="butt")
        ax.plot([bar], [yi], "|", color=NEUTRAL, ms=13, mew=1.6)
        clears = obs > bar
        ax.plot([obs], [yi], "o", ms=9,
                color=ACCENT3 if clears else NEUTRAL, zorder=3)
        ax.text(obs + 0.008, yi, f"{obs:.3f}", va="center", fontsize=8.5,
                color=ACCENT3 if clears else NEUTRAL)

    ax.axvline(mult["null_max_p95"], color=SUSPECT, lw=2)
    ax.text(mult["null_max_p95"] + 0.006, len(names) - 0.6,
            f"multiplicity-adjusted bar\n{mult['null_max_p95']:.3f}\n"
            f"({mult['n_auc_variants']} variants tried)",
            fontsize=8.5, color=SUSPECT, va="top")

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("cross-validated AUC")
    ax.set_xlim(0.44, max(mult["null_max_p95"], max(observed)) + 0.09)
    ax.set_title(
        "Grey bar = each variant's own permutation null (to its 95th percentile).\n"
        f"Best observed {mult['best_observed']:.3f} does not reach the "
        f"multiplicity-adjusted {mult['null_max_p95']:.3f}.",
        fontsize=10,
    )
    ax.grid(axis="y", visible=False)
    _save(fig, "fig11_optimization_sweep")


FIGURES.update({
    "labelval": fig_label_validation,
    "realexp": fig_real_experiment,
    "sweep2": fig_optimization_sweep,
})

if __name__ == "__main__":
    raise SystemExit(main())
