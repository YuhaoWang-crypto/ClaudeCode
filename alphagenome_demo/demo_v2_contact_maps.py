"""AlphaGenome v2 demo: CONTACT_MAPS (predicted 3D chromatin contacts).

The CONTACT_MAPS output is a (512, 512, tracks) matrix for a 1 Mb window at
2,048 bp resolution: predicted Hi-C contact frequency between every pair of
bins, as observed-over-expected.

This script:
  1. predicts the contact map and CTCF ChIP-seq for a 1 Mb chr22 window in
     transverse colon,
  2. calls domain boundaries from an insulation score (in skeletal muscle,
     where the predicted map has clear domains) and checks whether they
     coincide with predicted CTCF peaks,
  3. compares the map across three tissues,
  4. scrambles the CTCF site at the strongest boundary (predict_sequence) and
     measures how the predicted contacts change.

Requires ALPHAGENOME_API_KEY.

    python alphagenome_demo/demo_v2_contact_maps.py
"""

import json
import os
import pathlib
import urllib.request

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client

OUT_DIR = pathlib.Path(__file__).parent
pd.set_option("display.width", 200)
OT = dna_client.OutputType

interval = genome.Interval("chr22", 35677410, 36725986)
TISSUES = {
    "transverse colon": "UBERON:0001157",
    "gastrocnemius medialis": "UBERON:0011907",
    "dorsolateral prefrontal cortex": "UBERON:0009834",
}
PRIMARY = "gastrocnemius medialis"  # tissue used for boundary calling and edits
WINDOW = 15  # insulation window in bins (~31 kb each side)

model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)


def ctcf_profile(chip_tf, n_bins):
    """Mean of the observed CTCF tracks, re-binned to the contact-map bins."""
    md = chip_tf.metadata
    sel = ((md.transcription_factor == "CTCF") & (md.imputation_mode == "observed")).to_numpy()
    x = chip_tf.values[:, sel].mean(axis=1)
    return x.reshape(n_bins, -1).mean(axis=1), int(sel.sum())


def insulation(mat, w=WINDOW):
    """log2 mean contact across each bin boundary (low = insulating)."""
    n = mat.shape[0]
    score = np.full(n, np.nan)
    for i in range(w, n - w):
        score[i] = mat[i - w : i, i : i + w].mean()
    return np.log2(score / np.nanmean(score))


def boundaries(ins, min_depth=0.1, sep=10):
    """Local minima of the insulation score, strongest first."""
    idx = [
        i
        for i in range(1, len(ins) - 1)
        if np.isfinite(ins[i - 1 : i + 2]).all() and ins[i] < ins[i - 1] and ins[i] <= ins[i + 1]
    ]
    kept = []
    for i in sorted(idx, key=lambda i: ins[i]):
        lo, hi = max(0, i - 2 * WINDOW), min(len(ins), i + 2 * WINDOW)
        depth = np.nanmax(ins[lo:hi]) - ins[i]
        if depth >= min_depth and all(abs(i - j) >= sep for j in kept):
            kept.append(i)
    return kept


# ------------------------------------------------------------ 1. predictions
maps, ctcf = {}, {}
for tissue, condition in TISSUES.items():
    out = model.predict_interval(
        interval, requested_outputs=[OT.CONTACT_MAPS, OT.CHIP_TF], condition=condition
    )
    cm = out.contact_maps
    if not maps:
        print("CONTACT_MAPS:", cm.values.shape, f"at {cm.resolution} bp per bin")
    print(f"  {tissue}: track {cm.metadata.name[0]} ({cm.metadata.imputation_mode[0]}), "
          f"range {cm.values.min():.2f}-{cm.values.max():.2f}")
    maps[tissue] = cm.values[:, :, 0]
    ctcf[tissue], n_ctcf = ctcf_profile(out.chip_tf, cm.values.shape[0])
    print(f"    observed CTCF tracks averaged: {n_ctcf}")

RES = 2048
N = 512
colon = maps[PRIMARY]  # reference map of the primary tissue
print(f"Symmetric: max |M - M.T| = {np.abs(colon - colon.T).max():.3f}")
dist = np.abs(np.subtract.outer(np.arange(N), np.arange(N)))


def distance_normalise(mat):
    """Divide each diagonal by its mean, removing the decay with distance."""
    expected = np.array([mat[dist == d].mean() for d in range(N)])
    return mat / expected[dist]


print("Mean contact by genomic distance (the values still decay with distance):")
for t, m in maps.items():
    print(f"  {t}:", {f"{d * RES // 1000} kb": round(float(m[dist == d].mean()), 3) for d in [1, 10, 50, 200]})
norm = {t: distance_normalise(m) for t, m in maps.items()}
print("Spread of the distance-normalised map (std of log2, 20-400 kb), i.e. how much structure it has:")
band = (dist >= 10) & (dist <= 200)
for t in maps:
    print(f"  {t}: {np.log2(norm[t][band]).std():.3f}")

# --------------------------------------- 2. boundaries vs CTCF (primary)
ins = insulation(colon)
bnd = boundaries(ins)
coord = lambda i: interval.start + i * RES
c = ctcf[PRIMARY]
ctcf_thr = np.quantile(c, 0.95)
rows = []
for i in bnd:
    near = c[max(0, i - 2) : i + 3]
    rows.append({"bin": i, "position": coord(i), "insulation": ins[i], "max_ctcf_within_2_bins": near.max(),
                 "ctcf_peak_nearby": bool(near.max() >= ctcf_thr)})
bdf = pd.DataFrame(rows)
print(f"\n{len(bdf)} boundaries called (insulation minima); CTCF peak = top 5% of bins:")
print(bdf.round(3).to_string(index=False))
valid = np.flatnonzero(np.isfinite(ins))
bg = np.mean([c[max(0, i - 2) : i + 3].max() >= ctcf_thr for i in valid])
print(f"Boundaries with a CTCF peak within 2 bins: {bdf.ctcf_peak_nearby.mean():.0%} "
      f"(background over all bins: {bg:.0%})")
r = np.corrcoef(ins[valid], c[valid])[0, 1]
print(f"Correlation between insulation score and CTCF signal across bins: r = {r:.3f}")

# ------------------------------------------------------ 3. across tissues
names = list(maps)
print("\nCorrelation of distance-normalised log2 maps between tissues (20-400 kb):")
corr = pd.DataFrame(
    [[np.corrcoef(np.log2(norm[x][band]), np.log2(norm[y][band]))[0, 1] for y in names] for x in names],
    index=names, columns=names)
print(corr.round(3).to_string())

# --------------------------------- 4. scramble the CTCF site at a boundary
with_ctcf = bdf[bdf.ctcf_peak_nearby]
target = int((with_ctcf if len(with_ctcf) else bdf).sort_values("insulation").bin.iloc[0])
peak_bin = max(range(max(0, target - 2), target + 3), key=lambda i: c[i])
url = ("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;"
       f"chrom={interval.chromosome};start={interval.start};end={interval.end}")
with urllib.request.urlopen(url) as resp:
    ref_seq = json.load(resp)["dna"].upper()
rng = np.random.default_rng(0)
a, b = peak_bin * RES, (peak_bin + 1) * RES
ctrl_bin = int(valid[np.argmin(np.abs(c[valid] - np.median(c)) + 1e3 * (np.abs(valid - peak_bin) < 60))])


def scrambled(seq, bin_index):
    s, e = bin_index * RES, (bin_index + 1) * RES
    return seq[:s] + "".join(rng.permutation(list(seq[s:e]))) + seq[e:]


edits = {}
for name, bin_index in [("CTCF bin scrambled", peak_bin), ("control bin scrambled", ctrl_bin)]:
    out = model.predict_sequence(
        scrambled(ref_seq, bin_index), requested_outputs=[OT.CONTACT_MAPS, OT.CHIP_TF],
        condition=TISSUES[PRIMARY], interval=interval)
    edits[name] = (out.contact_maps.values[:, :, 0], ctcf_profile(out.chip_tf, N)[0], bin_index)

print(f"\nScrambling one 2 kb bin ({PRIMARY}). CTCF bin = {peak_bin} (chr22:{coord(peak_bin)}), "
      f"control bin = {ctrl_bin} (chr22:{coord(ctrl_bin)}):")
for name, (m, cc, bin_index) in edits.items():
    d = np.log2(m) - np.log2(colon)
    ins_new = insulation(m)
    cross = d[max(0, bin_index - 40) : bin_index, bin_index + 1 : bin_index + 41].mean()
    print(f"  {name}: CTCF signal at bin {c[bin_index]:.2f} -> {cc[bin_index]:.2f}; "
          f"insulation at bin {ins[bin_index]:.3f} -> {ins_new[bin_index]:.3f}; "
          f"mean log2 change of contacts spanning the bin (+/-80 kb) = {cross:+.3f}; "
          f"max |log2 change| anywhere = {np.abs(d).max():.3f}")

# ------------------------------------------------------------------ figure
ext = [interval.start / 1e6, interval.end / 1e6, interval.end / 1e6, interval.start / 1e6]
xs = (interval.start + (np.arange(N) + 0.5) * RES) / 1e6
fig = plt.figure(figsize=(15, 12))
gs = fig.add_gridspec(4, 3, height_ratios=[4, 0.8, 0.8, 4], hspace=0.35, wspace=0.25)
for k, t in enumerate(names):
    ax = fig.add_subplot(gs[0, k])
    im = ax.imshow(np.log2(norm[t]), cmap="RdBu_r", vmin=-1, vmax=1, extent=ext)
    ax.set_title(f"{t}\nlog2(contact / mean at that distance)", fontsize=10)
    if t == PRIMARY:
        for i in bnd:
            ax.plot(coord(i) / 1e6, coord(i) / 1e6, "k^", ms=4)
    fig.colorbar(im, ax=ax, shrink=0.7)
ax = fig.add_subplot(gs[1, :])
ax.plot(xs, ins, color="black", lw=1)
ax.set_ylabel("insulation", rotation=0, ha="right", va="center")
ax.set_title(f"{PRIMARY}: insulation score and predicted CTCF ChIP-seq (orange = called boundaries)", fontsize=10)
ax2 = fig.add_subplot(gs[2, :])
ax2.fill_between(xs, c, color="tab:green", lw=0)
ax2.set_ylabel("CTCF", rotation=0, ha="right", va="center")
for axis in (ax, ax2):
    axis.set_xlim(ext[0], ext[1])
    for i in bnd:
        axis.axvline(coord(i) / 1e6, color="tab:orange", lw=0.8)
panels = [(f"{names[0]} - {PRIMARY}\n(distance-normalised)", np.log2(norm[names[0]]) - np.log2(norm[PRIMARY]), None, 1.0)]
panels += [(f"{n} - reference", np.log2(m) - np.log2(colon), bi, 0.5) for n, (m, _, bi) in edits.items()]
for k, (title, d, bi, lim) in enumerate(panels):
    ax = fig.add_subplot(gs[3, k])
    im2 = ax.imshow(d, cmap="PuOr_r", vmin=-lim, vmax=lim, extent=ext)
    ax.set_title(f"{title}\nlog2 ratio", fontsize=10)
    if bi is not None:
        ax.plot(coord(bi) / 1e6, coord(bi) / 1e6, "k^", ms=5)
    ax.set_xlabel("chr22 position (Mb)")
    fig.colorbar(im2, ax=ax, shrink=0.7)
out_png = OUT_DIR / "v2_contact_maps.png"
fig.savefig(out_png, dpi=140, bbox_inches="tight")
print(f"\nSaved {out_png}")
