"""
LF5 — Nuclear fractal dimension & optical texture entropy as a cheap readout
of chromatin state ("replace scATAC-seq with a microscope photo").

Data: BBBC021 (Caie et al. 2010 Mol Cancer Ther; Ljosa et al. 2012 Nat
Methods), MCF-7 breast-cancer cells, DAPI channel, plate Week1_22123.  Only
the needed TIFFs are pulled out of the 840 MB plate zip by HTTP range reads.
Labels are the compound treatments — independent of any image feature:
    DMSO (control) | taxol (microtubule stabiliser -> mitotic arrest,
    condensed chromosomes) | nocodazole (destabiliser -> mitotic arrest) |
    AZ-A (Aurora-kinase inhibitor -> failed cytokinesis, polyploid nuclei) |
    cytochalasin B (actin disruptor -> binucleation).

Per nucleus:
    * intensity-normalised Shannon entropy of the DAPI histogram (I / mean_I)
    * GLCM (co-occurrence) entropy of the normalised texture
    * fractal dimension from the grey-level structure function (variogram,
      FD = 3 − H) inside the (2-px eroded) nuclear mask.  Differential box
      counting (Sarkar & Chaudhuri 1994) was tried first and FAILED the
      synthetic check (it compresses FD 2.2–2.8 into 1.9–2.1), so it is only
      reported as a cautionary benchmark
    * area and integrated DAPI (DNA content) as conventional references
Estimator check: both estimators on synthetic fractional-Brownian surfaces of known
FD = 3 − H at the same patch size as the nuclei.
"""
import numpy as np
from . import _data

PLATE = "Week1/Week1_22123"


# ---------------------------------------------------------------- estimators
def fbm_surface(n, H, rng):
    kx = np.fft.fftfreq(n)[:, None]; ky = np.fft.fftfreq(n)[None, :]
    k = np.hypot(kx, ky); k[0, 0] = 1
    amp = k ** (-(H + 1))
    amp[0, 0] = 0
    phase = np.exp(2j * np.pi * rng.random((n, n)))
    z = np.real(np.fft.ifft2(amp * phase))
    return (z - z.min()) / (z.max() - z.min())


def dbc_fd(patch, scales=(2, 3, 4, 6, 8)):
    """Differential box-counting fractal dimension of a grey-level patch."""
    M = min(patch.shape)
    p = patch[:M, :M].astype(float)
    p = (p - p.min()) / (p.max() - p.min() + 1e-12) * 255
    G = 256.0
    logs, logN = [], []
    for s in scales:
        if M // s < 2:
            continue
        h = s * G / M
        m = (M // s) * s
        blocks = p[:m, :m].reshape(m // s, s, m // s, s).swapaxes(1, 2)
        mx = blocks.max((2, 3)); mn = blocks.min((2, 3))
        nr = np.ceil((mx + 1e-9) / h) - np.ceil((mn + 1e-9) / h) + 1
        logs.append(np.log(M / s)); logN.append(np.log(nr.sum()))
    return np.polyfit(logs, logN, 1)[0]


def variogram_fd(img, mask=None, lags=(1, 2, 3, 4)):
    """Fractal dimension from the structure function inside a mask:
    γ(r) = <(I(x+r) − I(x))²> ∝ r^{2H},  FD = 3 − H  (grey-level surface)."""
    img = img.astype(float)
    if mask is None:
        mask = np.ones(img.shape, bool)
    g = []
    for r in lags:
        d = []
        for a, b, ma, mb in ((img[:, r:], img[:, :-r], mask[:, r:], mask[:, :-r]),
                             (img[r:, :], img[:-r, :], mask[r:, :], mask[:-r, :])):
            m = ma & mb
            d.append((a[m] - b[m]) ** 2)
        d = np.concatenate(d)
        if len(d) < 20:
            return np.nan
        g.append(d.mean())
    H = np.polyfit(np.log(lags), np.log(g), 1)[0] / 2
    return 3 - np.clip(H, 0, 1)


def hist_entropy(v, bins=32, rng=(0, 3)):
    h, _ = np.histogram(v, bins=bins, range=rng)
    p = h[h > 0] / h.sum()
    return -(p * np.log(p)).sum()


def glcm_entropy(patch, mask, levels=16):
    from skimage.feature import graycomatrix
    q = np.clip(patch / (patch[mask].mean() + 1e-9) / 3 * (levels - 1), 0, levels - 1)
    q = q.astype(np.uint8)
    q[~mask] = 0
    g = graycomatrix(q, [1, 2], [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4], levels=levels)
    g = g[1:, 1:].astype(float)        # drop pairs involving masked-out pixels (0)
    p = g.sum((2, 3)); p = p / p.sum()
    p = p[p > 0]
    return -(p * np.log(p)).sum()


def validate_fd(sizes=(24, 32, 64, 256), Hs=(0.2, 0.5, 0.8), reps=40, seed=0):
    """Both estimators on fBm patches of known FD; also the variogram with
    additive pixel noise (as in a real camera image, SNR ≈ 10)."""
    rng = np.random.default_rng(seed)
    out = {}
    for n in sizes:
        for H in Hs:
            e_dbc, e_var, e_noisy = [], [], []
            for _ in range(reps):
                big = fbm_surface(256, H, rng)
                pt = big[:n, :n]
                e_dbc.append(dbc_fd(pt)); e_var.append(variogram_fd(pt))
                e_noisy.append(variogram_fd(pt + rng.normal(0, 0.1 * pt.std(), pt.shape)))
            out[(n, H)] = dict(dbc=(np.mean(e_dbc), np.std(e_dbc)),
                               var=(np.mean(e_var), np.std(e_var)),
                               noisy=(np.mean(e_noisy), np.std(e_noisy)))
    return out


# ---------------------------------------------------------------- segmentation
def segment(img):
    from scipy import ndimage as ndi
    from skimage import filters, measure, morphology, segmentation
    from skimage.feature import peak_local_max
    sm = filters.gaussian(img, 1.5, preserve_range=True)
    bg = np.percentile(sm, 20)
    m = sm > filters.threshold_otsu(sm)
    m = morphology.remove_small_objects(m, 150)
    m = ndi.binary_fill_holes(m)
    dist = ndi.distance_transform_edt(m)
    pk = peak_local_max(filters.gaussian(dist, 2), min_distance=9, labels=measure.label(m))
    mk = np.zeros(m.shape, int); mk[tuple(pk.T)] = np.arange(1, len(pk) + 1)
    lab = segmentation.watershed(-dist, mk, mask=m)
    lab = segmentation.clear_border(lab)
    return lab, bg


def nucleus_features(img, lab, bg):
    from scipy import ndimage as ndi
    from skimage import measure
    rows = []
    img = img.astype(float) - bg
    for p in measure.regionprops(lab, intensity_image=img):
        if p.area < 150 or p.area > 6000 or p.solidity < 0.85:
            continue
        r0, c0, r1, c1 = p.bbox
        crop = img[r0:r1, c0:c1]; msk = p.image
        v = crop[msk]
        mu = v.mean()
        if mu <= 0:
            continue
        # erode by 2 px so the nuclear rim (a step edge) does not dominate γ(r)
        inner = ndi.binary_erosion(msk, iterations=2)
        fd = variogram_fd(crop, inner) if inner.sum() > 60 else np.nan
        rows.append(dict(area=p.area, dna=v.sum(), mean=mu,
                         H_hist=hist_entropy(v / mu), H_glcm=glcm_entropy(crop, msk),
                         fd=fd, cv=v.std() / mu))
    return rows


def collect(max_per_group=24, seed=0):
    import pandas as pd
    import tifffile
    m = _data.bbbc021_metadata()
    sub = m[m.Image_PathName_DAPI == PLATE].copy()
    groups = {"DMSO": sub.compound == "DMSO",
              "taxol": sub.compound == "taxol",
              "nocodazole": sub.compound == "nocodazole",
              "AZ-A (Aurora i.)": (sub.compound == "AZ-A") & (sub.concentration >= 1.0),
              "cytochalasin B": sub.compound == "cytochalasin B"}
    rng = np.random.default_rng(seed)
    rows = []
    for gname, sel in groups.items():
        files = sub[sel].Image_FileName_DAPI.tolist()
        if len(files) > max_per_group:
            files = list(rng.choice(files, max_per_group, replace=False))
        paths = _data.bbbc021_dapi(PLATE, files)
        for f, pth in zip(files, paths):
            img = tifffile.imread(pth)
            lab, bg = segment(img)
            for r in nucleus_features(img, lab, bg):
                r.update(group=gname, image=f); rows.append(r)
    return pd.DataFrame(rows)


def _auc(x, y):
    from scipy.stats import mannwhitneyu
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    u = mannwhitneyu(x, y).statistic / (len(x) * len(y))
    return max(u, 1 - u)


FEATS = [("H_hist", "intensity entropy"), ("H_glcm", "GLCM texture entropy"),
         ("fd", "fractal dim. (variogram)"), ("area", "nuclear area"), ("dna", "integrated DAPI")]


def run(make_fig=True):
    print("=" * 72); print("LF5  Nuclear texture entropy & fractal dimension (BBBC021 MCF-7)"); print("=" * 72)
    V = validate_fd()
    print("[check] FD estimators on fBm surfaces of known FD = 3 − H (mean ± sd, 40 reps):")
    for key, name in (("dbc", "box-counting (DBC)"), ("var", "variogram"),
                      ("noisy", "variogram + 10% pixel noise")):
        print(f"    {name}:")
        for n in (24, 32, 64, 256):
            print(f"      patch {n:3d}px: " + "  ".join(
                f"true {3 - H:.1f} -> {V[(n, H)][key][0]:.2f}±{V[(n, H)][key][1]:.2f}"
                for H in (0.2, 0.5, 0.8)))
    df = collect()
    ctrl = df[df.group == "DMSO"]
    print(f"\n[data] {len(df)} nuclei from {df.image.nunique()} images; per group:",
          dict(df.group.value_counts()))
    print("    median feature per group; AUC vs DMSO (0.5 = no signal, 1 = perfect):")
    hdr = "    " + f"{'group':18s}" + "".join(f"{lab:>27s}" for _, lab in FEATS)
    print(hdr)
    table = {}
    for gname in df.group.unique():
        g = df[df.group == gname]
        cells = []
        for f, _ in FEATS:
            med = np.nanmedian(g[f])
            auc = _auc(g[f].values, ctrl[f].values) if gname != "DMSO" else np.nan
            cells.append((med, auc))
        table[gname] = cells
        print("    " + f"{gname:18s}" + "".join(
            f"{med:>16.3g} (AUC {auc:.2f})" if np.isfinite(auc) else f"{med:>16.3g}           "
            for med, auc in cells))
    # image-level: does the median texture entropy of an image separate taxol from DMSO?
    img = df.groupby(["image", "group"]).median(numeric_only=True).reset_index()
    img_auc = {f: _auc(img[img.group == "taxol"][f].values, img[img.group == "DMSO"][f].values)
               for f, _ in FEATS}
    print("    image-level AUC taxol vs DMSO: " +
          ", ".join(f"{lab} {img_auc[f]:.2f}" for f, lab in FEATS))
    # does texture add information beyond size/DNA? partial out area & DNA
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    two = df[df.group.isin(["DMSO", "taxol"])].dropna()
    y = (two.group == "taxol").values
    base = two[["area", "dna"]].apply(np.log).values
    full = np.column_stack([base, two[["H_hist", "H_glcm", "fd"]].values])
    cv_base = cross_val_score(LogisticRegression(max_iter=2000), base, y, cv=5, scoring="roc_auc").mean()
    cv_full = cross_val_score(LogisticRegression(max_iter=2000), full, y, cv=5, scoring="roc_auc").mean()
    print(f"    nucleus-level taxol-vs-DMSO classifier (5-fold CV AUC): size+DNA {cv_base:.2f} "
          f"-> + entropy/FD {cv_full:.2f}")

    if make_fig:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import tifffile
        fig, ax = plt.subplots(1, 5, figsize=(22, 4.4))
        for k, gname in enumerate(["DMSO", "taxol"]):
            f = df[df.group == gname].image.iloc[0]
            im = tifffile.imread(_data.bbbc021_dapi(PLATE, [f])[0]).astype(float)
            ax[k].imshow(im[300:556, 400:656], cmap="gray",
                         vmin=np.percentile(im, 1), vmax=np.percentile(im, 99.7))
            ax[k].set_title(f"{gname} (DAPI, 256×256 crop)"); ax[k].axis("off")
        order = ["DMSO", "taxol", "nocodazole", "AZ-A (Aurora i.)", "cytochalasin B"]
        for k, (f, lab) in enumerate(FEATS[:3]):
            data = [df[df.group == g][f].dropna().values for g in order]
            ax[2 + k].boxplot(data, showfliers=False)
            ax[2 + k].set_xticks(range(1, 6))
            ax[2 + k].set_xticklabels([o.split(" ")[0] for o in order], rotation=30)
            ax[2 + k].set_title(lab)
        fig.tight_layout(); fig.savefig(_data.fig("lf5_nuclear_texture.png"), dpi=130)
        plt.close(fig)
    return dict(validation=V, table=table, img_auc=img_auc, cv=(cv_base, cv_full), df=df)


if __name__ == "__main__":
    run()
