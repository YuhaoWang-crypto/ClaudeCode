"""Digitize Dauphin-Ducharme et al. 2017 (Langmuir 33:4407), SI Figure S-7.

Source: 20-nt linear DNA, 5'-C6-thiol anchor, 3'-MB, 6-mercaptohexanol SAM, 2 mm Au,
SWV 0 -> -0.45 V vs Ag|AgCl, 1 mV step, f = 20/100/500 Hz x E_sw = 25/50/75 mV.
Solid lines = experiment (red forward, black backward, blue net); open circles = the
authors' COMSOL fit, which we deliberately DROP (morphological erosion removes the
thin circle outlines; the 3-4 px solid line survives).

Outputs (data/dauphin2017/):
  swv_digitized.csv   condition, f_Hz, Esw_V, curve, E_V, I_uA, sigma_read_uA
  calibration.json    per-panel tick pixel positions, linear map, tick residuals
  digitize_overlay.png  digitized points over the source image (visual audit)

The resulting numbers are FIGURE READINGS, not raw instrument data: one electrode
per condition is shown, samples along a curve are strongly correlated, and the
reading error below covers only pixel/calibration error, not experimental replicate
variance.
"""
import json
import os

import numpy as np
from PIL import Image
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
DDIR = os.path.join(HERE, "data", "dauphin2017")
SRC = os.path.join(DDIR, "figS7_source.jpeg")

FREQS = [20.0, 100.0, 500.0]
AMPS = [0.025, 0.050, 0.075]
Y_TICK_VALUES = {0: [0.1, 0.05, 0.0, -0.05], 1: [0.2, 0.1, 0.0, -0.1], 2: [0.4, 0.2, 0.0, -0.2]}
X_TICK_VALUES = [-0.2, -0.4]
E_GRID = np.round(np.arange(-0.060, -0.4501, -0.005), 4)   # 5 mV resampling grid


def _runs(idx):
    """Group sorted integer indices into runs; return run centers."""
    if len(idx) == 0:
        return []
    groups, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] <= 1:
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    groups.append(cur)
    return [float(np.mean(g)) for g in groups]


def find_axes(dark):
    """Return [(x_axis_row, [y_axis_col for each of 3 panels in that row])]."""
    rs = dark.sum(1)
    xrows = [int(round(r)) for r in _runs(np.where(rs > 1000)[0])]
    out = []
    for ay in xrows:
        line = np.r_[0, dark[ay, :].astype(int), 0]
        starts = np.where(np.diff(line) == 1)[0]
        ends = np.where(np.diff(line) == -1)[0]
        cols = []
        for s, e in zip(starts, ends):
            if e - s > 50:
                # y axis = column near the run start with the longest dark run above the x axis
                cand = range(s, s + 16)
                score = [dark[ay - 220:ay, c].sum() for c in cand]
                cols.append(int(list(cand)[int(np.argmax(score))]))
        out.append((ay, cols))
    return out


def panel_calibration(dark, ax_x, ax_y, row):
    # plot top = first y tick region; search ticks on the left of the y axis
    top = ax_y - 260
    col = dark[top:ax_y + 3, ax_x - 6:ax_x - 1].sum(1)
    yt = [top + r for r in _runs(np.where(col >= 4)[0])]
    # remove spurious hits (axis-label glyphs) by requiring near-uniform spacing
    yt = [y for y in yt if y >= top]
    vals = Y_TICK_VALUES[row]
    if len(yt) != len(vals):
        # keep the len(vals) ticks with the most regular spacing ending at the axis
        yt = sorted(yt)[-len(vals):]
    xr = dark[ax_y + 2:ax_y + 7, ax_x + 5:ax_x + 345].sum(0)
    xt = [ax_x + 5 + c for c in _runs(np.where(xr >= 4)[0])]
    assert len(xt) == 2, (row, ax_x, xt)
    py = np.polyfit(yt, vals, 1)
    yres = np.array(vals) - np.polyval(py, yt)
    px = np.polyfit(xt, X_TICK_VALUES, 1)
    return dict(y_ticks_px=[float(v) for v in yt], y_tick_values=vals, y_map=py.tolist(),
                y_tick_resid_uA=yres.tolist(), x_ticks_px=[float(v) for v in xt],
                x_map=px.tolist(), plot_x=[ax_x + 3, ax_x + 345], plot_y=[top, ax_y - 2])


def color_masks(a):
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    red = (r > 170) & (g < 110) & (b < 110)
    blue = (b > 150) & (r < 120) & (g > 90) & (g < 200)
    black = (r < 90) & (g < 90) & (b < 90)
    return {"forward": red, "backward": black, "net": blue}


def trace(mask, cal):
    x0, x1 = cal["plot_x"]
    y0, y1 = cal["plot_y"]
    sub = mask[y0:y1, x0:x1]
    er = ndimage.binary_erosion(sub, structure=np.ones((2, 2)))
    lab, n = ndimage.label(er, structure=np.ones((3, 3)))
    if n == 0:
        return np.array([]), np.array([]), np.array([])
    # keep components spanning >= 15 columns (line pieces); circles span ~8 px
    keep = np.zeros_like(er)
    for k in range(1, n + 1):
        cols = np.where((lab == k).any(0))[0]
        if cols.max() - cols.min() >= 15:
            keep |= lab == k
    xs, ys, sp = [], [], []
    for c in range(keep.shape[1]):
        rows = np.where(keep[:, c])[0]
        if len(rows):
            xs.append(x0 + c)
            ys.append(y0 + np.median(rows))
            sp.append(rows.max() - rows.min() + 1)
    return np.array(xs, float), np.array(ys, float), np.array(sp, float)


def main(make_plot=True):
    a = np.asarray(Image.open(SRC).convert("RGB")).astype(int)
    dark = a.sum(2) < 250
    axes = find_axes(dark)
    assert len(axes) == 3 and all(len(c) == 3 for _, c in axes), axes
    masks = color_masks(a)
    rows_out, calib = [], {}
    overlay = []
    for i, (ay, ycols) in enumerate(axes):
        for j, ax in enumerate(ycols):
            f, A = FREQS[i], AMPS[j]
            cond = f"f{int(f)}_A{int(A*1000)}"
            cal = panel_calibration(dark, ax, ay, i)
            calib[cond] = cal
            uA_per_px = abs(cal["y_map"][0])
            cal_err = float(np.max(np.abs(cal["y_tick_resid_uA"])))
            got = {}
            for curve, m in masks.items():
                xs, ys, sp = trace(m, cal)
                E = np.polyval(cal["x_map"], xs)
                I = np.polyval(cal["y_map"], ys)
                overlay.append((xs, ys, curve))
                order = np.argsort(E)
                E, I, sp = E[order], I[order], sp[order]
                pts = {}
                for Eg in E_GRID:
                    d = np.abs(E - Eg)
                    if d.size == 0 or d.min() > 0.006:   # no trace within 6 mV
                        continue
                    w = d < 0.003
                    if not w.any():
                        w = d == d.min()
                    # reading sigma: half the local line width (px) + calibration residual
                    s_ = float(np.sqrt((0.5 * np.median(sp[w]) * uA_per_px) ** 2 + cal_err ** 2))
                    pts[float(Eg)] = (float(np.median(I[w])), s_)
                got[curve] = pts
            net = got["net"]
            for curve in ("forward", "backward", "net"):
                pts = got[curve]
                vis = sorted(pts)
                for Eg in E_GRID:
                    Eg = float(Eg)
                    if Eg in pts:
                        rows_out.append((cond, f, A, curve, Eg, pts[Eg][0], pts[Eg][1], "visible"))
                    elif vis and min(vis) < Eg < max(vis):
                        # interior gap: linear interpolation between visible neighbours
                        lo = max(e for e in vis if e < Eg)
                        hi = min(e for e in vis if e > Eg)
                        t = (Eg - lo) / (hi - lo)
                        v = (1 - t) * pts[lo][0] + t * pts[hi][0]
                        rows_out.append((cond, f, A, curve, Eg, v, 2 * max(pts[lo][1], pts[hi][1]), "interp"))
                    elif Eg in net:
                        # occluded by the net (blue) trace drawn on top: all three coincide here
                        rows_out.append((cond, f, A, curve, Eg, net[Eg][0], 2 * net[Eg][1], "occluded_eq_net"))
    with open(os.path.join(DDIR, "calibration.json"), "w") as fh:
        json.dump(calib, fh, indent=1)
    with open(os.path.join(DDIR, "swv_digitized.csv"), "w") as fh:
        fh.write("condition,f_Hz,Esw_V,curve,E_V,I_uA,sigma_read_uA,flag\n")
        for r in rows_out:
            fh.write("%s,%g,%g,%s,%.4f,%.6f,%.6f,%s\n" % r)
    if make_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axp = plt.subplots(figsize=(13, 11))
        axp.imshow(a.astype(np.uint8))
        cmap = {"forward": "lime", "backward": "magenta", "net": "yellow"}
        for xs, ys, c in overlay:
            axp.plot(xs, ys, ".", ms=0.8, color=cmap[c])
        axp.set_axis_off()
        fig.tight_layout()
        fig.savefig(os.path.join(DDIR, "digitize_overlay.png"), dpi=110)
    return rows_out, calib


if __name__ == "__main__":
    rows, calib = main()
    print(len(rows), "grid points")
    for k, c in calib.items():
        print(k, "y ticks", [round(v, 1) for v in c["y_ticks_px"]],
              "x ticks", c["x_ticks_px"], "max tick resid uA %.4f" % max(map(abs, c["y_tick_resid_uA"])))
