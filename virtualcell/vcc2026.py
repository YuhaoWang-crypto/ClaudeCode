"""Predict the official Virtual Cell Challenge 2026 validation panel.

The 2026 task, as the downloaded bundle actually defines it: three anonymised
contexts A, B and C, each supplied as **18,400 non-targeting control cells and
nothing else**, on an 18,533-gene axis; predict 400 cells for each of 300 gene
knockdowns in every context.  The ground truth -- 138,400 cells per context --
is withheld and scored server-side, so nothing here can be scored locally; what
can be checked locally is that the model has anything to say about this panel at
all, and that the submission is well-formed.

Two facts decide the design, both established by
:mod:`virtualcell.gwps` and reproducible from the numbers this module prints:

* **The four-line atlas covers 0 of the 300 knockdowns.**  Replogle's and
  Nadig's essential-gene screens target what a cell needs to survive; the
  official panel targets regulators that it does not.  A model built on that
  atlas can only emit the average response for all 300, which is exactly the
  baseline the leaderboard normalises against.
* **The genome-wide K562 arm covers 272 of them**, so it becomes the source.
  It measures 8,246 genes, 7,681 of which are on the official axis; the
  remaining 10,852 official genes are unmeasured in any perturbation data
  available here and are predicted unchanged.

Usage::

    python -m virtualcell.vcc2026 --describe     # what is in the bundle
    python -m virtualcell.vcc2026 --build        # write prediction.h5ad
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

from dataclasses import replace

from . import metrics as M
from .data import CellLine
from .gwps import load_gwps, project
from .model import ContextTransferModel, Hyper, SourceBank
from .vcc_submit import CONTEXTS, build

OFFICIAL = Path("/home/user/vcc_official")


# --------------------------------------------------------------------------
# the official bundle
# --------------------------------------------------------------------------

def gene_axis() -> np.ndarray:
    return np.loadtxt(OFFICIAL / "gene_names.csv", dtype=str, skiprows=1)


def panel() -> np.ndarray:
    return np.loadtxt(OFFICIAL / "pert_counts.csv", dtype=str, skiprows=1)


def manifest() -> dict:
    return json.loads((OFFICIAL / "manifest.json").read_text())


def control_cells(context: str) -> sparse.csr_matrix:
    """Load one context's control cells as raw-count CSR."""
    with h5py.File(OFFICIAL / f"context_{context}.h5ad") as f:
        X = f["X"]
        shape = tuple(int(v) for v in X.attrs["shape"])
        return sparse.csr_matrix(
            (X["data"][:].astype(np.int32), X["indices"][:], X["indptr"][:]),
            shape=shape)


def control_mean(context: str) -> np.ndarray:
    """Pseudobulk log1p CP10K of one context's controls, without holding the cells.

    Streaming the CSR buffers keeps peak memory at a few hundred megabytes
    instead of the ~900 MB the full matrix costs, which matters because three
    contexts and a 9,866-knockdown source bank have to coexist.
    """
    with h5py.File(OFFICIAL / f"context_{context}.h5ad") as f:
        X = f["X"]
        n, g = (int(v) for v in X.attrs["shape"])
        ip = X["indptr"][:]
        tot = np.zeros(g)
        for s in range(0, n, 2000):
            e = min(s + 2000, n)
            a, b = ip[s], ip[e]
            tot += np.bincount(X["indices"][a:b],
                               weights=X["data"][a:b].astype(np.float64),
                               minlength=g)
    return np.log1p(tot / tot.sum() * 1e4)


# --------------------------------------------------------------------------
# the model, on the official axis
# --------------------------------------------------------------------------

def source_on_official(axis: np.ndarray) -> tuple[CellLine, np.ndarray]:
    """Genome-wide K562 screen restricted to genes shared with the official axis.

    Returns the source line and the column indices it occupies in ``axis``.
    The source is keyed by ENSG internally but the submission axis is keyed by
    symbol, so the join is by symbol with the first ENSG per symbol winning --
    the same rule :func:`virtualcell.vcc_submit.dedupe_genes` applies to the
    output.
    """
    src = load_gwps()
    first: dict[str, int] = {}
    for i, s in enumerate(src.symbols.astype(str)):
        first.setdefault(s, i)

    cols, rows = [], []
    for j, sym in enumerate(axis):
        i = first.get(sym)
        if i is not None:
            cols.append(j)
            rows.append(i)
    cols = np.array(cols)
    return project(src, src.genes[np.array(rows)]), cols


def predict_panel(hp: Hyper, targets: np.ndarray | None = None,
                  verbose: bool = True) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Predicted effects for every context, on the full official gene axis."""
    axis, perts = gene_axis(), panel()
    if targets is not None:
        perts = targets
    src, cols = source_on_official(axis)
    if verbose:
        covered = sum(p in set(src.names) for p in perts)
        print(f"  source: {src.names.size:,} knockdowns over {cols.size:,} of "
              f"{axis.size:,} official genes; {covered}/{perts.size} of the panel "
              f"measured directly")

    # With one source line the context weights are fixed at 1, so the bank -- the
    # consensus, the similarity graph, the response basis and the ~9,866-row SVD
    # behind it -- is the same for all three contexts and is built once.  Only
    # the modulation and the on-target term differ per context.
    bank = SourceBank.build([src], np.ones(1))

    out, model = {}, None
    for c in CONTEXTS:
        mu_full = control_mean(c)
        model = (model.retarget([src], mu_full[cols]) if model is not None else
                 ContextTransferModel(hp).fit([src], mu_full[cols],
                                              src.symbols, bank=bank))
        eff = np.zeros((perts.size, axis.size))
        eff[:, cols] = model.predict(perts)
        out[c] = (perts, eff)
        if verbose:
            nz = np.abs(eff).sum(1) > 0
            print(f"  context {c}: {nz.sum()}/{perts.size} knockdowns given a "
                  f"nonzero effect, median |effect| over predicted genes "
                  f"{np.median(np.abs(eff[:, cols])):.4f}, magnitude spread "
                  f"(norm CV) {M.norm_cv(eff[:, cols]):.3f}", flush=True)

    if verbose:
        # Predicts the discrimination score without needing the withheld truth:
        # a prediction whose magnitudes vary much less than the measured ones has
        # given away most of what L1 retrieval pays for.  See metrics.norm_cv.
        idx = {p: i for i, p in enumerate(src.names)}
        rows = [idx[p] for p in perts if p in idx]
        print(f"  reference: the same knockdowns measured in K562 have norm CV "
              f"{M.norm_cv(src.delta[rows]):.3f}")
    return out


HYPER_FILE = Path("results/gwps_single_source.json")
FRONTIER_FILE = Path("results/gwps_frontier.json")


def threshold_safe(hp: Hyper, path: Path = FRONTIER_FILE,
                   verbose: bool = True) -> Hyper:
    """Replace the tuned effect scale with one that meets every metric threshold.

    Two selection rules are in play and they disagree, for a reason worth being
    explicit about.  The tuner maximises the *unclipped* aggregate score, where a
    metric below baseline is a penalty that a large enough gain elsewhere can
    outweigh; in the single-source setting it duly buys discrimination and
    differential expression with mean absolute error.  The challenge does not
    score that way -- it "enforce[s] minimum thresholds on all metrics" -- so a
    submission tuned that way can be strong on two metrics and inadmissible on
    the third.

    :func:`virtualcell.gwps.frontier` measures where the crossing is.  This takes
    ``beta`` and ``gene_w`` from there and everything else from the tuner, which
    is the honest combination: the constraint comes from the metric the challenge
    enforces, the rest from cross-validation.
    """
    if not path.exists():
        if verbose:
            print(f"  {path} absent; keeping the tuned effect scale, which may "
                  f"fail the MAE threshold")
        return hp
    best = json.loads(path.read_text())["best"]
    usable = {k: v for k, v in best.items() if v.get("safe_beta") is not None}
    if not usable:
        if verbose:
            print("  no effect scale kept MAE at or below baseline on every "
                  "fold; keeping the tuned value and reporting the failure")
        return hp
    key = max(usable, key=lambda k: usable[k].get("score", 0.0))
    e = usable[key]
    if verbose:
        print(f"  effect scale from the frontier: beta {hp.beta:g} -> "
              f"{e['safe_beta']:g}, gene_w {hp.gene_w:g} -> {key} "
              f"(mean PDS {e['pds']:.3f}, ovl@100 {e['ovl']:.3f}, "
              f"score {e['score']:.3f})")
    return replace(hp, beta=float(e["safe_beta"]), gene_w=float(key))


def single_source_hyper(exclude: str | None = None,
                        path: Path = HYPER_FILE) -> Hyper:
    """Hyperparameters from the single-source benchmark in :mod:`virtualcell.gwps`.

    That benchmark tunes each fold on the *other* held-out lines, so the entry
    keyed by a line's name was selected without ever seeing it.  Passing
    ``exclude="Jurkat"`` therefore returns settings that are leak-free for
    scoring on Jurkat; passing nothing averages the folds, which is what the
    official submission uses since no official context is one of them.
    """
    if not path.exists():
        return Hyper(temp=np.inf)
    blob = json.loads(path.read_text())["hyper"]
    if exclude is not None:
        if exclude not in blob:
            raise KeyError(f"{path} has no fold for {exclude!r}; "
                           f"available: {', '.join(blob)}")
        folds = [blob[exclude]]
    else:
        folds = list(blob.values())

    agg = {}
    for k in folds[0]:
        vals = [f[k] for f in folds]
        d = getattr(Hyper(), k)
        if isinstance(d, bool):
            agg[k] = bool(round(float(np.mean([bool(v) for v in vals]))))
        elif isinstance(d, int):
            agg[k] = int(round(float(np.median([float(v) for v in vals]))))
        else:
            agg[k] = float(np.mean([float(v) for v in vals]))
    return Hyper(**agg)


# --------------------------------------------------------------------------

# Enough lineage markers to place a context, not to type it definitively.  The
# contexts are anonymised on purpose and the model never uses their identity --
# it conditions on the observed control profile -- but knowing what they are
# tells you how far the source screen is being asked to reach.
MARKERS = {
    "T lymphoid":            ["CD3E", "CD3D", "CD2", "LCK", "TRAC", "CD7", "IL7R"],
    "erythroid / CML":       ["HBG1", "HBG2", "GATA1", "HBZ", "KLF1", "GYPA"],
    "hepatocyte":            ["ALB", "APOA1", "APOB", "TTR", "SERPINA1", "FGB", "AFP"],
    "epithelial (RPE-like)": ["KRT8", "KRT18", "CDH1", "EPCAM", "FOLR1", "CLU"],
    "squamous epithelial":   ["KRT5", "KRT13", "KRT15", "TACSTD2", "CSTA", "FXYD3"],
    "mesenchymal":           ["COL1A1", "COL1A2", "FN1", "VIM", "THY1", "ACTA2"],
    "myeloid":               ["SPI1", "LYZ", "CD33", "ITGAM", "MPO", "CSF1R"],
    "neural":                ["SOX2", "NES", "TUBB3", "MAP2", "STMN2"],
}


def identify() -> None:
    """Place each anonymised context by lineage marker expression."""
    axis = gene_axis()
    idx = {g: i for i, g in enumerate(axis)}
    prof = {c: control_mean(c) for c in CONTEXTS}

    print(f"mean log1p CP10K over marker sets\n\n{'lineage':<24}"
          + "".join(f"{c:>8}" for c in CONTEXTS) + "   found")
    print("-" * 56)
    for label, genes in MARKERS.items():
        rows = [idx[g] for g in genes if g in idx]
        if not rows:
            continue
        print(f"{label:<24}"
              + "".join(f"{prof[c][rows].mean():>8.2f}" for c in CONTEXTS)
              + f"   {len(rows)}/{len(genes)}")

    P = np.stack([prof[c] for c in CONTEXTS])
    print("\nmost distinctive genes per context (excess over the other two)")
    for j, c in enumerate(CONTEXTS):
        sc = P[j] - np.delete(P, j, 0).max(0)
        top = np.argsort(sc)[::-1][:10]
        print(f"  {c}: " + ", ".join(f"{axis[i]}({sc[i]:+.1f})" for i in top))

    print("\npairwise correlation of the three control profiles")
    for a in range(3):
        for b in range(a + 1, 3):
            print(f"  {CONTEXTS[a]}-{CONTEXTS[b]}: "
                  f"{np.corrcoef(P[a], P[b])[0, 1]:.3f}")


def describe() -> None:
    """Print what the bundle contains and what the atlas can reach in it."""
    m = manifest()
    axis, perts = gene_axis(), panel()
    print(f"season {m['season']} / partition {m['partition']} / panel "
          f"{m['panel_id']}")
    print(f"  {len(m['contexts'])} contexts, {m['n_constructs']} knockdowns, "
          f"{m['n_genes']:,} genes, {m['cells_per_pert']} cells per knockdown")
    a = m["per_context"][m["contexts"][0]]
    print(f"  per context: {a['control_cells']:,} control cells given, "
          f"{a['ground_truth_cells']:,} ground-truth cells withheld")

    from .data import load_all
    lines, _, _ = load_all()
    print("\npanel coverage by source screen")
    for cl in lines:
        n = sum(p in set(cl.names) for p in perts)
        print(f"  {cl.name + ' (essential panel)':<32} {n:>3}/300")
    gw = load_gwps()
    tgt, meas = set(gw.names.astype(str)), set(gw.symbols.astype(str))
    print(f"  {'K562 genome-wide (Replogle)':<32} {sum(p in tgt for p in perts):>3}/300")

    # A gene not covered by the source is predicted unchanged, so the raw
    # fraction of the axis understates the damage badly: what is missing is the
    # long tail, and what is covered is most of what a cell actually transcribes.
    cov = np.array([g in meas for g in axis])
    print(f"\ngene-axis overlap: four-line atlas "
          f"{len(set(lines[0].symbols.astype(str)) & set(axis)):,}, "
          f"genome-wide screen {cov.sum():,} of {axis.size:,}")
    print(f"{'context':<9}{'share of expression':>21}{'of top-2000 expressed':>24}")
    for c in CONTEXTS:
        lin = np.expm1(control_mean(c))
        top = np.argsort(lin)[::-1][:2000]
        print(f"{c:<9}{lin[cov].sum() / lin.sum():>20.1%}{cov[top].mean():>23.1%}")

    miss = [p for p in perts if p not in tgt]
    routed = [m for m in miss if m in meas]
    print(f"\n{len(miss)} panel targets are not knocked down in the screen: "
          f"{len(routed)} are measured genes and route through the "
          f"functional-neighbour path ({', '.join(routed)}), the other "
          f"{len(miss) - len(routed)} fall back to the mean response")

    # The panel is not the essential-gene panel this project benchmarks on, and
    # not only in identity: it is measurably weaker, which is the harder regime.
    mag = np.linalg.norm(gw.delta, axis=1)
    idx = {p: i for i, p in enumerate(gw.names)}
    on_panel = np.array([mag[idx[p]] for p in perts if p in idx])
    essential = np.array([mag[idx[p]] for p in lines[0].names if p in idx])
    print(f"\nmedian L2 effect in K562: official panel {np.median(on_panel):.2f}, "
          f"essential panel {np.median(essential):.2f}, "
          f"all {gw.names.size:,} screened knockdowns {np.median(mag):.2f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--describe", action="store_true")
    ap.add_argument("--identify", action="store_true",
                    help="place the anonymised contexts by lineage markers")
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--out", type=Path,
                    default=Path("/home/user/vcc_submission/prediction.h5ad"))
    ap.add_argument("--n-cells", type=int, default=400)
    ap.add_argument("--hyper", type=Path, default=HYPER_FILE,
                    help="tuned hyperparameters from the single-source benchmark")
    ap.add_argument("--keep-tuned-scale", action="store_true",
                    help="do not override the tuned effect scale with the "
                         "threshold-safe one from the frontier")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.identify:
        identify()
        if not args.build:
            return
    if args.describe or not args.build:
        describe()
        if not args.build:
            return

    if args.hyper.exists():
        hp = single_source_hyper(path=args.hyper)
        n = len(json.loads(args.hyper.read_text())["hyper"])
        print(f"hyperparameters (mean over {n} held-out lines): {hp}")
    else:
        hp = Hyper(temp=np.inf)
        print(f"WARNING: {args.hyper} absent, using untuned defaults")
    if not args.keep_tuned_scale:
        hp = threshold_safe(hp)

    predictions = predict_panel(hp)
    build(predictions, {c: (lambda c=c: control_cells(c)) for c in CONTEXTS},
          gene_axis(), args.out, args.n_cells, args.seed)
    print("\nvalidate with:")
    print(f"  vcc prep {args.out} -g {OFFICIAL / 'gene_names.csv'} "
          f"--perts {OFFICIAL / 'pert_counts.csv'} -o prediction.vcc")


if __name__ == "__main__":
    main()
