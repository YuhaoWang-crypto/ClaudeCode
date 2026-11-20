"""Measure the one number the over-calling argument leaves unmeasured.

:mod:`virtualcell.wrapup_analysis` establishes from the published 2026 results
that ``sig_jaccard`` is capped at ``n_real / n_pred`` by arithmetic, that the
winning entry sits at ~80% of the cap its own ~500-gene call set imposes, and
that calling the right number instead would be worth more than the published
first-to-second margin.  All of that is certain.

What is *not* certain is the recall that survives the thresholding.  Keeping a
perturbation's top ``n_real`` predicted genes instead of its top 500 raises the
cap, but only helps if the genes kept are the right ones -- and that depends on
how well a real model ranks genes in a context it has never seen, which no
amount of arithmetic settles.  This module measures it.

**Setup.**  Source is Replogle's K562 CRISPRi screen, the reference screen 277
of the 300 official test targets appear in.  Target is Nadig's Jurkat, which
the 2026 wrap-up de-anonymised as one of the three validation lines, so the
measurement happens in a context the challenge actually scored.  The model sees
the source's measured effects and the target's control cells, and nothing else
-- the zero-shot setting the challenge defines.

**Three arms per prediction, and the middle one is the deployable one.**

* *uncalibrated* -- apply the predicted effect to every gene, as the field does;
* *calibrated* -- sparsify to counts rank-matched from the **source** context's
  own measured response sizes, which uses no target ground truth and is what
  one could actually submit;
* *oracle counts* -- sparsify to the target's own measured counts.  This is not
  submittable and is reported only as the ceiling the deployable arm is trying
  to reach, so that the cost of transferring the count distribution is visible
  rather than hidden.

**Two documented deviations from the official scoring.**  Nadig's Jurkat holds
262,956 cells over ~2,000 targets, so the official 400 cells per perturbation
is not available; each perturbation is scored on all the cells it has and the
prediction is given the *same* number, so the test's power is matched between
reference and prediction.  Fewer cells means less power, which calls *fewer*
genes -- so this understates over-calling rather than flattering the result.
And the control pool is thousands rather than 18,400: with 400 cells on the
other side the detectable effect scales as ``sqrt(1/n1 + 1/n2)``, which moves
about 3% between those two, far less than anything measured here.

    python -m virtualcell.de_recall
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

from . import vcc2026_metrics as M
from .data import CellLine
from .de_calibration import rank_matched_targets
from .model import ContextTransferModel, Hyper, SourceBank
from .vcc_submit import sample_cells

#: Hyperparameters from this project's single-source benchmark, Jurkat fold.
#: That fold tuned on RPE1 and HepG2 only, so using it on a Jurkat target is
#: leakage-free -- but it was tuned against the **2025** metrics, where L1
#: discrimination rewarded magnitude spread.  Re-tuning against the 2026 suite
#: is a separate job; this arm measures the model as it stands.
JURKAT_FOLD_HYPER = dict(
    temp=float("inf"), smooth=0.5, n_neighbors=5, shrink=0.5, gamma=0.5,
    mod_clip=2.0, beta=0.8, use_global=0.2, on_target=True, rank=80,
    rank_mix=0.25, unseen_k=25, renorm=0.0)

DATA = Path("/home/user/vcc_data/singlecell")
CONTROL = "non-targeting"


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------

def _load(path: Path) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    a = ad.read_h5ad(path)
    X = sparse.csr_matrix(a.X)
    label = np.asarray(a.obs["gene"]).astype(str)
    sym = np.asarray(a.var["gene_name"] if "gene_name" in a.var
                     else a.var_names).astype(str)
    return X, label, sym


def ctrl_cpm_mean(X: sparse.csr_matrix) -> np.ndarray:
    """Per-gene mean CPM over a group's cells, without densifying it.

    The expression gate needs this over ~8,000 control cells x 18,000 genes;
    densifying that to take a column mean costs over a gigabyte for a number
    that sparse arithmetic gives directly.
    """
    tot = np.asarray(X.sum(axis=1)).ravel()
    inv = np.zeros_like(tot)
    np.divide(1.0, tot, out=inv, where=tot > 0)
    scaled = X.multiply(inv[:, None]).tocsr()
    return np.asarray(scaled.sum(axis=0)).ravel() / X.shape[0] * M.TARGET_SUM_CELL


def log1p_cp10k(X: sparse.csr_matrix) -> np.ndarray:
    """Pseudobulk profile in the log1p CP10K space the model's effects live in."""
    total = np.asarray(X.sum(axis=0)).ravel()
    return np.log1p(total / max(total.sum(), 1.0) * 1e4)


def source_cellline(X: sparse.csr_matrix, label: np.ndarray, sym: np.ndarray,
                    name: str = "K562", n_splits: int = 10, min_cells: int = 10,
                    seed: int = 0) -> CellLine:
    """Pseudobulk a single-cell screen into the model's source representation.

    The control cells are split into ``n_splits`` disjoint pseudo-replicates
    rather than pooled, because :class:`CellLine` uses the spread *between*
    control replicates to estimate what a null difference looks like, and a
    single pooled control would leave that at zero.
    """
    rng = np.random.default_rng(seed)
    idx = np.flatnonzero(label == CONTROL)
    rng.shuffle(idx)
    groups = [g for g in np.array_split(idx, n_splits) if g.size]
    ctrl = np.stack([log1p_cp10k(X[np.sort(g)]) for g in groups])
    ctrl_ncells = np.array([g.size for g in groups])

    names, pert, ncells = [], [], []
    for n in sorted(set(label.tolist()) - {CONTROL}):
        rows = np.flatnonzero(label == n)
        if rows.size < min_cells:
            continue
        names.append(n)
        pert.append(log1p_cp10k(X[rows]))
        ncells.append(rows.size)
    return CellLine(name=name, ctrl=ctrl, pert=np.stack(pert),
                    names=np.array(names), ncells=np.array(ncells),
                    ctrl_ncells=ctrl_ncells, genes=sym, symbols=sym)


# ---------------------------------------------------------------------------
# the DE table, for one group against the reference control
# ---------------------------------------------------------------------------

def de_table(cells: sparse.csr_matrix, ctrl_dense: np.ndarray,
             ctrl_mean: np.ndarray, gate: np.ndarray
             ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(significant, log2fc, p_adjusted_proxy)`` over the gated genes.

    The comparison group is the reference's control cells on both sides, as the
    specification requires: the prediction is tested against the same cells the
    reference was.
    """
    g = M.cpm(np.asarray(cells[:, gate].todense(), dtype=np.float32))
    p = M.wilcoxon_p(g, ctrl_dense)
    sig = M.bh_reject(p)
    lfc = M.log2fc(g, ctrl_mean)
    # Returned as the third value and handed to `direction_reach` for both its
    # adjusted-p and raw-p ordering tiers.  Within one perturbation BH is
    # monotone in the raw p-value, so ordering by either gives the same
    # sequence, and the tiers only ever break ties.
    return sig, lfc, p


# ---------------------------------------------------------------------------
# the measurement
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", type=Path, default=DATA / "jurkat_de.h5ad")
    ap.add_argument("--source", type=Path, default=DATA / "k562_de.h5ad")
    ap.add_argument("--n-perts", type=int, default=100)
    ap.add_argument("--min-cells", type=int, default=80)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/de_recall.json"))
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    print("loading")
    Xt, lt, st = _load(args.target)
    Xs, ls, ss = _load(args.source)
    print(f"  target {Xt.shape[0]:,} cells x {Xt.shape[1]:,} genes, "
          f"{np.unique(lt).size} labels")
    print(f"  source {Xs.shape[0]:,} cells x {Xs.shape[1]:,} genes, "
          f"{np.unique(ls).size} labels")

    # common gene axis, by symbol
    shared_g = np.intersect1d(st, ss)
    ti = {s: i for i, s in enumerate(st)}
    si = {s: i for i, s in enumerate(ss)}
    tcol = np.array([ti[s] for s in shared_g])
    scol = np.array([si[s] for s in shared_g])
    Xt = Xt[:, tcol].tocsr()
    Xs = Xs[:, scol].tocsr()
    sym = shared_g
    print(f"  shared genes: {sym.size:,}")

    # source effects, pseudobulk, in log1p CP10K
    s_ctrl = log1p_cp10k(Xs[ls == CONTROL])
    src_effect: dict[str, np.ndarray] = {}
    for name in np.unique(ls):
        if name == CONTROL:
            continue
        rows = ls == name
        if rows.sum() < 10:
            continue
        src_effect[name] = log1p_cp10k(Xs[rows]) - s_ctrl

    # target: controls, and the perturbations that are also in the source
    t_ctrl = Xt[lt == CONTROL]
    counts_t = {n: int((lt == n).sum()) for n in np.unique(lt) if n != CONTROL}
    usable = [n for n, c in counts_t.items()
              if c >= args.min_cells and n in src_effect and n in set(sym)]
    usable = sorted(usable, key=lambda n: -counts_t[n])[:args.n_perts]
    print(f"  target controls: {t_ctrl.shape[0]:,} cells")
    print(f"  scoring {len(usable)} perturbations present in both, "
          f"{min(counts_t[n] for n in usable)}-{max(counts_t[n] for n in usable)}"
          f" cells each")
    if not usable:
        raise SystemExit("no shared perturbations with enough cells")

    # gate, from the target's control cells alone
    cmean_full = ctrl_cpm_mean(t_ctrl)
    gate = cmean_full > M.MIN_CPM
    print(f"  genes passing the 5-CPM gate: {int(gate.sum()):,}")
    ctrl_dense = M.cpm(np.asarray(t_ctrl[:, gate].todense(), dtype=np.float32))
    ctrl_mean = ctrl_dense.mean(axis=0)
    gated_sym = sym[gate]
    t_ctrl_mean_log = log1p_cp10k(t_ctrl)

    # --- reference ------------------------------------------------------
    print("\nreference DE on measured target cells")
    real: dict[str, dict] = {}
    for k, name in enumerate(usable):
        cells = Xt[lt == name]
        sig, lfc, _ = de_table(cells, ctrl_dense, ctrl_mean, gate)
        own = np.flatnonzero(gated_sym == name)
        if own.size:
            sig[own[0]] = False
        real[name] = {"sig": sig, "lfc": lfc, "n": int(sig.sum()),
                      "ncells": cells.shape[0]}
        if (k + 1) % 20 == 0:
            print(f"  {k + 1}/{len(usable)}", flush=True)
    n_real = np.array([real[n]["n"] for n in usable])
    print(f"  significant genes per target: median {np.median(n_real):.0f}, "
          f"IQR {np.percentile(n_real, 25):.0f}-{np.percentile(n_real, 75):.0f}, "
          f"range {n_real.min()}-{n_real.max()}")

    # the source's own response sizes, which is what a submission may use
    print("\nsource response sizes (the deployable count estimate)")
    s_ctrl_sp = Xs[ls == CONTROL]
    s_cmean = ctrl_cpm_mean(s_ctrl_sp)
    s_gate = s_cmean > M.MIN_CPM
    s_ctrl_dense = M.cpm(np.asarray(s_ctrl_sp[:, s_gate].todense(),
                                    dtype=np.float32))
    s_ctrl_m = s_ctrl_dense.mean(axis=0)
    src_counts = []
    for name in usable[: min(40, len(usable))]:
        rows = ls == name
        if rows.sum() < 10:
            continue
        sig, _, _ = de_table(Xs[rows], s_ctrl_dense, s_ctrl_m, s_gate)
        src_counts.append(int(sig.sum()))
    src_counts = np.array(src_counts) if src_counts else n_real.copy()
    print(f"  median {np.median(src_counts):.0f}, "
          f"range {src_counts.min()}-{src_counts.max()} "
          f"(from {src_counts.size} source perturbations)")

    # --- prediction arms ------------------------------------------------
    models: dict[str, np.ndarray] = {
        "naive transfer": np.stack([src_effect[n] for n in usable])}

    print("\nfitting ContextTransfer on the source screen + target controls only")
    try:
        line = source_cellline(Xs, ls, sym, seed=args.seed)
        bank = SourceBank.build([line], np.ones(1))
        hp = Hyper(**JURKAT_FOLD_HYPER)
        model = ContextTransferModel(hp).fit([line], t_ctrl_mean_log, sym,
                                             bank=bank)
        models["ContextTransfer"] = model.predict(np.array(usable))
        print(f"  fitted on {line.names.size:,} source perturbations")
    except Exception as exc:                       # noqa: BLE001
        print(f"  SKIPPED: {type(exc).__name__}: {exc}")

    oracle_k = n_real.copy()
    print(f"\ncount targets: oracle median {np.median(oracle_k):.0f}")

    def run(full_delta: np.ndarray, label: str,
            k_per_pert: np.ndarray | None) -> dict:
        rows = []
        for i, name in enumerate(usable):
            eff = full_delta[i].copy()
            if k_per_pert is not None:
                # Sparsify only *within* the gated set.  Genes below the 5-CPM
                # gate are never tested, but they do carry counts, so zeroing
                # them would change the generated cells' library sizes and
                # therefore every gene's CPM -- making the arms differ by more
                # than the one thing under test.  They are left alone.
                keep_idx = np.flatnonzero(gate)
                scores = np.abs(eff[keep_idx])
                k = int(min(max(k_per_pert[i], 1), keep_idx.size))
                top = keep_idx[np.argpartition(-scores, k - 1)[:k]]
                mask = ~gate.copy()
                mask[top] = True
                eff = np.where(mask, eff, 0.0)
            cells = sample_cells(t_ctrl, eff, real[name]["ncells"], rng,
                                 control_mean=t_ctrl_mean_log)
            sig, lfc, p = de_table(cells, ctrl_dense, ctrl_mean, gate)
            own = np.flatnonzero(gated_sym == name)
            if own.size:
                sig[own[0]] = False
            r, rl = real[name]["sig"], real[name]["lfc"]
            inter = int((sig & r).sum())
            rows.append({
                "pert": name, "n_pred": int(sig.sum()), "n_real": int(r.sum()),
                "recall": inter / max(int(r.sum()), 1),
                "jaccard": M.sig_jaccard(sig, r),
                "fidelity": M.direction_fidelity(sig, lfc, r, rl),
                "reach": M.direction_reach(sig, p, p, lfc, r, rl, gated_sym),
                "nmae": M.lfc_nmae(lfc, rl, r)})
        return {"label": label, "rows": rows}

    arms = []
    for mname, full_delta in models.items():
        deployable_k = rank_matched_targets(full_delta[:, gate], src_counts)
        print(f"  {mname}: deployable count median "
              f"{np.median(deployable_k):.0f}")
        arms.append(run(full_delta, f"{mname}: uncalibrated", None))
        arms.append(run(full_delta, f"{mname}: calibrated (source counts)",
                        deployable_k))
        arms.append(run(full_delta, f"{mname}: oracle counts [not submittable]",
                        oracle_k))

    def agg(rows, key):
        v = np.array([r[key] for r in rows], dtype=float)
        v = v[np.isfinite(v)]
        return float(np.mean(v)) if v.size else float("nan")

    print("\n" + "=" * 94)
    print("K562 -> Jurkat, zero-shot: source effects plus the target's control "
          "cells only")
    print(f"{len(usable)} perturbations, cell counts matched between reference "
          f"and prediction\n")
    print(f"  {'arm':<46}{'calls':>8}{'recall':>9}{'jaccard':>9}"
          f"{'fidelity':>10}{'reach':>8}{'nmae':>8}")
    print("  " + "-" * 92)
    out = {"n_perts": len(usable), "n_genes_gated": int(gate.sum()),
           "n_real_median": float(np.median(n_real)),
           "source_counts": src_counts.tolist(), "arms": {}}
    for a in arms:
        r = a["rows"]
        rec = {k: agg(r, k) for k in
               ("n_pred", "recall", "jaccard", "fidelity", "reach", "nmae")}
        out["arms"][a["label"]] = {"mean": rec, "rows": r}
        print(f"  {a['label']:<46}{rec['n_pred']:>8.0f}{rec['recall']:>9.3f}"
              f"{rec['jaccard']:>9.4f}{rec['fidelity']:>10.4f}"
              f"{rec['reach']:>8.4f}{rec['nmae']:>8.4f}")
    print(f"  {'reference (measured, = the 1.0 anchor)':<46}"
          f"{np.median(n_real):>8.0f}{1.0:>9.3f}{1.0:>9.4f}{1.0:>10.4f}"
          f"{1.0:>8.4f}{0.0:>8.4f}")

    # the number this module exists to produce
    b, r_anchor = 0.029, 0.399          # published jaccard anchors, midpoints
    print("\n" + "=" * 94)
    print("The unmeasured number, measured: recall retained when the call count "
          "is corrected")
    for mname in models:
        unc = out["arms"][f"{mname}: uncalibrated"]["mean"]
        cal = out["arms"][f"{mname}: calibrated (source counts)"]["mean"]
        orc = out["arms"][f"{mname}: oracle counts [not submittable]"]["mean"]
        print(f"\n  {mname}")
        for label, v in (("uncalibrated", unc), ("deployable", cal),
                         ("oracle counts", orc)):
            print(f"    {label:<16} {v['n_pred']:>6.0f} calls  "
                  f"recall {v['recall']:>6.1%}  raw J {v['jaccard']:.4f}  "
                  f"scaled J {(v['jaccard'] - b) / (r_anchor - b):+.3f}  "
                  f"composite {((v['jaccard'] - b) / (r_anchor - b)) / 6:+.4f}")
        out["arms"][f"{mname}: uncalibrated"]["scaled_jaccard"] = \
            (unc["jaccard"] - b) / (r_anchor - b)
        out["arms"][f"{mname}: calibrated (source counts)"]["scaled_jaccard"] = \
            (cal["jaccard"] - b) / (r_anchor - b)
    print("\n  Scaling uses the published anchors for the official contexts, so "
          "it is an\n  indicative translation onto the leaderboard's units, not "
          "a leaderboard score:\n  these are different cells, a different panel "
          "and fewer cells per group.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
