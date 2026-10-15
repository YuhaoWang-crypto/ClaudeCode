"""
Template: run the state-restoration pipeline on a NEW indication.

Copy this file, fill in the CONFIG block, run:
    python3 run_state_pipeline.py            (from the repository root)

Two input routes:
  A  a GEO series (expression + sample characteristics are read from the
     series matrix; you map characteristics to the pipeline's columns);
  B  two CSV files: expression (features x samples, log scale) and a sample
     sheet with columns  sample, subject, kind (ref|state), time, arm,
     outcome (optional), group (optional).

Everything below CONFIG is generic.
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from netjac import statepipe as sp          # noqa: E402
from netjac.data import geo                 # noqa: E402

CONFIG = {
    # ---- route A: GEO --------------------------------------------------
    "gse": None,                      # e.g. "GSE85034"
    "gpl": None,                      # e.g. "GPL10558" (probe -> gene symbol)
    "log2": False,                    # True if the matrix is linear (MAS5): log2(max(x, 16))
    # map GEO sample metadata (DataFrame M) to pipeline columns
    "subject": lambda M: M["title"].str.extract(r"Subject (\d+)")[0],
    "kind": lambda M: np.where(M["timepoint"] == "NL", "ref", "state"),
    "time": lambda M: M["timepoint"].map({"LS": 0, "WK 1": 1, "WK 2": 2, "WK 4": 4, "WK 16": 16}),
    "arm": lambda M: M["treatment"],
    "outcome": None,                  # function M -> per-sample bool (constant within subject) or None
    "group": None,                    # function M -> per-sample replication group or None
    # ---- route B: CSV ----------------------------------------------------
    "expr_csv": None,                 # features x samples
    "samples_csv": None,
    # ---- modules ---------------------------------------------------------
    # curated modules {name: [genes]}; None -> data-driven auto_modules(k)
    "modules": None,
    "auto_k": 6,
    "control_arms": ("placebo", "Placebo", "control"),
}


def load(cfg):
    if cfg["gse"]:
        X, M = geo.series_matrix(cfg["gse"])
        if cfg["log2"]:
            X = np.log2(X.clip(lower=16))
        G = geo.collapse_to_genes(geo.quantile_normalise(X), geo.gpl_symbols(cfg["gpl"]))
        subj = np.asarray(cfg["subject"](M))
        kind = np.asarray(cfg["kind"](M))
        time = np.asarray(cfg["time"](M), float)
        arm = dict(zip(subj, np.asarray(cfg["arm"](M)))) if cfg["arm"] else {}
        out = dict(zip(subj, np.asarray(cfg["outcome"](M)))) if cfg["outcome"] else {}
        grp = dict(zip(subj, np.asarray(cfg["group"](M)))) if cfg["group"] else {}
        keep = pd.notna(subj) & (kind == "ref") | pd.notna(time)
        Xs = G[M.index[keep]].T.values
        feats = np.array(G.index)
        subj, kind, time = subj[keep], kind[keep], time[keep]
    else:
        E = pd.read_csv(cfg["expr_csv"], index_col=0)
        Sm = pd.read_csv(cfg["samples_csv"])
        Xs = E[Sm["sample"]].T.values
        feats = np.array(E.index)
        subj, kind, time = Sm.subject.astype(str).values, Sm.kind.values, Sm.time.astype(float).values
        first = Sm.drop_duplicates("subject").set_index("subject")
        arm = first["arm"].to_dict() if "arm" in Sm else {}
        out = first["outcome"].astype(bool).to_dict() if "outcome" in Sm else {}
        grp = first["group"].to_dict() if "group" in Sm else {}
    return sp.StateData(X=Xs, features=feats, subject=subj, kind=kind, time=time,
                        arm=arm, outcome=out, group=grp, control_arms=cfg["control_arms"])


if __name__ == "__main__":
    D = load(CONFIG)
    modules = CONFIG["modules"]
    if modules is None:
        modules = sp.auto_modules(D, sp.steady_state(D), k=CONFIG["auto_k"])
    card = sp.run(D, modules)
    sp.print_card(card, CONFIG["gse"] or CONFIG["expr_csv"])
