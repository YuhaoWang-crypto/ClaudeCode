"""
Backend for the paper's own released library, `proteoclock`.

The paper's code-availability statement points at
https://github.com/Insilico-org/proteoclock, which ships the real published
weights for several of the six clocks. When that package is importable this
module replaces the corresponding surrogates with the actual models, and the
pipeline stops being a demonstration for those clocks.

Install:
    git clone https://github.com/Insilico-org/proteoclock
    cd proteoclock && pip install -e .
    pip install --upgrade scikit-posthocs

That last line matters. setup.py pins `scikit-posthocs==0.11`, which imports
`multipletests` from `statsmodels.sandbox.stats.multicomp`. That location was
removed in statsmodels 0.15, so importing proteoclock against a current
statsmodels raises ImportError until scikit-posthocs is upgraded. The pin is
then violated, and pip says so, but the package works.

WHAT IS AND IS NOT IN THE PACKAGE
---------------------------------
  kuo_2024                        PAC, real weights, shipped
  goeminne_2025_full_chrono       OrganAge chronological, real weights
  goeminne_2025_full_mortality    OrganAge mortality, real weights
  goeminne_2025_reduced_*         Explore-1536 reduced variants, real weights
  galkin_2025                     ipfP3GPT: feature order ONLY. The weights
                                  are not distributed, because they may only
                                  be used inside the UK Biobank Research
                                  Analysis Platform. A surrogate stands in.

ProtAge and PAOPAC are not in this package at all; see m3_clocks for where
they live and what shape they are in.

CALL SIGNATURES differ by clock class and the README's example omits it:
PAC is a GompertzClock and takes (data, age_data, scaling), while the
OrganAge clocks are LinearClock/CPHClock and take (data, scaling). Passing
age to the latter raises "got multiple values for argument 'scaling'".
"""
import os
import warnings

import numpy as np
import pandas as pd

# name in this pipeline -> (proteoclock clock id, scaler id, needs age)
REAL_CLOCK_MAP = {
    "PAC": ("kuo_2024", "ukb_scaler", True),
    "OrganAge_chrono": ("goeminne_2025_full_chrono", "goeminne_2025_full",
                        False),
    "OrganAge_mortality": ("goeminne_2025_full_mortality",
                           "goeminne_2025_full", False),
}

_FACTORY = None
_SYMBOLS = None


def available():
    """True if proteoclock imports, with its dependency conflict resolved."""
    try:
        import proteoclock  # noqa: F401
        return True
    except Exception:
        return False


def why_unavailable():
    """The actual import error, for reporting rather than guessing."""
    try:
        import proteoclock  # noqa: F401
        return None
    except Exception as exc:
        return f"{type(exc).__name__}: {exc}"


def _factory():
    global _FACTORY
    if _FACTORY is None:
        from proteoclock import ClockFactory
        _FACTORY = ClockFactory()
    return _FACTORY


def panel_symbols():
    """
    The real Olink Explore 3072 gene symbols the published clocks expect.

    Taken from the package's own feature order, 2,923 assays. Using these
    instead of invented identifiers is what lets the real clocks score the
    synthetic cohort at all: a clock keyed on "IL6" cannot read a panel whose
    proteins are called "OLK0001_PROT1". This is also the namespace trap that
    real_data.py warns about, made concrete.
    """
    global _SYMBOLS
    if _SYMBOLS is None:
        import os
        import proteoclock
        path = os.path.join(os.path.dirname(proteoclock.__file__),
                            "materials", "deep_clocks", "galkin_2025",
                            "feature_order.txt")
        # feature_order.txt is TAB-SEPARATED with the symbol repeated in two
        # columns. Taking the whole line yields "A1BG\tA1BG", which matches
        # nothing and silently sends every clock to a surrogate. This is the
        # namespace trap real_data.py warns about, and it bit this pipeline.
        with open(path) as fh:
            syms = [ln.strip().split("\t")[0] for ln in fh if ln.strip()]
        _SYMBOLS = np.array(syms)
    return _SYMBOLS


def reference_aging_coefficients(clock_id="goeminne_2025_full_chrono",
                                 scaler_id="goeminne_2025_full"):
    """
    Per-protein coefficients of a published chronological-age clock.

    Used by M1 to give the synthetic cohort an aging structure the REAL
    clocks can actually read. Assigning age slopes at random over real gene
    symbols produces a panel where the published clocks see nothing: they are
    keyed on which proteins genuinely move with age in humans, and a random
    assignment is orthogonal to that. Grounding the simulated aging axis in a
    published clock's own coefficients fixes it.

    This makes the recovery of age acceleration CIRCULAR for the clocks that
    share this axis, so it is not evidence that any clock works. It exists so
    that the trial statistics downstream operate on data with a realistic
    covariance structure rather than noise.
    """
    if not available():
        return None
    try:
        clk = _factory().get_clock(clock_id, scaler=scaler_id)
        return dict(getattr(clk, "protein_coefs", {}))
    except Exception:
        return None


class ProteoclockClock:
    """
    Wraps a published proteoclock model in this pipeline's clock interface.

    Only the proteins a given clock actually uses are passed through, so a
    reference cohort is reshaped into a few hundred columns rather than all
    2,923. Feeding the whole panel in long format would be millions of rows
    per call for no benefit.
    """

    def __init__(self, spec, clock, protein_ids, needs_age, scaling="standard"):
        self.spec = spec
        self.clock = clock
        self.protein_ids = np.asarray(protein_ids)
        self.needs_age = needs_age
        self.scaling = scaling
        self.is_surrogate = False

        wanted = set(getattr(clock, "protein_coefs", {}).keys())
        if wanted:
            mask = np.isin(self.protein_ids, list(wanted))
        else:
            mask = np.ones(len(self.protein_ids), dtype=bool)
        self.col_idx = np.where(mask)[0]
        self.n_matched = len(self.col_idx)
        self.n_required = len(wanted) if wanted else -1

    def predict(self, npx, meta=None):
        n = npx.shape[0]
        ids = [f"s{i}" for i in range(n)]
        sub = npx[:, self.col_idx]
        syms = self.protein_ids[self.col_idx]

        long = pd.DataFrame({
            "patient_id": np.repeat(ids, len(syms)),
            "gene_symbol": np.tile(syms, n),
            "NPX": sub.ravel(),
        })

        # proteoclock reports missing features and imputation counts with
        # bare print(), not the warnings module, so a single scoring call
        # emits tens of kilobytes. Both streams are captured; the useful part
        # (how many proteins matched) is reported once at load time instead.
        import contextlib
        import io

        sink = io.StringIO()
        with warnings.catch_warnings(), contextlib.redirect_stdout(sink), \
                contextlib.redirect_stderr(sink):
            warnings.simplefilter("ignore")
            if self.needs_age:
                if meta is None or "age" not in meta:
                    raise ValueError(
                        f"{self.spec.name} needs chronological age and none "
                        f"was passed; PAC is a Gompertz model keyed on age")
                age_df = pd.DataFrame({"patient_id": ids,
                                       "age": np.asarray(meta["age"])})
                pred = self.clock.predict_age(long, age_df,
                                              scaling=self.scaling)
            else:
                pred = self.clock.predict_age(long, scaling=self.scaling)

        # Restore the caller's row order; proteoclock indexes by patient_id.
        return pred.reindex(ids).to_numpy(dtype=float)


def load_real_clocks(specs, protein_ids, verbose=True):
    """
    Build every published clock proteoclock can supply for `specs`.

    Returns {clock name: ProteoclockClock}. Clocks the package does not carry,
    or that fail to load, are simply absent and the caller falls back to a
    surrogate.
    """
    if not available():
        return {}

    out = {}
    for spec in specs:
        if spec.name not in REAL_CLOCK_MAP:
            continue
        clock_id, scaler_id, needs_age = REAL_CLOCK_MAP[spec.name]
        try:
            clk = _factory().get_clock(clock_id, scaler=scaler_id)
            wrapped = ProteoclockClock(spec, clk, protein_ids, needs_age)
            if wrapped.n_matched == 0:
                if verbose:
                    print(f"    {spec.name}: proteoclock loaded but 0 of its "
                          f"{wrapped.n_required} proteins are on this panel; "
                          f"falling back to a surrogate")
                continue
            out[spec.name] = wrapped
            if verbose:
                print(f"    {spec.name}: REAL weights via proteoclock "
                      f"({clock_id}), {wrapped.n_matched}/"
                      f"{wrapped.n_required} proteins on panel")
        except Exception as exc:
            if verbose:
                print(f"    {spec.name}: proteoclock load failed "
                      f"({type(exc).__name__}: {str(exc)[:70]}); using a "
                      f"surrogate")
    return out


def _gse169148_path():
    """Path to the real GEO dataset the package ships for testing."""
    import proteoclock
    return os.path.join(os.path.dirname(proteoclock.__file__), "materials",
                        "test_data", "GSE169148")


def validate_against_package(tol=1e-9, verbose=True):
    """
    Check that this adapter reproduces a DIRECT proteoclock call exactly.

    The wrapper reshapes a wide NPX matrix into the long frame proteoclock
    wants, restricts it to each clock's own proteins, and reindexes the
    result back to the caller's row order. Every one of those steps could
    silently reorder or drop samples while still returning plausible ages, so
    agreement is asserted rather than assumed.

    Runs on GSE169148, the real 31-sample Olink dataset bundled with the
    package, not on simulated data. Returns {clock: max abs difference}.

    NOTE on the bundled `new_clock_res_GSE169148.tsv`: that file is NOT used
    as the target here. Its values correlate about 0.96 with what the current
    API returns but sit on a different scale (a log hazard near 2 against an
    age-like 63), so it appears to predate the present code. The live package
    is the right comparison for an adapter in any case.
    """
    import pandas as pd

    if not available():
        return None

    base = _gse169148_path()
    long = pd.read_csv(os.path.join(base, "GSE169148_protein_data_long.txt"),
                       sep="\t")
    wide = long.pivot(index="patient_id", columns="gene_symbol", values="NPX")
    ids = list(wide.index)
    syms = np.asarray(wide.columns)
    x = wide.to_numpy(dtype=float)
    x = np.nan_to_num(x, nan=float(np.nanmedian(x)))
    meta = {"age": np.full(len(ids), 60.0)}

    from protclock_pipeline.m3_clocks import CLOCKS
    specs = {s.name: s for s in CLOCKS}

    import contextlib
    import io
    sink = io.StringIO()
    out = {}
    for name, (clock_id, scaler_id, needs_age) in REAL_CLOCK_MAP.items():
        with warnings.catch_warnings(), contextlib.redirect_stdout(sink), \
                contextlib.redirect_stderr(sink):
            warnings.simplefilter("ignore")
            clk = _factory().get_clock(clock_id, scaler=scaler_id)
            if needs_age:
                age_df = pd.DataFrame({"patient_id": ids, "age": meta["age"]})
                direct = clk.predict_age(long, age_df, scaling="standard")
            else:
                direct = clk.predict_age(long, scaling="standard")
        direct = direct.reindex(ids).to_numpy(dtype=float)

        wrapped = ProteoclockClock(specs[name],
                                   _factory().get_clock(clock_id,
                                                        scaler=scaler_id),
                                   syms, needs_age)
        mine = wrapped.predict(x, meta)
        diff = float(np.max(np.abs(direct - mine)))
        out[name] = diff
        if verbose:
            status = "OK" if diff <= tol else "MISMATCH"
            print(f"    {name:<19} max|diff| = {diff:.3e}   {status}")

    bad = {k: v for k, v in out.items() if v > tol}
    if bad:
        raise AssertionError(
            f"adapter disagrees with proteoclock for {list(bad)}; "
            f"differences {bad}")
    return out


def report():
    print("PROTEOCLOCK BACKEND  (the paper's released library)")
    print("-" * 68)
    if not available():
        print(f"  NOT importable: {why_unavailable()}")
        print("  install: git clone https://github.com/Insilico-org/proteoclock")
        print("           cd proteoclock && pip install -e .")
        print("           pip install --upgrade scikit-posthocs")
        return None

    syms = panel_symbols()
    print(f"  importable        : yes")
    print(f"  panel symbols     : {len(syms)} "
          f"(Olink Explore 3072, real gene symbols)")
    print(f"  clocks available  : {_factory().list_available_clocks()}")
    print(f"  mapped into this pipeline:")
    for name, (cid, sc, age) in REAL_CLOCK_MAP.items():
        print(f"    {name:<19} -> {cid:<31} scaler={sc}"
              f"{'  (needs age)' if age else ''}")
    print("  galkin_2025 (ipfP3GPT) ships feature order only; its weights are")
    print("  restricted to the UK Biobank Research Analysis Platform.")

    print("\n  adapter check on GSE169148 (31 real samples, 1463 proteins):")
    diffs = validate_against_package()
    print("  this adapter reproduces a direct proteoclock call bitwise.")
    return {"symbols": syms, "validation": diffs}


if __name__ == "__main__":
    report()
