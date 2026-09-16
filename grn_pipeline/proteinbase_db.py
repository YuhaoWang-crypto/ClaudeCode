"""
Proteinbase data layer — a local, queryable mirror of the Adaptyv Bio
Proteinbase release (experimental protein-design data with measured Kd).

WHY THIS EXISTS
---------------
Every binding claim in this pipeline (M5–M10) rested on a handful of ChEMBL
small-molecule IC50s: M10 validates the Boltz screen against n=5 measured
potencies. That is far too small to calibrate a co-folding score.

Proteinbase is the missing truth set. It is the only open release that pairs
*designed protein binders* with (a) Kd measured in one lab under standardised
protocols, (b) the computational scores that were used to pick them, and
crucially (c) the NEGATIVES — designs that expressed but did not bind, which
are the denominator every hit-rate claim needs and which almost nobody
publishes.

Note "one lab" is not "one instrument": Kd here comes from BOTH SPR and BLI
(roughly half each), so the `platform` column travels with every measurement
and anything sensitive to a cross-platform offset must stratify on it.

Snapshot 28_01_2026 contains:
    5,253 designs · 2,630 (design, target) assay pairs · 1,332 Kd measurements
    2,160 measured non-binders vs 470 binders · 3,796 Boltz-2 complex scores

DATA MODEL
----------
The raw CSV is one row per design with all measurements crammed into a nested
`evaluations` JSON array. That shape is unusable for analysis, so it is
normalised into two tidy tables:

  designs   one row per design      — sequence-level: sequence, design method,
                                      author, ESMFold pLDDT, ProteinMPNN score,
                                      expression outcome and yield.
  complexes one row per (design, target) — interface-level: Boltz-2 confidence
                                      metrics AND the experimental outcome
                                      (Kd, kon, koff, binder label).

Splitting on that axis is lossless here: no design in the snapshot carries the
same metric against two different targets (verified in `qc_report`), while 113
designs *are* assayed against two targets — so (design, target) is the correct
key for anything interface-level, and collapsing to the design alone would
silently merge two different experiments.

The derived tables are committed gzipped under `data/proteinbase/`, so the
analysis runs offline. The 40 MB raw CSV is cached (gitignored) and only
re-fetched with `--rebuild`.

LICENCE / ATTRIBUTION (required, do not strip)
---------------------------------------------
Proteinbase data is published by Adaptyv Bio under the Open Data Commons
Attribution Licence (ODC-By). Any work using these tables must carry:

    "This work used Proteinbase by Adaptyv Bio under ODC-BY license"

Usage
-----
    python3 -m grn_pipeline.proteinbase_db            # QC report from cache
    python3 -m grn_pipeline.proteinbase_db --rebuild  # re-fetch + rebuild
    python3 -m grn_pipeline.proteinbase_db --sqlite   # emit a queryable .db

    from grn_pipeline import proteinbase_db as pb
    designs, complexes = pb.load()
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import statistics
import subprocess
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------- constants

SNAPSHOT = "28_01_2026"
RAW_URL = f"https://storage.proteinbase.com/proteinbase_all_data_{SNAPSHOT}.csv"
ATTRIBUTION = "This work used Proteinbase by Adaptyv Bio under ODC-BY license"
LICENSE = "ODC-By 1.0 (Open Data Commons Attribution)"

_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = _ROOT / "data" / "proteinbase"
CACHE_DIR = _ROOT / "figures" / "_pb_cache"
DESIGNS_CSV = DATA_DIR / "pb_designs.csv.gz"
COMPLEXES_CSV = DATA_DIR / "pb_complexes.csv.gz"
PROVENANCE_JSON = DATA_DIR / "provenance.json"
RAW_CSV = CACHE_DIR / f"proteinbase_all_data_{SNAPSHOT}.csv"

# Sequence-level computational metrics -> one value per design.
DESIGN_METRICS = (
    "esmfold_plddt",
    "proteinmpnn_score",
    "redesigned_proteinmpnn_score",
    "proteinmpnn_seq_recovery",
    "molecular_weight",
    "isoelectric_point",
    "ted_confidence",
)
# Interface-level computational metrics -> one value per (design, target).
COMPLEX_METRICS = (
    "boltz2_ipsae",
    "boltz2_min_ipsae",
    "boltz2_iptm",
    "boltz2_ptm",
    "boltz2_plddt",
    "boltz2_complex_plddt",
    "boltz2_complex_iplddt",
    "boltz2_complex_pde",
    "boltz2_lis",
    "boltz2_pdockq",
    "boltz2_pdockq2",
    "shape_complimentarity_boltz2_binder_ss",
)

DESIGN_COLUMNS = (
    ["design_id", "name", "author", "design_method", "design_class",
     "is_control", "seq_len", "sequence", "expressed", "expression_yield_ug_ml"]
    + list(DESIGN_METRICS)
)
COMPLEX_COLUMNS = (
    ["design_id", "target", "design_method", "binds", "binding_strength",
     "kd_M", "kd_nM", "pkd", "kd_n_replicates", "kd_log10_spread",
     "kon_M1s1", "koff_s1", "kon_n_replicates", "koff_n_replicates",
     "fit_model", "platform"]
    + list(COMPLEX_METRICS)
)

# Binder labels. "None" is Proteinbase's label for a measured NON-binder;
# it is data, not a missing value, and dropping it destroys the negatives.
BINDER_LABELS = ("Strong", "Medium", "Weak")
NONBINDER_LABEL = "None"


# ------------------------------------------------------------------- fetch

def fetch_raw(force: bool = False, retries: int = 4) -> Path:
    """Download the snapshot CSV into the (gitignored) cache, with backoff."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if RAW_CSV.exists() and not force:
        return RAW_CSV
    delay = 2
    for attempt in range(1, retries + 1):
        tmp = RAW_CSV.with_suffix(".part")
        r = subprocess.run(
            ["curl", "-sSL", "--fail", "-o", str(tmp), RAW_URL],
            capture_output=True, text=True,
        )
        if r.returncode == 0 and tmp.exists() and tmp.stat().st_size > 1_000_000:
            tmp.replace(RAW_CSV)
            return RAW_CSV
        tmp.unlink(missing_ok=True)
        if attempt == retries:
            raise RuntimeError(
                f"could not fetch {RAW_URL} after {retries} attempts: "
                f"{r.stderr.strip()[:200]}"
            )
        time.sleep(delay)
        delay *= 2
    raise RuntimeError("unreachable")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ------------------------------------------------------------------- parse

def _gmean(values):
    """Geometric mean — the right average for a dissociation constant.

    Kd is log-distributed and spans 5 orders of magnitude here; an arithmetic
    mean of replicates would be dragged to the weakest one.
    """
    return 10 ** statistics.fmean(math.log10(v) for v in values if v > 0)


def _iter_raw(raw_path: Path):
    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    with open(raw_path, encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            yield row


def build_tables(raw_path: Path | None = None):
    """Normalise the raw nested CSV into (designs, complexes) row lists."""
    raw_path = raw_path or fetch_raw()
    designs, complexes = [], []

    for row in _iter_raw(raw_path):
        pid = row["id"]
        name = row.get("name") or ""
        method = row.get("designMethod") or ""
        seq = row.get("sequence") or ""

        d_num: dict[str, list[float]] = {}
        c_num: dict[tuple[str, str], dict[str, list[float]]] = {}
        kd: dict[str, list[float]] = {}
        kon: dict[str, list[float]] = {}
        koff: dict[str, list[float]] = {}
        label: dict[str, str] = {}
        binds: dict[str, bool] = {}
        fit: dict[str, str] = {}
        # Which biophysical platform produced the trace. Kd in this snapshot
        # comes from BOTH SPR and BLI, and cross-platform Kd offsets are a
        # known systematic, so the platform has to travel with the number.
        platform: dict[str, set[str]] = {}
        design_class = ""
        expressed = None
        yields: list[float] = []

        for ev in json.loads(row["evaluations"] or "[]"):
            metric, target = ev.get("metric"), ev.get("target")
            value, etype = ev.get("value"), ev.get("type")
            numeric = isinstance(value, (int, float)) and not isinstance(value, bool)

            if etype == "computational":
                if metric == "design_class":
                    design_class = value or ""
                elif numeric and metric in DESIGN_METRICS:
                    d_num.setdefault(metric, []).append(float(value))
                elif numeric and metric in COMPLEX_METRICS:
                    key = (pid, target or "")
                    c_num.setdefault(key, {}).setdefault(metric, []).append(float(value))

            elif etype == "experimental":
                tgt = target or ""
                if metric == "kd" and numeric:
                    kd.setdefault(tgt, []).append(float(value))
                elif metric == "kon" and numeric:
                    kon.setdefault(tgt, []).append(float(value))
                elif metric == "koff" and numeric:
                    koff.setdefault(tgt, []).append(float(value))
                elif metric == "binding_strength":
                    label[tgt] = value
                elif metric == "binding":
                    binds[tgt] = bool(value)
                elif metric == "selected_binding_fit_model":
                    fit[tgt] = value
                elif metric == "spr_kinetic_curves":
                    platform.setdefault(tgt, set()).add("SPR")
                elif metric == "bli_kinetic_curves":
                    platform.setdefault(tgt, set()).add("BLI")
                elif metric == "expressed":
                    expressed = bool(value)
                elif metric == "expression-yield" and numeric:
                    yields.append(float(value))

        designs.append({
            "design_id": pid,
            "name": name,
            "author": row.get("author") or "",
            "design_method": method,
            "design_class": design_class,
            # Adaptyv's own reference/benchmark constructs, not submitted designs.
            "is_control": int("control" in name.lower()),
            "seq_len": len(seq),
            "sequence": seq,
            "expressed": "" if expressed is None else int(expressed),
            "expression_yield_ug_ml": _fmt(statistics.fmean(yields)) if yields else "",
            **{m: _fmt(statistics.fmean(v)) for m, v in d_num.items()},
        })

        # One complex row per (design, target) seen by EITHER side, so a design
        # scored but never assayed, and one assayed but never scored, both survive.
        targets = set(kd) | set(label) | set(binds) | {t for (_, t) in c_num}
        for tgt in sorted(t for t in targets if t):
            cm = c_num.get((pid, tgt), {})
            reps = kd.get(tgt, [])
            kd_g = _gmean(reps) if reps else None
            spread = (max(map(math.log10, reps)) - min(map(math.log10, reps))
                      if len(reps) > 1 else "")
            lbl = label.get(tgt, "")
            complexes.append({
                "design_id": pid,
                "target": tgt,
                "design_method": method,
                "binds": "" if tgt not in binds else int(binds[tgt]),
                "binding_strength": lbl,
                "kd_M": _fmt(kd_g, 6) if kd_g else "",
                "kd_nM": _fmt(kd_g * 1e9, 4) if kd_g else "",
                "pkd": _fmt(-math.log10(kd_g), 4) if kd_g else "",
                "kd_n_replicates": len(reps),
                "kd_log10_spread": _fmt(spread, 4) if spread != "" else "",
                "kon_M1s1": _fmt(_gmean(kon[tgt]), 6) if kon.get(tgt) else "",
                "koff_s1": _fmt(_gmean(koff[tgt]), 6) if koff.get(tgt) else "",
                "kon_n_replicates": len(kon.get(tgt, [])),
                "koff_n_replicates": len(koff.get(tgt, [])),
                "fit_model": fit.get(tgt, ""),
                "platform": "+".join(sorted(platform.get(tgt, ()))),
                **{m: _fmt(statistics.fmean(v)) for m, v in cm.items()},
            })

    return designs, complexes


def _fmt(x, sig: int = 6):
    if x is None or x == "":
        return ""
    return f"{float(x):.{sig}g}"


# -------------------------------------------------------------- write/load

def _write_gz(path: Path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in columns})


def rebuild(force_fetch: bool = False):
    """Re-fetch the snapshot and regenerate the committed derived tables."""
    raw = fetch_raw(force=force_fetch)
    designs, complexes = build_tables(raw)
    _write_gz(DESIGNS_CSV, DESIGN_COLUMNS, designs)
    _write_gz(COMPLEXES_CSV, COMPLEX_COLUMNS, complexes)
    PROVENANCE_JSON.parent.mkdir(parents=True, exist_ok=True)
    PROVENANCE_JSON.write_text(json.dumps({
        "source": "Proteinbase (Adaptyv Bio)",
        "url": RAW_URL,
        "snapshot": SNAPSHOT,
        "license": LICENSE,
        "attribution": ATTRIBUTION,
        "raw_bytes": raw.stat().st_size,
        "raw_sha256": _sha256(raw),
        "n_designs": len(designs),
        "n_complexes": len(complexes),
        "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "built_by": "grn_pipeline.proteinbase_db",
    }, indent=2) + "\n")
    return designs, complexes


def _read_gz(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def load(rebuild_if_missing: bool = True):
    """Return (designs, complexes) as lists of dicts; strings stay strings."""
    if DESIGNS_CSV.exists() and COMPLEXES_CSV.exists():
        return _read_gz(DESIGNS_CSV), _read_gz(COMPLEXES_CSV)
    if not rebuild_if_missing:
        raise FileNotFoundError(f"{DESIGNS_CSV} missing; run --rebuild")
    return rebuild()


def num(row, key):
    """Parse a numeric cell; '' -> None (missing is not zero)."""
    v = row.get(key, "")
    if v == "" or v is None:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def provenance():
    if PROVENANCE_JSON.exists():
        return json.loads(PROVENANCE_JSON.read_text())
    return {}


# ------------------------------------------------------------------ sqlite

def to_sqlite(path: Path | None = None):
    """Build a queryable SQLite db (gitignored) from the derived tables."""
    import sqlite3

    path = Path(path or CACHE_DIR / "proteinbase.db")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.unlink(missing_ok=True)
    designs, complexes = load()

    con = sqlite3.connect(path)
    con.execute(f"CREATE TABLE designs ({','.join(DESIGN_COLUMNS)})")
    con.execute(f"CREATE TABLE complexes ({','.join(COMPLEX_COLUMNS)})")
    con.executemany(
        f"INSERT INTO designs VALUES ({','.join('?' * len(DESIGN_COLUMNS))})",
        [[d.get(c, "") for c in DESIGN_COLUMNS] for d in designs])
    con.executemany(
        f"INSERT INTO complexes VALUES ({','.join('?' * len(COMPLEX_COLUMNS))})",
        [[c.get(k, "") for k in COMPLEX_COLUMNS] for c in complexes])
    con.execute("CREATE INDEX ix_c_design ON complexes(design_id)")
    con.execute("CREATE INDEX ix_c_target ON complexes(target)")
    con.execute("CREATE INDEX ix_d_method ON designs(design_method)")
    con.execute("""
        CREATE VIEW measured AS
        SELECT c.design_id, c.target, c.design_method, c.binding_strength,
               c.binds, c.kd_nM, c.pkd, c.platform,
               c.boltz2_ipsae, c.boltz2_min_ipsae,
               c.boltz2_iptm, c.boltz2_complex_iplddt, c.boltz2_lis,
               d.esmfold_plddt, d.proteinmpnn_score, d.seq_len, d.is_control
        FROM complexes c JOIN designs d ON d.design_id = c.design_id
        WHERE c.binding_strength != ''
    """)
    con.commit()
    con.close()
    return path


# ---------------------------------------------------------------------- QC

def qc_report(verbose: bool = True):
    """Integrity checks. Two of these found real defects in the snapshot."""
    designs, complexes = load()
    issues, notes = [], []

    # 1. Constant computational columns — a score that never varies cannot
    #    rank anything, and silently scores AUROC exactly 0.500.
    constant = {}
    for m in COMPLEX_METRICS:
        vals = [num(c, m) for c in complexes]
        vals = [v for v in vals if v is not None]
        if vals and len(set(vals)) == 1:
            constant[m] = (vals[0], len(vals))
            issues.append(f"{m} is constant ({vals[0]:g}) across all {len(vals)} rows")

    # 2. Kd must equal koff/kon. Compare only where the replicate counts agree:
    #    these columns are per-replicate geometric means, so a pair with 5 Kd
    #    but 4 kon readings disagrees by construction, not because the data is
    #    wrong. Mixing those in would manufacture a defect that isn't there.
    dis, unequal = [], 0
    for c in complexes:
        kd, kon, koff = num(c, "kd_M"), num(c, "kon_M1s1"), num(c, "koff_s1")
        if not (kd and kon and koff):
            continue
        n = {num(c, "kd_n_replicates"), num(c, "kon_n_replicates"),
             num(c, "koff_n_replicates")}
        if len(n) > 1:
            unequal += 1
            continue
        dis.append(abs(math.log10((koff / kon) / kd)))
    kin = (statistics.median(dis), max(dis), len(dis)) if dis else (0, 0, 0)
    if kin[0] > 0.01:
        issues.append(f"Kd != koff/kon (median log10 dev {kin[0]:.2f})")
    if unequal:
        notes.append(f"{unequal} pair(s) have unequal kd/kon/koff replicate "
                     f"counts — kinetics check skipped for those")

    # 3. Replicate spread = the measurement noise floor, i.e. the ceiling on
    #    any correlation a predictor can achieve against these numbers.
    spreads = [num(c, "kd_log10_spread") for c in complexes]
    spreads = [s for s in spreads if s is not None]
    noise = statistics.median(spreads) if spreads else float("nan")

    # 4. Label vs Kd consistency.
    bands, mislabelled = {}, []
    for c in complexes:
        kd_nM, lbl = num(c, "kd_nM"), c["binding_strength"]
        if kd_nM and lbl:
            bands.setdefault(lbl, []).append(kd_nM)
            if lbl == NONBINDER_LABEL:
                mislabelled.append((c["design_id"], c["target"], kd_nM))
    if mislabelled:
        notes.append(f"{len(mislabelled)} pair(s) labelled non-binder yet carry a Kd")

    # 5. A metric must not appear against two targets for one design, or the
    #    (design, target) split would be lossy.
    seen = {}
    for c in complexes:
        for m in COMPLEX_METRICS:
            if num(c, m) is not None:
                seen.setdefault((c["design_id"], m), set()).add(c["target"])
    multi = sum(1 for t in seen.values() if len(t) > 1)
    if multi:
        issues.append(f"{multi} design/metric pairs span >1 target")

    n_lab = sum(1 for c in complexes if c["binding_strength"])
    n_pos = sum(1 for c in complexes if c["binding_strength"] in BINDER_LABELS)
    n_kd = sum(1 for c in complexes if c["kd_M"])

    if verbose:
        print("=" * 68)
        print(f"PROTEINBASE DB — snapshot {SNAPSHOT}  ({LICENSE})")
        print("=" * 68)
        prov = provenance()
        if prov:
            print(f"source   : {prov.get('url')}")
            print(f"retrieved: {prov.get('retrieved_utc')}  "
                  f"sha256={str(prov.get('raw_sha256'))[:16]}…")
        print(f"designs  : {len(designs)}   complexes (design,target): {len(complexes)}")
        print(f"assayed  : {n_lab} pairs  ->  {n_pos} binders / {n_lab - n_pos} "
              f"measured NON-binders")
        print(f"Kd       : {n_kd} pairs with a measured dissociation constant")
        print(f"\nKd replicate spread (noise floor): median {noise:.3f} log10 "
              f"= {10 ** noise:.2f}x fold")
        print(f"Kd vs koff/kon: median log10 deviation {kin[0]:.4f} "
              f"(max {kin[1]:.4f}, n={kin[2]})")
        print("\nbinder label -> measured Kd:")
        for lbl in list(BINDER_LABELS) + [NONBINDER_LABEL]:
            v = sorted(bands.get(lbl, []))
            if v:
                print(f"  {lbl:7s} n={len(v):4d}  median {statistics.median(v):9.2f} nM"
                      f"   [{v[0]:.2f}, {v[-1]:.1f}]")
        if issues:
            print("\nDEFECTS FOUND (do not use these columns):")
            for i in issues:
                print(f"  ! {i}")
        for n in notes:
            print(f"  ~ {n}")
        print(f"\n{ATTRIBUTION}")

    return {
        "n_designs": len(designs), "n_complexes": len(complexes),
        "n_labelled": n_lab, "n_binders": n_pos, "n_kd": n_kd,
        "constant_columns": constant, "noise_floor_log10": noise,
        "kinetics_consistency": kin, "label_bands": bands,
        "mislabelled": mislabelled, "issues": issues,
    }


if __name__ == "__main__":
    args = set(sys.argv[1:])
    if "--rebuild" in args:
        d, c = rebuild(force_fetch="--force" in args)
        print(f"rebuilt: {len(d)} designs, {len(c)} complexes -> {DATA_DIR}")
    qc_report()
    if "--sqlite" in args:
        p = to_sqlite()
        print(f"\nsqlite written to {p} "
              f"({p.stat().st_size / 1e6:.1f} MB) — try: "
              f'sqlite3 {p} "SELECT * FROM measured LIMIT 5"')
