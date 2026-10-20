"""Build cell state library v0.2.

    python -m cell_state_library.build [--source PATH]

Downloads the author matrix at a pinned commit (or uses --source), verifies its
sha256, builds data/cell_state_v0_2.sqlite, runs identifiability checks, selects
gene-pair panels and writes manifest / validation / identifiability reports.
The source CSV and the SQLite file stay out of git (see README: no explicit
license on the author repository).
"""
import argparse
import collections
import csv
import hashlib
import json
import math
import sqlite3
import urllib.request
from pathlib import Path

from . import VERSION, identifiability as ident, panel_select
from .nnls import score_nnls

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
DB = DATA / "cell_state_v0_2.sqlite"
REPORTS = ROOT / "reports"

SOURCE_ID = "HIRISA_GITHUB"
REPO = "https://github.com/aifimmunology/IFN_manuscript"
COMMIT = "3f4b18740c96adc0643be2a64c73f484eea8d270"
CSV_PATH = "Inputs/IFNa_IFNg_Scoring_ISGs_Mat.csv"
RAW_URL = f"https://raw.githubusercontent.com/aifimmunology/IFN_manuscript/{COMMIT}/{CSV_PATH}"
SHA256 = "813fa4aa0b69686a9cfeff1c62e0eca7e2b21df02460be55a4b07cd29c9ba584"
METHODS_URL = "https://apps.allenimmunology.org/aifi/resources/ifn-response/methods/"
ANALYSIS_URL = "https://apps.allenimmunology.org/aifi/resources/ifn-response/analysis/"

PANEL_K = 8                                  # pairs per panel -> 16 genes
SIGMA = ident.SIGMA_REF
TAU = SIGMA * math.sqrt(2)                   # dead zone = 1 SD of a pair difference under noise

# cell_type_id: (display, level, parent, culture context when stimulated)
CELL_TYPES = {
    "PBMC":       ("PBMC", "bulk", None, "total PBMC"),
    "Monocyte":   ("Monocyte (CD14-enriched)", "L1", "PBMC", "enriched CD14 monocytes (magnetic negative selection)"),
    "Tcell":      ("T cell", "L1", "PBMC", "enriched T cells (magnetic negative selection)"),
    "Bcell":      ("B cell", "L1", "PBMC", "enriched B cells (magnetic negative selection)"),
    "NK":         ("NK cell", "L1", "PBMC", "enriched NK cells (magnetic negative selection)"),
    "CD4_Naive":  ("CD4 naive T", "L2", "Tcell", None), "CD4_Memory": ("CD4 memory T", "L2", "Tcell", None),
    "CD8_Naive":  ("CD8 naive T", "L2", "Tcell", None), "CD8_Memory": ("CD8 memory T", "L2", "Tcell", None),
    "Treg":       ("Treg", "L2", "Tcell", None), "MAIT": ("MAIT", "L2", "Tcell", None), "gdT": ("gamma-delta T", "L2", "Tcell", None),
    "B_Naive":    ("B naive", "L2", "Bcell", None), "B_Memory": ("B memory", "L2", "Bcell", None),
    "Plasma":     ("Plasma cell", "L2", "Bcell", None),
    "NK.CD56hi":  ("NK CD56bright", "L2", "NK", None), "NK.CD56dim": ("NK CD56dim", "L2", "NK", None),
}
# The 13 subsets activated in v0.2: the 12 L2 labels accepted by the author R code + Monocyte.
ACTIVE = ["Monocyte", "CD4_Naive", "CD4_Memory", "CD8_Naive", "CD8_Memory", "Treg", "MAIT", "gdT",
          "B_Naive", "B_Memory", "Plasma", "NK.CD56hi", "NK.CD56dim"]

CONDITIONS = [
    ("CTRL", "none", None, None, 21.0, "source_protocol_reference"),
    ("IFNA", "IFN-alpha2a", 1000.0, "U/mL", 21.0, "source_protocol_reference"),
    ("IFNG", "IFN-gamma", 25.0, "ng/mL", 21.0, "source_protocol_reference"),
    ("IFNB", "IFN-beta", 100.0, "U/mL", 21.0, "protocol_only_no_reference_in_matrix"),
    ("IFNL1", "IFN-lambda1", 100.0, "ng/mL", 21.0, "protocol_only_no_reference_in_matrix"),
]
SUMMARY_STAT = ("median across 5 donors of per-donor MAST log2FC (1,000-cell subsamples, 20 iterations) "
                "vs matched unstimulated 21 h control")
GENE_SELECTION = ("author selection: padj<0.01 and |log2FC|>0.2 in >=10/20 iterations, shared by all 5 donors; "
                  "selected per cell type, so gene lists differ between cell types")


def fetch_source(path):
    if path is None:
        path = DATA / "source" / "IFNa_IFNg_Scoring_ISGs_Mat.csv"
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(RAW_URL, timeout=60) as r:
                path.write_bytes(r.read())
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if digest != SHA256:
        raise SystemExit(f"sha256 mismatch for {path}: {digest} != pinned {SHA256}")
    return Path(path), digest


def load_groups(path):
    records = list(csv.DictReader(open(path, newline="")))
    keys = [(r["celltype_level"], r["subset"], r["Gene"]) for r in records]
    assert len(set(keys)) == len(keys), "duplicate (level, subset, gene) rows"
    groups = collections.defaultdict(dict)
    for r in records:
        a, g = float(r["IFNa"]), float(r["IFNg"])
        assert math.isfinite(a) and math.isfinite(g)
        groups[(r["celltype_level"], r["subset"])][r["Gene"]] = (a, g)
    assert groups[("L1", "Monocyte")] == groups[("L2", "Monocyte")], "L1/L2 Monocyte expected identical"
    return records, groups


def alias_target(level, label):
    return "PBMC" if level == "PBMC" else label


def culture_context(ct):
    c = ct
    while CELL_TYPES[c][3] is None:
        c = CELL_TYPES[c][2]
    ctx = CELL_TYPES[c][3]
    if c != ct:
        ctx += f"; '{ct}' assigned after scRNA-seq clustering"
        if ct == "Plasma":
            ctx += " (B-cell culture assumed from lineage; verify)"
    return ctx


def build(source=None):
    path, digest = fetch_source(source)
    records, groups = load_groups(path)
    DATA.mkdir(exist_ok=True)
    REPORTS.mkdir(exist_ok=True)
    if DB.exists():
        DB.unlink()
    con = sqlite3.connect(DB)
    con.executescript((ROOT / "schema.sql").read_text())

    con.execute("INSERT INTO sources VALUES(?,?,?,?,?,?,?)", (
        SOURCE_ID, f"{REPO}/blob/{COMMIT}/{CSV_PATH}", COMMIT, digest, "published_aggregate_response_matrix",
        "downloaded_public_github",
        "No explicit license found in the author repository; keep attribution and confirm terms before commercial redistribution."))
    for ct, (name, level, parent, _) in CELL_TYPES.items():
        note = "duplicated in source as L1 and L2 with identical values; stored once" if ct == "Monocyte" else None
        con.execute("INSERT INTO cell_types VALUES(?,?,?,?,?)", (ct, name, level, parent, note))
    for (level, label) in groups:
        con.execute("INSERT INTO cell_type_aliases VALUES(?,?,?,?)", (SOURCE_ID, level, label, alias_target(level, label)))
    con.executemany("INSERT INTO conditions VALUES(?,?,?,?,?,?,?)",
                    [(c, s, d, u, t, METHODS_URL, st) for c, s, d, u, t, st in CONDITIONS])

    # profiles + long-format values
    loaded = {}
    for (level, label), genes in sorted(groups.items()):
        ct = alias_target(level, label)
        if ct in loaded:
            continue  # L2 Monocyte duplicate of L1 Monocyte
        loaded[ct] = genes
        status = "active" if ct in ACTIVE else "staged_parent_level"
        for cond, idx in (("IFNA", 0), ("IFNG", 1)):
            pid = f"{SOURCE_ID}:{ct}:{cond}"
            con.execute("INSERT INTO reference_profiles VALUES(?,?,?,?,?,?,?,?,?,?)", (
                pid, SOURCE_ID, ct, cond, "CTRL", culture_context(ct), 5, SUMMARY_STAT, GENE_SELECTION, status))
            con.executemany("INSERT INTO reference_values VALUES(?,?,?,NULL)",
                            [(pid, gene, v[idx]) for gene, v in sorted(genes.items())])

    # identifiability for the 13 active subsets
    ident_rows = {}
    for ct in ACTIVE:
        genes = sorted(loaded[ct])
        a = [loaded[ct][x][0] for x in genes]
        g = [loaded[ct][x][1] for x in genes]
        r = ident.check(a, g)
        sens = {f"sigma_{s}": ident.check(a, g, sigma=s)["status"] for s in (0.25, 0.5, 1.0)}
        author_rule = int(ct in ident.AUTHOR_IFNG_NONRESPONDING)
        notes = []
        if ct == "CD8_Memory":
            notes.append("author R list says 'Memory_CD8', which never matches this label, so the author zero rule does not fire")
        if author_rule and r["status"] == "both":
            notes.append("author forces IFNg_score=0 here, but the data pass the identifiability gate")
        if not author_rule and r["status"] == "ifna_only" and ct != "CD8_Memory":
            notes.append("author does not zero IFNg here, but IFNg fails the identifiability gate")
        r["sensitivity"] = sens
        r["author_ifng_zero_rule"] = author_rule
        r["notes"] = "; ".join(notes) or None
        ident_rows[ct] = r
        con.execute("INSERT INTO identifiability VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
            ct, SIGMA, r["n_genes"], r["rms_ifna"], r["rms_ifng"], r["cosine"], r["condition_number"],
            r["se_ifna"], r["se_ifng"], r["z_ifna"], r["z_ifng"], r["status"], author_rule, r["notes"]))
        ok = ident.identifiable_stimuli(r["status"])
        for cond in ("IFNA", "IFNG"):
            if cond not in ok:
                con.execute("UPDATE reference_profiles SET status='active_reference_not_identifiable' WHERE profile_id=?",
                            (f"{SOURCE_ID}:{ct}:{cond}",))

    # gene-pair panels
    panels = {}
    for ct in ACTIVE:
        genes = sorted(loaded[ct])
        classes = ident.identifiable_stimuli(ident_rows[ct]["status"])
        panel = panel_select.select_panel(genes, [loaded[ct][x][0] for x in genes], [loaded[ct][x][1] for x in genes],
                                          classes, PANEL_K, SIGMA, TAU)
        panels[ct] = panel
        con.executemany("INSERT INTO pair_panels VALUES(?,?,?,?,?,?,?,?,?)", [
            (f"PAIR_v0_2:{ct}", ct, p["rank"], p["gene_i"], p["gene_j"], p["d_ifna"], p["d_ifng"], TAU, SIGMA) for p in panel])

    # v0.1 planned experiments carried over (still not performed)
    plans = []
    for i in range(1, 6):
        for c in ("CTRL", "IFNA", "IFNG"):
            plans.append((f"P{i:02d}_{c}", f"UNASSIGNED_DONOR_SLOT_{i}", "Monocyte", c,
                          None if c == "CTRL" else f"P{i:02d}_CTRL", "planned_not_performed"))
    con.executemany("INSERT INTO planned_experiments VALUES(?,?,?,?,?,?)", plans)
    con.commit()
    assert con.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
    assert not con.execute("PRAGMA foreign_key_check").fetchall()
    assert con.execute("SELECT COUNT(*) FROM samples").fetchone()[0] == 0
    n_values = con.execute("SELECT COUNT(*) FROM reference_values").fetchone()[0]
    expected = 2 * sum(len(v) for v in loaded.values())
    assert n_values == expected
    con.close()

    validation = self_checks(loaded, ident_rows)
    write_reports(records, groups, loaded, ident_rows, panels, digest, n_values, validation)
    return loaded, ident_rows, panels


def self_checks(loaded, ident_rows):
    """Software checks only: reconstruction and error handling, not biological accuracy."""
    out = {"scope": "Software integrity and self-reconstruction only. Not independent biological accuracy, "
                    "not leave-donor-out validation."}
    # v0.1 expected values on Monocyte (from v0.1 validation.json)
    mono = loaded["Monocyte"]
    refs = {"IFNA": {k: v[0] for k, v in mono.items()}, "IFNG": {k: v[1] for k, v in mono.items()}}
    for name, key in (("IFNA", 0), ("IFNG", 1)):
        r = score_nnls(refs, {k: v[key] for k, v in mono.items()})
        assert r["relative_residual_norm"] < 1e-8
        assert abs(r["scores"]["IFNA"] - (name == "IFNA")) < 1e-8 and abs(r["scores"]["IFNG"] - (name == "IFNG")) < 1e-8
        assert abs(r["reference_cosines"]["IFNA~IFNG"] - 0.6771438363501308) < 1e-12
        assert abs(r["reference_condition_number"] - 2.2791904382226296) < 1e-9
    out["v0_1_monocyte_compatibility"] = "passed (scores, cosine 0.67714, condition number 2.27919 match v0.1)"
    recon = {}
    for ct in ACTIVE:
        g = loaded[ct]
        refs = {"IFNA": {k: v[0] for k, v in g.items()}, "IFNG": {k: v[1] for k, v in g.items()}}
        res = {}
        for name, key in (("IFNA", 0), ("IFNG", 1)):
            r = score_nnls(refs, {k: v[key] for k, v in g.items()})
            res[name] = round(r["scores"][name], 6)
            assert r["relative_residual_norm"] < 1e-6
        recon[ct] = res
    out["self_reconstruction_13_subsets"] = recon
    errors = {}
    for name, refs_, target in (
            ("missing_genes", {"IFNA": {k: v[0] for k, v in mono.items()}, "IFNG": {k: v[1] for k, v in mono.items()}}, {"ISG15": 1.0}),
            ("zero_signal", {"IFNA": {k: v[0] for k, v in mono.items()}, "IFNG": {k: v[1] for k, v in mono.items()}}, {k: 0.0 for k in mono}),
            ("non_finite", {"IFNA": {k: v[0] for k, v in mono.items()}, "IFNG": {k: v[1] for k, v in mono.items()}}, {"ISG15": float("nan")})):
        try:
            score_nnls(refs_, target)
        except ValueError as e:
            errors[name] = str(e)
        else:
            raise AssertionError(name)
    out["rejected_inputs"] = errors
    return out


def write_reports(records, groups, loaded, ident_rows, panels, digest, n_values, validation):
    manifest = {
        "version": VERSION, "created": "2026-10-20",
        "schema": "long format: reference_profiles (cell type x condition vs control) + reference_values (profile x gene)",
        "active_cell_types": ACTIVE,
        "staged_cell_types": sorted(set(loaded) - set(ACTIVE)),
        "identifiability_summary": {ct: r["status"] for ct, r in ident_rows.items()},
        "reference_values_rows": n_values,
        "pair_panel": {"pairs_per_panel": PANEL_K, "genes_per_panel": 2 * PANEL_K, "sigma": SIGMA, "tau": TAU},
        "own_measured_samples": 0, "planned_samples": 15,
        "source_rows": len(records), "source_label_groups": len(groups),
        "source_url": REPO, "source_commit": COMMIT, "source_sha256": digest,
        "methods_url": METHODS_URL, "analysis_url": ANALYSIS_URL,
        "license_note": "No explicit license found in the author repository. Source CSV and built SQLite are not committed; "
                        "build.py downloads the pinned file. Confirm reuse rights before commercial redistribution.",
        "status": "13 subsets active as research references; no own experiments; no patient validation",
    }
    (REPORTS / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    (REPORTS / "validation.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2))
    (REPORTS / "identifiability.json").write_text(json.dumps(
        {"assumed_sigma": SIGMA, "z_min": ident.Z_MIN, "min_genes": ident.MIN_GENES, "cell_types": ident_rows},
        ensure_ascii=False, indent=2))
    (REPORTS / "pair_panels.json").write_text(json.dumps(
        {"sigma": SIGMA, "tau": TAU, "panels": panels}, ensure_ascii=False, indent=2))
    (ROOT / "own_sample_template.json").write_text(json.dumps({
        "sample_id": "", "donor_id": "", "cell_type_id": "Monocyte", "condition_id": "", "paired_control_id": "",
        "protocol_version": "", "culture_context": "enriched_subset_or_PBMC_specify", "medium": "", "serum_type": "",
        "serum_lot": "", "stimulus": "", "stimulus_dose": None, "stimulus_unit": "", "stimulus_lot": "",
        "duration_h": None, "temperature_C": None, "cell_density": None, "cell_density_unit": "", "purity_pct": None,
        "viability_pct": None, "assay": "", "batch_id": "", "replicate_kind": "biological", "data_uri": "",
        "qc_status": "not_measured", "functional_readouts": []}, ensure_ascii=False, indent=2))
    (ROOT / "expansion_catalog.json").write_text(json.dumps([
        {"step": "done_in_v0_2", "cells": "12 author L2 subsets + Monocyte", "resource": "HIRISA author matrix", "url": REPO,
         "data_status": "activated with per-subset identifiability; IFNg withheld in 5 subsets"},
        {"step": 1, "cells": "same 13 subsets, donor level", "resource": "HIRISA per-donor DEG tables (HISE, login)",
         "url": ANALYSIS_URL, "data_status": "not accessible without sign-in",
         "entry_gate": "replaces the assumed sigma with measured donor noise; enables leave-one-donor-out validation"},
        {"step": 2, "cells": "Monocyte, T, B, NK in PBMC; 90 cytokines", "resource": "Parse 10M PBMC cytokines",
         "url": "https://apps.allenimmunology.org/aifi/resources/parse-10m-cytokines/",
         "data_status": "metadata verified, expression not downloaded",
         "entry_gate": "new profiles in the same long tables (condition per cytokine); 24 h PBMC context differs from 21 h enriched cultures; CC BY-NC 4.0"},
        {"step": 3, "cells": "Mouse immune cells, 86 cytokines in vivo", "resource": "Immune Dictionary",
         "url": "https://www.immune-dictionary.org/", "data_status": "not imported",
         "entry_gate": "species and in vivo context recorded on the profile; cross-species use needs validation"},
        {"step": 4, "cells": "Cancer cell lines x drugs", "resource": "Tahoe-100M / LINCS L1000",
         "url": "https://huggingface.co/datasets/tahoebio/Tahoe-100M", "data_status": "not imported",
         "entry_gate": "separate source_id and evidence type; match line, compound, dose, time, control"},
    ], ensure_ascii=False, indent=2))
    (ROOT / "target_query_template.json").write_text(json.dumps({
        "cell_type": "CD4_Naive", "celltype_level": "L2",
        "contrast": "REPLACE: target vs matched control, same cell type",
        "log2fc": {}}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, help="local copy of IFNa_IFNg_Scoring_ISGs_Mat.csv (sha256 is checked)")
    args = ap.parse_args()
    loaded, ident_rows, panels = build(args.source)
    print(f"built {DB.relative_to(ROOT.parent)}")
    for ct, r in ident_rows.items():
        print(f"  {ct:11s} n={r['n_genes']:4d} cos={r['cosine']:.2f} z_a={r['z_ifna']:6.1f} z_g={r['z_ifng']:6.1f} "
              f"-> {r['status']:9s} author_zero={r['author_ifng_zero_rule']}  panel={len(panels[ct])} pairs")
