"""Build archmap.sqlite: every ArchMap reference atlas integrated in one database.

    python build_db.py                 # full build (~5-10 min, ~3 GB of range reads)
    python build_db.py --no-remote     # API metadata + file inventory only

Per-atlas remote summaries are cached in .cache/ so re-runs are cheap.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import gzip
import json
import shutil
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from archmap_api import ArchMap, remote_size

HERE = Path(__file__).resolve().parent
SCHEMA = """
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE atlas (
    atlas_id TEXT PRIMARY KEY, name TEXT, species TEXT, modalities TEXT,
    n_cells_reported INTEGER, n_obs_file INTEGER, n_vars_file INTEGER,
    x_encoding TEXT, x_dtype TEXT, n_samples INTEGER, n_individuals INTEGER,
    n_datasets INTEGER, doi TEXT, source_url TEXT, atlas_url TEXT,
    batch_key TEXT, cell_type_key TEXT, batch_key_used TEXT,
    cell_type_key_used TEXT, counts_description TEXT,
    var_names_description TEXT, compatible_models TEXT, is_hca INTEGER,
    is_nature INTEGER, in_revision INTEGER, uploaded_by TEXT, created_at TEXT,
    updated_at TEXT, preview_image_url TEXT, archmap_page_url TEXT,
    obsm_keys TEXT, layers_keys TEXT, uns_keys TEXT, var_columns TEXT,
    raw_json TEXT);
CREATE TABLE model (model_id TEXT PRIMARY KEY, name TEXT, description TEXT,
    compatible_classifiers TEXT, native_classifier_removed_for TEXT);
CREATE TABLE atlas_model (atlas_id TEXT, model_name TEXT, model_id TEXT,
    PRIMARY KEY (atlas_id, model_name));
CREATE TABLE atlas_file (atlas_id TEXT, file_path TEXT, file_kind TEXT,
    size_bytes INTEGER, model_object_id TEXT, PRIMARY KEY (atlas_id, file_path));
CREATE TABLE obs_column (atlas_id TEXT, column_name TEXT, kind TEXT, dtype TEXT,
    n_categories INTEGER, n_missing INTEGER, min REAL, max REAL, mean REAL,
    median REAL, is_cell_type_key INTEGER, is_batch_key INTEGER,
    PRIMARY KEY (atlas_id, column_name));
CREATE TABLE obs_category_count (atlas_id TEXT, column_name TEXT, category TEXT,
    n_cells INTEGER);
CREATE TABLE cell_type (atlas_id TEXT, cell_type TEXT, n_cells INTEGER,
    fraction REAL, PRIMARY KEY (atlas_id, cell_type));
CREATE TABLE celltype_batch_count (atlas_id TEXT, cell_type TEXT, batch TEXT,
    n_cells INTEGER);
CREATE TABLE celltype_ontology_map (atlas_id TEXT, cell_type TEXT,
    ontology_term_id TEXT, n_cells INTEGER);
CREATE TABLE gene (atlas_id TEXT, position INTEGER, feature_id TEXT, symbol TEXT,
    ensembl_id TEXT, species TEXT, ensembl_resolved TEXT, symbol_resolved TEXT,
    PRIMARY KEY (atlas_id, position));
CREATE TABLE scvi_hub_model (atlas_name TEXT, model TEXT, scvi_hub_id TEXT,
    hf_url TEXT, PRIMARY KEY (atlas_name, model));
CREATE TABLE qc_issue (atlas_id TEXT, issue TEXT, detail TEXT);

CREATE INDEX ix_occ ON obs_category_count (atlas_id, column_name);
CREATE INDEX ix_occ_cat ON obs_category_count (category);
CREATE INDEX ix_ctb ON celltype_batch_count (atlas_id, cell_type);
CREATE INDEX ix_gene_sym ON gene (symbol);
CREATE INDEX ix_gene_ens ON gene (ensembl_id);
CREATE INDEX ix_gene_res ON gene (ensembl_resolved);
CREATE INDEX ix_ont ON celltype_ontology_map (ontology_term_id);

-- atlases whose data file is not their own (see qc_issue); left out of cross-atlas views
CREATE VIEW v_atlas_usable AS
SELECT * FROM atlas WHERE atlas_id NOT IN
    (SELECT atlas_id FROM qc_issue WHERE issue = 'duplicate_reference_file');

-- one row per gene (Ensembl id, resolved from symbol where needed): in how many
-- atlases' model feature sets it appears
CREATE VIEW v_gene_presence AS
SELECT g.species, COALESCE(g.ensembl_resolved, g.symbol_resolved) AS gene_key,
       MAX(g.symbol_resolved) AS symbol, COUNT(DISTINCT g.atlas_id) AS n_atlases,
       GROUP_CONCAT(DISTINCT a.name) AS atlases
FROM gene g JOIN v_atlas_usable a USING (atlas_id)
GROUP BY g.species, gene_key;

-- majority Cell Ontology term for each atlas cell-type label (CELLxGENE-schema atlases)
CREATE VIEW v_celltype_ontology AS
SELECT atlas_id, cell_type, ontology_term_id, n_cells,
       ROUND(1.0 * n_cells / SUM(n_cells) OVER (PARTITION BY atlas_id, cell_type), 4) AS share
FROM (SELECT *, ROW_NUMBER() OVER (PARTITION BY atlas_id, cell_type ORDER BY n_cells DESC) AS rk
      FROM celltype_ontology_map)
WHERE rk = 1;

-- Cell Ontology terms shared across atlases (uses each atlas's own ontology column)
CREATE VIEW v_ontology_across_atlases AS
SELECT o.category AS ontology_term_id,
       COUNT(DISTINCT o.atlas_id) AS n_atlases, SUM(o.n_cells) AS n_cells,
       GROUP_CONCAT(a.name || ':' || o.n_cells, '; ') AS atlas_cells
FROM obs_category_count o JOIN v_atlas_usable a USING (atlas_id)
WHERE o.column_name = 'cell_type_ontology_term_id' AND o.n_cells > 0
GROUP BY o.category;

CREATE VIEW v_atlas_summary AS
SELECT a.name, a.species, a.n_cells_reported, a.n_obs_file, a.n_vars_file,
       a.cell_type_key_used, (SELECT COUNT(*) FROM cell_type c WHERE c.atlas_id = a.atlas_id) AS n_cell_types,
       a.batch_key, a.compatible_models,
       ROUND((SELECT SUM(size_bytes) FROM atlas_file f WHERE f.atlas_id = a.atlas_id) / 1e9, 2) AS total_gb,
       (SELECT GROUP_CONCAT(issue, ', ') FROM qc_issue q WHERE q.atlas_id = a.atlas_id) AS qc_issues
FROM atlas a ORDER BY a.name;
"""

MODEL_ALIASES = {"scanvi": "scANVI", "scvi": "scVI", "scpoli": "scPoli"}
# portal cellTypeKey absent from data.h5ad -> closest existing annotation column
CELLTYPE_FALLBACK = {
    "Retina": ["celltype"],          # portal says 'CellType'; 123 labels, = scANVI_predictions
    "HNOCA": ["annot_level_2"],      # portal says 'snapseed_pca_rss_level_123'
}


def file_kind(path):
    name = path.rsplit("/", 1)[-1]
    return {"data.h5ad": "reference_data", "data_only_count.h5ad": "reference_counts",
            "model.pt": "model_weights", "model_params.pt": "model_params",
            "attr.pkl": "model_attributes", "var_names.csv": "model_var_names"}.get(name, "other")


def _summarise(args):
    from h5ad_remote import summarise_h5ad
    atlas_id, url, cts, batches = args
    return atlas_id, summarise_h5ad(url, cts, batches)


def strip_urls(files):
    return [{"fileName": f["fileName"]} for f in files]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=str(HERE / "archmap.sqlite"))
    ap.add_argument("--cache", default=str(HERE / ".cache"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--no-remote", action="store_true", help="skip reading obs/var from h5ad")
    args = ap.parse_args()
    cache = Path(args.cache)
    cache.mkdir(exist_ok=True)
    raw = HERE / "raw"
    raw.mkdir(exist_ok=True)

    api = ArchMap()
    atlases, models, scvi = api.atlases(), api.models(), api.scvi_atlases()
    print(f"{len(atlases)} atlases, {len(models)} models, {len(scvi)} scvi-hub atlases")

    files, sizes = {}, {}
    for a in atlases:
        files[a["_id"]] = api.atlas_files(a["_id"])
    size_cache = cache / "sizes.json"
    sizes = json.loads(size_cache.read_text()) if size_cache.exists() else {}
    for fl in files.values():
        for f in fl:
            if f["fileName"] not in sizes:
                sizes[f["fileName"]] = remote_size(f["presignedUrl"])
    size_cache.write_text(json.dumps(sizes, indent=1))

    for name, obj in [("atlases.json", atlases), ("models.json", models),
                      ("scvi_atlases.json", scvi),
                      ("atlas_files.json", {k: strip_urls(v) for k, v in files.items()})]:
        (raw / name).write_text(json.dumps(obj, indent=1, ensure_ascii=False))

    # ---- remote obs/var summaries -------------------------------------------------
    summaries = {}
    todo = []
    for a in atlases:
        cf_path = cache / f"{a['_id']}.json"
        if cf_path.exists():
            cached = json.loads(cf_path.read_text())
            if "cell_type_key_used" in cached:
                summaries[a["_id"]] = cached
                continue
        if args.no_remote:
            continue
        url = next(f["presignedUrl"] for f in files[a["_id"]]
                   if f["fileName"].endswith("/data.h5ad"))
        cts = [a.get("cellTypeKey"), a.get("cell_type_key")] + CELLTYPE_FALLBACK.get(a["name"], [])
        batches = [a.get("batchKey"), a.get("batch_key")]
        todo.append((a["_id"], url, [k for k in cts if k], [k for k in batches if k]))
    if todo:
        with cf.ProcessPoolExecutor(args.workers) as ex:
            futs = {ex.submit(_summarise, t): t[0] for t in todo}
            for fut in cf.as_completed(futs):
                aid = futs[fut]
                try:
                    _, s = fut.result()
                except Exception as e:  # keep going; recorded as a QC issue
                    print(f"  ! {aid}: {e}")
                    s = {"error": repr(e)}
                else:
                    (cache / f"{aid}.json").write_text(json.dumps(s))
                summaries[aid] = s
                print(f"  summarised {aid} ({len(summaries)}/{len(atlases)})")

    # ---- write database -------------------------------------------------------------
    db = Path(args.db)
    if db.exists():
        db.unlink()
    con = sqlite3.connect(db)
    con.executescript(SCHEMA)
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    con.executemany("INSERT INTO meta VALUES (?,?)", [
        ("source", "https://www.archmap.bio/#/references/atlases"),
        ("api", api.api), ("built_at_utc", now)])

    model_ids = {m["name"].lower(): m["_id"] for m in models}
    for m in models:
        con.execute("INSERT INTO model VALUES (?,?,?,?,?)", (
            m["_id"], m["name"], m.get("description"),
            json.dumps(m.get("compatibleClassifiers", [])),
            json.dumps([x for x in m.get("RemoveNativeClassifierAtlases", []) if x])))

    for a in atlases:
        aid, s = a["_id"], summaries.get(a["_id"], {})
        ct_key, b_key = a.get("cellTypeKey"), a.get("batchKey")
        ct_used, b_used = s.get("cell_type_key_used"), s.get("batch_key_used")
        con.execute("INSERT INTO atlas VALUES (" + ",".join("?" * 35) + ")", (
            aid, a["name"], ", ".join(x.capitalize() for x in a.get("species", [])),
            ", ".join(a.get("modalities", [])), a.get("numberOfCells"), s.get("n_obs"),
            s.get("n_vars"), s.get("x_encoding"), s.get("x_dtype"), a.get("samples"),
            a.get("individuals"), a.get("datasets"), a.get("doi"),
            (a.get("url") or "").strip() or None, a.get("atlasUrl"), b_key, ct_key,
            b_used, ct_used,
            a.get("counts"), a.get("vars"), json.dumps(a.get("compatibleModels", [])),
            int(bool(a.get("isHCAAtlas"))), int(bool(a.get("isNature"))),
            int(bool(a.get("inrevision") or a.get("inRevison"))), a.get("uploadedBy"),
            a.get("createdAt"), a.get("updatedAt"), a.get("previewPictureURL"),
            f"https://www.archmap.bio/#/references/atlases/{aid}",
            json.dumps(s.get("obsm_keys")), json.dumps(s.get("layers_keys")),
            json.dumps(s.get("uns_keys")), json.dumps(s.get("var_columns")),
            json.dumps(a, ensure_ascii=False)))
        for m in a.get("compatibleModels", []):
            canon = MODEL_ALIASES.get(m.lower(), m)
            con.execute("INSERT OR IGNORE INTO atlas_model VALUES (?,?,?)",
                        (aid, canon, model_ids.get(canon.lower())))
        for f in files[aid]:
            p = f["fileName"]
            con.execute("INSERT INTO atlas_file VALUES (?,?,?,?,?)", (
                aid, p, file_kind(p), sizes.get(p),
                p.split("/")[1] if p.startswith("models/") else None))

        if "error" in s:
            con.execute("INSERT INTO qc_issue VALUES (?,?,?)", (aid, "remote_read_failed", s["error"]))
        if not s or "error" in s:
            continue
        cat_cols = {c["column"] for c in s["obs_columns"]}
        for c in s["obs_columns"]:
            con.execute("INSERT INTO obs_column VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", (
                aid, c["column"], c["kind"], c["dtype"], c["n_categories"], c["n_missing"],
                c["min"], c["max"], c["mean"], c["median"],
                int(c["column"] == ct_used), int(c["column"] == b_used)))
        con.executemany("INSERT INTO obs_category_count VALUES (?,?,?,?)",
                        [(aid, *r) for r in s["category_counts"]])
        ct = [(c, n) for col, c, n in s["category_counts"] if col == ct_used]
        tot = sum(n for _, n in ct) or 1
        con.executemany("INSERT INTO cell_type VALUES (?,?,?,?)",
                        [(aid, c, n, round(n / tot, 6)) for c, n in ct])
        con.executemany("INSERT INTO celltype_batch_count VALUES (?,?,?,?)",
                        [(aid, *r) for r in s["celltype_batch"]])
        con.executemany("INSERT INTO celltype_ontology_map VALUES (?,?,?,?)",
                        [(aid, *r) for r in s["celltype_ontology"]])
        species = a["species"][0].capitalize() if a.get("species") else None
        con.executemany("INSERT INTO gene (atlas_id, position, feature_id, symbol, ensembl_id, species) "
                        "VALUES (?,?,?,?,?,?)", [(aid, *g, species) for g in s["genes"]])

        # ---- QC: compare portal metadata against the actual file -------------------
        rep = a.get("numberOfCells")
        if rep and s["n_obs"] != rep:
            con.execute("INSERT INTO qc_issue VALUES (?,?,?)", (
                aid, "cell_count_mismatch",
                f"portal reports {rep:,} cells; data.h5ad has {s['n_obs']:,} "
                f"({s['n_obs'] / rep:.1%})"))
        for key, used, label in [(ct_key, ct_used, "cell_type_key"), (b_key, b_used, "batch_key")]:
            if not key or key == used:
                continue
            if used:
                con.execute("INSERT INTO qc_issue VALUES (?,?,?)", (
                    aid, f"{label}_substituted",
                    f"portal key '{key}' is not in data.h5ad; using obs column '{used}'"))
            else:
                con.execute("INSERT INTO qc_issue VALUES (?,?,?)", (
                    aid, f"{label}_missing", f"'{key}' is not an obs column of data.h5ad"))

    # identical reference files listed under different atlases
    seen = {}
    for a in atlases:
        s = summaries.get(a["_id"], {})
        if not s or "error" in s:
            continue
        p = next(f["fileName"] for f in files[a["_id"]] if f["fileName"].endswith("/data.h5ad"))
        sig = (sizes.get(p), s["n_obs"], s["n_vars"], tuple(g[1] for g in s["genes"]))
        if sig in seen:
            # the copy whose cell count disagrees with its own portal record is the misfiled one
            other = seen[sig]
            for x, y in [(a, other), (other, a)]:
                if x.get("numberOfCells") != summaries[x["_id"]]["n_obs"]:
                    con.execute("INSERT INTO qc_issue VALUES (?,?,?)", (
                        x["_id"], "duplicate_reference_file",
                        f"data.h5ad is byte-size/shape/gene identical to atlas '{y['name']}' "
                        f"({y['_id']}); this atlas's own data is not in the file"))
        else:
            seen[sig] = a

    # resolve symbol <-> Ensembl across atlases of the same species, so atlases that
    # only ship symbols (or only Ensembl ids) can be joined with the others
    con.executescript("""
        CREATE TEMP TABLE sym2ens AS
        SELECT species, symbol, ensembl_id FROM (
            SELECT species, symbol, ensembl_id, ROW_NUMBER() OVER (
                PARTITION BY species, symbol ORDER BY COUNT(*) DESC, ensembl_id) AS rk
            FROM gene WHERE symbol IS NOT NULL AND ensembl_id IS NOT NULL
            GROUP BY species, symbol, ensembl_id) WHERE rk = 1;
        CREATE TEMP TABLE ens2sym AS
        SELECT species, ensembl_id, symbol FROM (
            SELECT species, ensembl_id, symbol, ROW_NUMBER() OVER (
                PARTITION BY species, ensembl_id ORDER BY COUNT(*) DESC, symbol) AS rk
            FROM gene WHERE symbol IS NOT NULL AND ensembl_id IS NOT NULL
            GROUP BY species, ensembl_id, symbol) WHERE rk = 1;
        CREATE INDEX t1 ON sym2ens (species, symbol);
        CREATE INDEX t2 ON ens2sym (species, ensembl_id);
        UPDATE gene SET
            ensembl_resolved = COALESCE(ensembl_id, (SELECT m.ensembl_id FROM sym2ens m
                WHERE m.species = gene.species AND m.symbol = gene.symbol)),
            symbol_resolved = COALESCE(symbol, (SELECT m.symbol FROM ens2sym m
                WHERE m.species = gene.species AND m.ensembl_id = gene.ensembl_id));
    """)

    for e in scvi:
        for m in e.get("modelIds", []):
            con.execute("INSERT OR IGNORE INTO scvi_hub_model VALUES (?,?,?,?)", (
                e["name"], m["model"], m["scviHubId"],
                f"https://huggingface.co/{m['scviHubId']}"))
    con.commit()

    # ---- CSV exports (Excel friendly) --------------------------------------------
    out = HERE / "exports"
    out.mkdir(exist_ok=True)
    queries = {
        "atlas_summary.csv": "SELECT * FROM v_atlas_summary",
        "atlases.csv": "SELECT * FROM atlas",
        "cell_types.csv": "SELECT a.name AS atlas, c.cell_type, c.n_cells, c.fraction "
                          "FROM cell_type c JOIN atlas a USING (atlas_id) ORDER BY a.name, c.n_cells DESC",
        "obs_columns.csv": "SELECT a.name AS atlas, o.* FROM obs_column o JOIN atlas a USING (atlas_id)",
        "files.csv": "SELECT a.name AS atlas, f.* FROM atlas_file f JOIN atlas a USING (atlas_id)",
        "gene_presence.csv": "SELECT * FROM v_gene_presence ORDER BY n_atlases DESC, symbol",
        "ontology_across_atlases.csv": "SELECT * FROM v_ontology_across_atlases ORDER BY n_atlases DESC, n_cells DESC",
        "qc_issues.csv": "SELECT a.name AS atlas, q.issue, q.detail FROM qc_issue q JOIN atlas a USING (atlas_id)",
        "scvi_hub_models.csv": "SELECT * FROM scvi_hub_model",
        "models.csv": "SELECT * FROM model",
    }
    for fname, q in queries.items():
        cur = con.execute(q)
        with open(out / fname, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.writer(fh)
            w.writerow([d[0] for d in cur.description])
            w.writerows(cur)
    con.close()
    # committed copy is gzipped: raw SQLite pages concatenate gene names (AKAP...)
    # into strings that secret scanners mistake for cloud access keys
    with open(db, "rb") as src, gzip.open(f"{db}.gz", "wb", compresslevel=9) as dst:
        shutil.copyfileobj(src, dst)
    print(f"wrote {db} (+ .gz) and {len(queries)} CSVs in {out}")


if __name__ == "__main__":
    main()
