"""Build db/fda_ivd_markers.sqlite from the FDA marker catalog workbook.

Layers (each table carries a `source` so the provenance is explicit):
  marker, marker_submission          <- catalog workbook (CLIA-derived, catalogue level)
  denovo_catalog, pma_catalog, cdx, nat_catalog <- catalog supplementary sheets
  submission                          <- union of every K/DEN/P number seen; enriched later by openfda_enrich.py
  product_code                        <- openFDA classification (openfda_enrich.py)
  document, extraction                <- decision-summary fetch + parse (fetch_ds.py, ds_parse.py)
  curated_target, curated_section     <- the 25-analyte hand-curated atlas (load_curated.py)

Usage: python3 db/build_db.py <catalog.xlsx> [db path]
"""
import sys, os, re, sqlite3, datetime
import openpyxl

XLSX = sys.argv[1]
DB = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fda_ivd_markers.sqlite')

SCHEMA = """
PRAGMA journal_mode=WAL;
CREATE TABLE IF NOT EXISTS marker (
  marker_id TEXT PRIMARY KEY, category_cn TEXT, label TEXT, tier TEXT,
  n_510k INTEGER, n_denovo INTEGER, n_pma INTEGER, link_status TEXT,
  ex_510k TEXT, ex_denovo TEXT, ex_pma TEXT, clia_analyte_id TEXT, n_clia INTEGER, clia_doc_examples TEXT,
  notes TEXT, official_link TEXT, source TEXT DEFAULT 'catalog_v1');
CREATE TABLE IF NOT EXISTS marker_submission (
  marker_id TEXT, marker_label TEXT, pathway TEXT, submission_no TEXT, product_name TEXT,
  qualifier1 TEXT, qualifier2 TEXT, fda_link TEXT, source TEXT DEFAULT 'catalog_v1',
  PRIMARY KEY (marker_id, pathway, submission_no, product_name, qualifier1, qualifier2));
CREATE TABLE IF NOT EXISTS denovo_catalog (den_no TEXT PRIMARY KEY, product_name TEXT, product_code TEXT, classification_name TEXT,
  classification_def TEXT, specialty TEXT, decision_date TEXT, tier_note TEXT, fda_link TEXT);
CREATE TABLE IF NOT EXISTS pma_catalog (pma_no TEXT PRIMARY KEY, product_name TEXT, generic_name TEXT, product_code TEXT, classification_def TEXT,
  decision_date TEXT, decision_code TEXT, n_records INTEGER, status_note TEXT, approval_order_excerpt TEXT, fda_link TEXT);
CREATE TABLE IF NOT EXISTS cdx (record_id TEXT PRIMARY KEY, biomarker TEXT, variant_detail TEXT, device TEXT, indication_sample TEXT,
  submission_no TEXT, pathway TEXT, ref_no_date TEXT, scope_note TEXT, fda_link TEXT);
CREATE TABLE IF NOT EXISTS nat_catalog (record_id TEXT PRIMARY KEY, use_or_organism TEXT, product_name TEXT, manufacturer TEXT,
  submission_no TEXT, pathway TEXT, status_note TEXT, fda_link TEXT);
CREATE TABLE IF NOT EXISTS category_overview (category_cn TEXT PRIMARY KEY, n_labels INTEGER, n_510k INTEGER, n_denovo INTEGER, n_pma INTEGER, n_pending INTEGER, examples TEXT);
CREATE TABLE IF NOT EXISTS method_notes (item TEXT, handling TEXT, official_source TEXT);
CREATE TABLE IF NOT EXISTS submission (
  submission_no TEXT PRIMARY KEY, pathway TEXT, year INTEGER,
  decision_date TEXT, decision TEXT, applicant TEXT, device_name TEXT, product_code TEXT, regulation_number TEXT,
  advisory_committee TEXT, statement_or_summary TEXT, clearance_type TEXT, third_party TEXT, expedited TEXT,
  in_catalog INTEGER DEFAULT 0, openfda_found INTEGER DEFAULT 0, has_decision_summary INTEGER, has_510k_summary INTEGER, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS product_code (product_code TEXT PRIMARY KEY, device_name TEXT, regulation_number TEXT, device_class TEXT,
  medical_specialty TEXT, review_panel TEXT, n_510k_total INTEGER, definition TEXT);
CREATE TABLE IF NOT EXISTS document (submission_no TEXT, doc_type TEXT, url TEXT, n_chars INTEGER, txt_path TEXT, fetched_at TEXT,
  PRIMARY KEY (submission_no, doc_type));
CREATE TABLE IF NOT EXISTS extraction (submission_no TEXT, field TEXT, value TEXT, section TEXT, template TEXT, confidence TEXT,
  PRIMARY KEY (submission_no, field));
CREATE TABLE IF NOT EXISTS curated_target (target_id TEXT PRIMARY KEY, group_letter TEXT, group_name TEXT, title TEXT, product_codes TEXT,
  archetype TEXT, representative_cutoff TEXT, marker_ids TEXT, source TEXT DEFAULT 'atlas_2026-09-22');
CREATE TABLE IF NOT EXISTS curated_section (target_id TEXT, section_no INTEGER, section_title TEXT, body_md TEXT, PRIMARY KEY (target_id, section_no));
CREATE INDEX IF NOT EXISTS ix_ms_sub ON marker_submission(submission_no);
CREATE INDEX IF NOT EXISTS ix_ms_marker ON marker_submission(marker_id);
CREATE INDEX IF NOT EXISTS ix_sub_pc ON submission(product_code);
CREATE INDEX IF NOT EXISTS ix_sub_year ON submission(year);
CREATE INDEX IF NOT EXISTS ix_ext_field ON extraction(field);
"""

def year_of(no):
    m = re.match(r'^K(\d\d)\d{4}$', no or '')
    if m:
        y = int(m.group(1)); return 1900 + y if y > 50 else 2000 + y
    m = re.match(r'^DEN(\d\d)', no or '')
    if m:
        return 2000 + int(m.group(1))
    m = re.match(r'^P(\d\d)', no or '')
    if m:
        y = int(m.group(1)); return 1900 + y if y > 50 else 2000 + y
    return None

def s(v):
    if v is None: return None
    if isinstance(v, datetime.datetime): return v.date().isoformat()
    return str(v).strip()

def main():
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    def rows(name, hdr=4):
        r = list(wb[name].iter_rows(values_only=True)); return [x for x in r[hdr+1:] if x and x[0] is not None]
    if os.path.exists(DB): os.remove(DB)
    con = sqlite3.connect(DB); con.executescript(SCHEMA); cur = con.cursor()

    for r in rows('分类概览'):
        cur.execute('INSERT OR REPLACE INTO category_overview VALUES (?,?,?,?,?,?,?)', [s(x) for x in r[:7]])
    for r in rows('方法与来源'):
        cur.execute('INSERT INTO method_notes VALUES (?,?,?)', [s(x) for x in (list(r)+[None]*3)[:3]])
    for r in rows('Marker清单'):
        r = (list(r) + [None]*16)[:16]
        cur.execute('INSERT OR REPLACE INTO marker (marker_id,category_cn,label,tier,n_510k,n_denovo,n_pma,link_status,ex_510k,ex_denovo,ex_pma,clia_analyte_id,n_clia,clia_doc_examples,notes,official_link) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                    [s(r[0]), s(r[1]), s(r[2]), s(r[3]), int(r[4] or 0), int(r[5] or 0), int(r[6] or 0), s(r[7]), s(r[8]), s(r[9]), s(r[10]), s(r[11]), int(r[12] or 0), s(r[13]), s(r[14]), s(r[15])])
    subs = {}
    for r in rows('申报关联'):
        r = (list(r) + [None]*8)[:8]
        cur.execute('INSERT OR IGNORE INTO marker_submission (marker_id,marker_label,pathway,submission_no,product_name,qualifier1,qualifier2,fda_link) VALUES (?,?,?,?,?,?,?,?)',
                    [s(x) or '' for x in r])
        no = s(r[3]) or ''
        base = no.split('/')[0]
        if base: subs.setdefault(base, s(r[2]))
    for r in rows('DeNovo补充'):
        r = (list(r) + [None]*9)[:9]
        cur.execute('INSERT OR REPLACE INTO denovo_catalog VALUES (?,?,?,?,?,?,?,?,?)', [s(x) for x in r])
        subs.setdefault(s(r[0]), 'De Novo')
    for r in rows('PMA补充'):
        r = (list(r) + [None]*11)[:11]
        cur.execute('INSERT OR REPLACE INTO pma_catalog VALUES (?,?,?,?,?,?,?,?,?,?,?)', [s(r[0]), s(r[1]), s(r[2]), s(r[3]), s(r[4]), s(r[5]), s(r[6]), int(r[7] or 0), s(r[8]), s(r[9]), s(r[10])])
        subs.setdefault(s(r[0]), 'PMA')
    for r in rows('CDx标志物补充'):
        r = (list(r) + [None]*10)[:10]
        cur.execute('INSERT OR REPLACE INTO cdx VALUES (?,?,?,?,?,?,?,?,?,?)', [s(x) for x in r])
    for r in rows('核酸检测补充'):
        r = (list(r) + [None]*8)[:8]
        cur.execute('INSERT OR REPLACE INTO nat_catalog VALUES (?,?,?,?,?,?,?,?)', [s(x) for x in r])
        no = (s(r[4]) or '').split('/')[0]
        if re.match(r'^(K|DEN|P)\d+', no): subs.setdefault(no, s(r[5]))
    for no, path in subs.items():
        if not re.match(r'^(K|DEN|P)\d+', no): continue
        pw = 'De Novo' if no.startswith('DEN') else ('PMA' if no.startswith('P') else '510(k)')
        cur.execute('INSERT OR IGNORE INTO submission (submission_no, pathway, year, in_catalog) VALUES (?,?,?,1)', (no, pw, year_of(no)))
    con.commit()
    for t in ['marker','marker_submission','submission','denovo_catalog','pma_catalog','cdx','nat_catalog']:
        print(t, cur.execute(f'SELECT COUNT(*) FROM {t}').fetchone()[0])
    con.close()

if __name__ == '__main__':
    main()
