"""Export the database to CSV (one file per table) and to a compact JSON bundle used by the HTML browser.

Usage: python3 db/export.py [db path] [out dir]
"""
import sys, os, csv, json, sqlite3, re
HERE = os.path.dirname(os.path.abspath(__file__))
DB = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'fda_ivd_markers.sqlite')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'export')
os.makedirs(OUT, exist_ok=True)
# The browser bundle carries a clipped subset of fields (full text stays in SQLite / CSV); sizes are tuned to keep the single-file page well under 16 MB.
KEY_FIELDS = ['measurand', 'type_of_test', 'intended_use', 'specimen_types', 'assay_cutoff', 'clinical_cutoff', 'cutoff_numbers',
              'reference_range', 'clsi_codes', 'detection_limit', 'mc_slope', 'mc_r', 'clinical_studies', 'sens_pct', 'spec_pct', 'sample_n', 'predicate']
SHORT = {'intended_use': 320, 'assay_cutoff': 320, 'clinical_cutoff': 320, 'reference_range': 320, 'detection_limit': 160, 'clinical_studies': 320,
         'predicate': 120, 'measurand': 80, 'type_of_test': 80, 'clsi_codes': 160, 'cutoff_numbers': 160, 'sens_pct': 80, 'spec_pct': 80, 'sample_n': 60, 'mc_slope': 60, 'mc_r': 60, 'specimen_types': 120}

def clip(s, n):
    s = re.sub(r'\s+', ' ', s or '').strip()
    return s if len(s) <= n else s[:n - 1] + '…'

con = sqlite3.connect(DB); con.row_factory = sqlite3.Row; cur = con.cursor()
tables = [r[0] for r in cur.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()]
for t in tables:
    rows = cur.execute(f'SELECT * FROM {t}').fetchall()
    with open(os.path.join(OUT, f'{t}.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        if rows:
            w.writerow(rows[0].keys())
            for r in rows: w.writerow([r[k] for k in r.keys()])
    print(f'{t}: {len(rows)} rows')

# JSON bundle for the browser
markers = [dict(r) for r in cur.execute('SELECT marker_id, category_cn, label, tier, n_510k, n_denovo, n_pma, link_status, clia_analyte_id FROM marker ORDER BY category_cn, label')]
links = {}
for r in cur.execute('SELECT marker_id, pathway, submission_no FROM marker_submission'):
    links.setdefault(r['marker_id'], []).append([r['pathway'], r['submission_no']])
subs = {}
for r in cur.execute('SELECT * FROM submission'):
    subs[r['submission_no']] = {'p': r['pathway'], 'y': r['year'], 'd': r['decision_date'], 'a': clip(r['applicant'], 60), 'n': clip(r['device_name'], 110), 'pc': r['product_code'],
                                'reg': r['regulation_number'], 'cat': r['in_catalog'], 'ds': r['has_decision_summary'], 'sm': r['has_510k_summary']}
ext = {}
for r in cur.execute('SELECT submission_no, field, value, template FROM extraction'):
    if r['field'] in KEY_FIELDS:
        ext.setdefault(r['submission_no'], {'tpl': r['template']})[r['field']] = clip(r['value'], SHORT.get(r['field'], 300))
pcs = {r['product_code']: dict(r) for r in cur.execute('SELECT * FROM product_code')}
curated = []
for r in con.execute('SELECT * FROM curated_target ORDER BY target_id').fetchall():
    secs = {s['section_no']: s['body_md'] for s in con.execute('SELECT section_no, body_md FROM curated_section WHERE target_id=?', (r['target_id'],)).fetchall()}
    curated.append({**dict(r), 'sections': secs})
cats = [dict(r) for r in cur.execute('SELECT * FROM category_overview')]
denovo = [dict(r) for r in cur.execute('SELECT den_no, product_name, product_code, classification_name, specialty, decision_date FROM denovo_catalog')]
pma = [dict(r) for r in cur.execute('SELECT pma_no, product_name, generic_name, product_code, decision_date, decision_code, status_note FROM pma_catalog')]
cdx = [dict(r) for r in cur.execute('SELECT record_id, biomarker, variant_detail, device, indication_sample, submission_no, pathway FROM cdx')]
bundle = {'generated': __import__('datetime').date.today().isoformat(), 'markers': markers, 'links': links, 'submissions': subs, 'extraction': ext,
          'product_codes': pcs, 'curated': curated, 'categories': cats, 'denovo': denovo, 'pma': pma, 'cdx': cdx}
p = os.path.join(OUT, 'bundle.json')
json.dump(bundle, open(p, 'w', encoding='utf-8'), ensure_ascii=False, separators=(',', ':'))
print('bundle.json', round(os.path.getsize(p) / 1e6, 1), 'MB; submissions', len(subs), 'with extraction', len(ext))
con.close()
