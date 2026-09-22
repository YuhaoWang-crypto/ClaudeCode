"""Enrich `submission` from openFDA and fill the post-2014 gap by product code.

Stage A: for every K number already in `submission`, batch-query openFDA device/510k (50 per call)
         and store decision_date / applicant / device_name / product_code / regulation / etc.
Stage B: for every product_code seen in stage A, pull the classification record and ALL 510(k)s
         with decision_date >= 2003 (paged, 1000 per call); add any K number not yet in `submission`
         (in_catalog=0) so the database is current to the openFDA snapshot.

Usage: python3 db/openfda_enrich.py [db path]
"""
import sys, os, json, time, sqlite3, urllib.request, urllib.parse, re
DB = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fda_ivd_markers.sqlite')
UA = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36'
BASE = 'https://api.fda.gov/device/'

def api(path, params, retries=4):
    url = BASE + path + '?' + urllib.parse.urlencode(params)
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA})
            return json.load(urllib.request.urlopen(req, timeout=60))
        except urllib.error.HTTPError as e:
            if e.code == 404: return {'results': [], 'meta': {'results': {'total': 0}}}
            if e.code == 429: time.sleep(10 * (i + 1)); continue
            time.sleep(2 * (i + 1))
        except Exception:
            time.sleep(2 * (i + 1))
    return {'results': [], 'meta': {'results': {'total': 0}}}

def upsert_510k(cur, r, in_catalog=None):
    no = r.get('k_number')
    if not no: return
    d = r.get('decision_date') or ''
    year = int(d[:4]) if d[:4].isdigit() else None
    cur.execute('INSERT OR IGNORE INTO submission (submission_no, pathway, year, in_catalog) VALUES (?,?,?,?)', (no, 'De Novo' if no.startswith('DEN') else '510(k)', year, 0 if in_catalog is None else in_catalog))
    cur.execute('''UPDATE submission SET decision_date=?, decision=?, applicant=?, device_name=?, product_code=?, regulation_number=?,
                   advisory_committee=?, statement_or_summary=?, clearance_type=?, third_party=?, expedited=?, openfda_found=1, year=COALESCE(?, year)
                   WHERE submission_no=?''',
                (d, r.get('decision_description'), r.get('applicant'), r.get('device_name'), r.get('product_code'),
                 (r.get('openfda') or {}).get('regulation_number'), r.get('advisory_committee_description'), r.get('statement_or_summary'),
                 r.get('clearance_type'), r.get('third_party_flag'), r.get('expedited_review_flag'), year, no))

def stage_a(con):
    cur = con.cursor()
    ks = [x[0] for x in cur.execute("SELECT submission_no FROM submission WHERE pathway='510(k)' AND openfda_found=0").fetchall()]
    print('stage A: K numbers to enrich', len(ks)); found = 0
    for i in range(0, len(ks), 50):
        batch = ks[i:i+50]
        q = 'k_number:(' + ' OR '.join(batch) + ')'
        d = api('510k.json', {'search': q, 'limit': 100})
        for r in d.get('results', []):
            upsert_510k(cur, r, in_catalog=1); found += 1
        con.commit()
        if (i // 50) % 20 == 0: print(f'  {i+len(batch)}/{len(ks)} found={found}', flush=True)
        time.sleep(0.25)
    print('stage A done, found', found)

def stage_b(con, min_year=2003):
    cur = con.cursor()
    codes = [x[0] for x in cur.execute("SELECT DISTINCT product_code FROM submission WHERE product_code IS NOT NULL AND product_code<>''").fetchall()]
    print('stage B: product codes', len(codes)); added = 0
    for j, code in enumerate(codes):
        c = api('classification.json', {'search': f'product_code:{code}', 'limit': 1})
        for r in c.get('results', []):
            cur.execute('INSERT OR REPLACE INTO product_code VALUES (?,?,?,?,?,?,?,?)',
                        (code, r.get('device_name'), r.get('regulation_number'), r.get('device_class'), r.get('medical_specialty_description'),
                         r.get('review_panel'), len((r.get('openfda') or {}).get('k_number', []) or []), r.get('definition')))
        skip = 0
        while True:
            d = api('510k.json', {'search': f'product_code:{code} AND decision_date:[{min_year}-01-01 TO 2030-12-31]', 'limit': 1000, 'skip': skip})
            res = d.get('results', [])
            for r in res:
                before = cur.execute('SELECT openfda_found FROM submission WHERE submission_no=?', (r.get('k_number'),)).fetchone()
                upsert_510k(cur, r)
                if before is None: added += 1
            total = d.get('meta', {}).get('results', {}).get('total', 0)
            skip += len(res)
            if not res or skip >= total or skip >= 5000: break
            time.sleep(0.25)
        con.commit()
        if j % 25 == 0: print(f'  {j+1}/{len(codes)} codes, added={added}', flush=True)
        time.sleep(0.25)
    print('stage B done, added', added)

if __name__ == '__main__':
    con = sqlite3.connect(DB, timeout=120)
    stage_a(con)
    stage_b(con)
    cur = con.cursor()
    print('submissions total', cur.execute('SELECT COUNT(*) FROM submission').fetchone()[0],
          'with openfda', cur.execute('SELECT COUNT(*) FROM submission WHERE openfda_found=1').fetchone()[0],
          'product codes', cur.execute('SELECT COUNT(*) FROM product_code').fetchone()[0])
    con.close()
