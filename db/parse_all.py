"""Parse every fetched decision summary in a txt directory into the `extraction` table
(one row per submission_no × field) and register the documents that exist on disk.

Usage: python3 db/parse_all.py <db path> <txt dir>
"""
import sys, os, re, glob, sqlite3, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ds_parse import parse_file, split_docs

DB, TXT = sys.argv[1], sys.argv[2]
NA = re.compile(r'^(N/?A|Not applicable\.?|None\.?|Not Applicable\.?)$', re.I)

def main():
    con = sqlite3.connect(DB, timeout=120); cur = con.cursor()
    files = sorted(glob.glob(os.path.join(TXT, '*.txt')))
    n_parsed = n_docs = 0
    for f in files:
        no = os.path.basename(f)[:-4]
        if not re.match(r'^(K|DEN)\d+$', no): continue
        txt = open(f, encoding='utf-8', errors='ignore').read()
        if txt.startswith('[nothing found]'): continue
        now = datetime.datetime.utcnow().isoformat(timespec='seconds')
        for m in re.finditer(r'^===== (DECISION_SUMMARY|510K_SUMMARY|DENOVO_SUMMARY) (\S+) \((\d+) chars\) =====', txt, flags=re.M):
            cur.execute('INSERT OR IGNORE INTO document VALUES (?,?,?,?,?,?)', (no, m.group(1).lower(), m.group(2), int(m.group(3)), os.path.join('txt', no + '.txt'), now)); n_docs += 1
        cur.execute('INSERT OR IGNORE INTO submission (submission_no, pathway, year) VALUES (?,?,?)', (no, 'De Novo' if no.startswith('DEN') else '510(k)', None))
        cur.execute('UPDATE submission SET has_decision_summary=?, has_510k_summary=? WHERE submission_no=?',
                    (1 if 'DECISION_SUMMARY' in txt[:200000] else 0, 1 if ('510K_SUMMARY' in txt or 'DENOVO_SUMMARY' in txt) else 0, no))
        tpl, fields, derived = parse_file(f)
        if tpl == 'none': continue
        cur.execute('DELETE FROM extraction WHERE submission_no=?', (no,))
        for k, v in fields.items():
            v = (v or '').strip()
            if not v or NA.match(v): continue
            cur.execute('INSERT OR REPLACE INTO extraction VALUES (?,?,?,?,?,?)', (no, k, v[:20000], k, tpl, 'section'))
        for k, v in derived.items():
            if v: cur.execute('INSERT OR REPLACE INTO extraction VALUES (?,?,?,?,?,?)', (no, k, v, 'derived', tpl, 'regex'))
        n_parsed += 1
        if n_parsed % 500 == 0: con.commit(); print('  parsed', n_parsed, flush=True)
    con.commit()
    print('documents registered', n_docs, 'summaries parsed', n_parsed,
          'extraction rows', cur.execute('SELECT COUNT(*) FROM extraction').fetchone()[0])
    con.close()

if __name__ == '__main__':
    main()
