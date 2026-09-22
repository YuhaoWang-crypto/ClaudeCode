"""Paced, resumable bulk fetch of FDA 510(k) decision summaries (and 510(k) summaries).

Reads the list of submission numbers to fetch from the DB (510(k) with year >= 2003, plus De Novo),
skips numbers already present in txt/, downloads with a browser User-Agent, backs off when the
FDA abuse-detection page is returned, and records every attempt in `document` / fetch_log.tsv.

Usage: python3 db/fetch_ds.py <db path> <txt dir> [min_year] [pace_seconds]
"""
import sys, os, time, sqlite3, subprocess, datetime, re
DB, TXT = sys.argv[1], sys.argv[2]
MIN_YEAR = int(sys.argv[3]) if len(sys.argv) > 3 else 2003
PACE = float(sys.argv[4]) if len(sys.argv) > 4 else 1.2
SKIP = int(sys.argv[6]) if len(sys.argv) > 6 else 0          # start this worker at an offset into the todo list
EXTRACT_ALL = os.environ.get('FETCH_EXTRACT_ALL') == '1'     # by default only the decision summary is text-extracted now
UA = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36'
os.makedirs(TXT, exist_ok=True)
LOG = os.path.join(TXT, 'fetch_log.tsv')

def pdf_text(path):
    from pypdf import PdfReader
    try:
        return '\n'.join((p.extract_text() or '') for p in PdfReader(path).pages)
    except Exception as e:
        return ''

def curl(url, out):
    r = subprocess.run(['curl', '-sSL', '--max-time', '90', '-A', UA, '-H', 'Accept: application/pdf,*/*', '-o', out, '-w', '%{http_code}', url], capture_output=True, text=True)
    code = r.stdout.strip()
    head = open(out, 'rb').read(5) if os.path.exists(out) else b''
    is_pdf = code == '200' and head == b'%PDF-'
    abuse = (not is_pdf) and os.path.exists(out) and b'abuse' in open(out, 'rb').read(4000).lower()
    if not is_pdf and os.path.exists(out): os.remove(out)
    return is_pdf, code, abuse

def fetch_one(no):
    yy = no[1:3] if no.startswith('K') else None
    cands = [('decision_summary', f'https://www.accessdata.fda.gov/cdrh_docs/reviews/{no}.pdf')]
    if yy: cands.append(('510k_summary', f'https://www.accessdata.fda.gov/cdrh_docs/pdf{yy}/{no}.pdf'))
    else: cands.append(('denovo_summary', f'https://www.accessdata.fda.gov/cdrh_docs/pdf{no[3:5]}/{no}.pdf'))
    parts, got = [], {}
    for label, url in cands:
        if label != 'decision_summary' and 'decision_summary' in got and not EXTRACT_ALL:
            continue  # the decision summary is what gets parsed; skip the second download to halve server load
        pdf = os.path.join(TXT, f'{no}_{label.upper()}.pdf')
        for attempt in range(3):
            ok, code, abuse = curl(url, pdf)
            if abuse:
                print(f'  abuse page on {no} {label}; backing off 90s', flush=True); time.sleep(90); continue
            break
        if ok:
            if label == 'decision_summary' or EXTRACT_ALL:
                t = pdf_text(pdf)
            else:
                t = '[text extraction deferred; PDF kept on disk]'
            parts.append(f'===== {label.upper()} {url} ({len(t)} chars) =====\n{t}')
            got[label] = (url, len(t))
        time.sleep(PACE)
    open(os.path.join(TXT, no + '.txt'), 'w', encoding='utf-8').write('\n\n'.join(parts) if parts else '[nothing found]')
    return got

def main():
    con = sqlite3.connect(DB, timeout=120); cur = con.cursor()
    order = 'ASC' if (len(sys.argv) > 5 and sys.argv[5].upper() == 'ASC') else 'DESC'
    todo = [x[0] for x in cur.execute(f"SELECT submission_no FROM submission WHERE (pathway='510(k)' AND year>=?) OR pathway='De Novo' ORDER BY year {order}, submission_no {order}", (MIN_YEAR,)).fetchall()]
    todo = [n for n in todo if not os.path.exists(os.path.join(TXT, n + '.txt'))][SKIP:]
    print('to fetch:', len(todo), 'order', order, 'skip', SKIP, flush=True)
    for i, no in enumerate(todo):
        if os.path.exists(os.path.join(TXT, no + '.txt')):  # another worker got there first
            continue
        got = fetch_one(no)
        now = datetime.datetime.utcnow().isoformat(timespec='seconds')
        for label, (url, n) in got.items():
            cur.execute('INSERT OR REPLACE INTO document VALUES (?,?,?,?,?,?)', (no, label, url, n, os.path.join('txt', no + '.txt'), now))
        cur.execute('UPDATE submission SET has_decision_summary=?, has_510k_summary=?, fetched_at=? WHERE submission_no=?',
                    (1 if 'decision_summary' in got else 0, 1 if ('510k_summary' in got or 'denovo_summary' in got) else 0, now, no))
        con.commit()
        with open(LOG, 'a') as f: f.write(f'{now}\t{no}\t{",".join(got.keys()) or "-"}\n')
        if i % 50 == 0: print(f'  {i+1}/{len(todo)} {no} {list(got.keys())}', flush=True)
    print('done', flush=True)

if __name__ == '__main__':
    main()
