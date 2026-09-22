"""Helpers for FDA 510(k) research.
Usage:
  python3 fda_fetch.py list <PRODUCT_CODE> [limit]      -> recent clearances (K number, date, device, applicant)
  python3 fda_fetch.py get  <KNUMBER> [KNUMBER ...]     -> download decision summary + 510(k) summary, save text to txt/<K>.txt, print char count
  python3 fda_fetch.py show <KNUMBER> [maxchars]        -> print saved text
"""
import json, sys, os, urllib.request, urllib.parse, subprocess, re
UA = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36'
HERE = os.path.dirname(os.path.abspath(__file__)); TXT = os.path.join(HERE, 'txt'); os.makedirs(TXT, exist_ok=True)

def api(url):
    req = urllib.request.Request(url, headers={'User-Agent': UA})
    return json.load(urllib.request.urlopen(req, timeout=60))

def list_code(code, limit=40):
    q = urllib.parse.quote(f'product_code:{code}')
    d = api(f'https://api.fda.gov/device/510k.json?search={q}&sort=decision_date:desc&limit={limit}')
    print(f"product_code {code}: total={d['meta']['results']['total']}")
    for r in d['results']:
        print(f"{r['k_number']} | {r['decision_date']} | {r.get('decision_description','')[:22]} | {r['device_name'][:70]} | {r['applicant'][:35]} | {r.get('statement_or_summary','')}")

def curl(url, out):
    r = subprocess.run(['curl','-sSL','--max-time','90','-A',UA,'-H','Accept: application/pdf,*/*','-o',out,'-w','%{http_code}',url], capture_output=True, text=True)
    ok = r.stdout.strip()=='200' and os.path.getsize(out)>1000 and open(out,'rb').read(5)==b'%PDF-'
    if not ok and os.path.exists(out): os.remove(out)
    return ok

def pdf_text(path):
    from pypdf import PdfReader
    try:
        return '\n'.join((p.extract_text() or '') for p in PdfReader(path).pages)
    except Exception as e:
        return f'[pdf error {e}]'

def get(k):
    k = k.upper(); yy = k[1:3]
    out = os.path.join(TXT, k+'.txt'); parts = []
    cands = [('DECISION_SUMMARY', f'https://www.accessdata.fda.gov/cdrh_docs/reviews/{k}.pdf'),
             ('510K_SUMMARY', f'https://www.accessdata.fda.gov/cdrh_docs/pdf{yy}/{k}.pdf'),
             ('510K_SUMMARY', f'https://www.accessdata.fda.gov/cdrh_docs/pdf/{k}.pdf')]
    seen=set()
    for label, url in cands:
        if label in seen: continue
        pdf = os.path.join(TXT, f'{k}_{label}.pdf')
        if curl(url, pdf):
            seen.add(label); t = pdf_text(pdf)
            parts.append(f'===== {label} {url} ({len(t)} chars) =====\n{t}')
    open(out,'w').write('\n\n'.join(parts) if parts else '[nothing found]')
    print(f'{k}: {[p.split()[1] for p in parts]} -> {out} ({sum(len(p) for p in parts)} chars)')

if __name__=='__main__':
    cmd=sys.argv[1]
    if cmd=='list': list_code(sys.argv[2], int(sys.argv[3]) if len(sys.argv)>3 else 40)
    elif cmd=='get':
        for k in sys.argv[2:]: get(k)
    elif cmd=='show':
        t=open(os.path.join(TXT, sys.argv[2].upper()+'.txt')).read(); n=int(sys.argv[3]) if len(sys.argv)>3 else 6000; print(t[:n])
