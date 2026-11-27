"""Harvest cis-pQTL leads from UKB-PPP (Olink), parallelised (remote tabix is latency-bound)."""
import os, sys, re, csv, json, subprocess, urllib.request
from concurrent.futures import ThreadPoolExecutor
Q = os.environ["QB_DATA"]; OUT = "results"
gene38 = json.load(open(f"{Q}/gene38.json"))
def key(s): return " ".join(re.sub(r"[^a-z0-9 ]", " ", s.lower()).split())
nm = {}
for x in csv.DictReader(open(f"{Q}/hgnc.tsv"), delimiter="\t"):
    for v in [x["name"], x["symbol"]] + x.get("alias_name", "").split("|") + x.get("prev_name", "").split("|"):
        if v.strip(): nm.setdefault(key(v), x["symbol"])
def folder(n): lo = (n - 1) // 1000 * 1000 + 1; return f"GCST{lo:08d}-GCST{(lo+999):08d}"
BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
def fetch(url, timeout=25):
    for _ in range(3):
        try: return urllib.request.urlopen(url, timeout=timeout).read()
        except Exception: pass
    return None
def work(n):
    acc = f"GCST{n:08d}"; fol = folder(n)
    meta = fetch(f"{BASE}/{fol}/{acc}/{acc}.tsv.gz-meta.yaml", 20)
    if not meta: return None
    m = re.search(r"trait_description:\s*\n\s*-\s*(.+)", meta.decode("utf-8", "ignore"))
    if not m: return None
    prot = m.group(1).strip(); g = nm.get(key(re.sub(r"\s+levels?$", "", prot, flags=re.I)))
    if not g or g not in gene38: return None
    c, tss, strand = gene38[g]; H = f"{BASE}/{fol}/{acc}/harmonised/{acc}.h.tsv.gz"
    try:
        r = subprocess.run(["tabix", H, f"{c}:{max(1,tss-500000)}-{tss+500000}"], capture_output=True, text=True, timeout=150).stdout
    except Exception: return None
    best = None
    for line in r.strip().split("\n"):
        if not line: continue
        f = line.split("\t")
        try: p = float(f[7])
        except Exception: continue
        if best is None or p < best[0]: best = (p, f)
    if best is None: return None
    p, f = best
    return f"{acc}\t{prot}\t{g}\t{c}\t{tss}\t{f[1]}\t{f[9]}\t{f[8]}\t{f[4]}\t{f[5]}\t{f[6]}\t{f[7]}\t{r.count(chr(10))}"
start = int(sys.argv[1]); count = int(sys.argv[2])
out = open(f"{OUT}/cis_leads.tsv", "w")
out.write("acc\tprotein\tgene\tchr\ttss\tlead_pos\tlead_rsid\tvariant_id\tbeta\tse\teaf\tp\tn_cis\n")
done = 0
with ThreadPoolExecutor(max_workers=14) as ex:
    for res in ex.map(work, range(start, start + count)):
        done += 1
        if res: out.write(res + "\n"); out.flush()
        if done % 50 == 0: print(f"processed {done}/{count}, leads {sum(1 for _ in open(f'{OUT}/cis_leads.tsv'))-1}", flush=True)
print("DONE", flush=True)
