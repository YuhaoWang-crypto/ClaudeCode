#!/bin/bash
PQ=$1; GENE38=$2
python3 - "$PQ" "$GENE38" <<'PY'
import json,sys,subprocess,os,re
pq,gene38f=sys.argv[1],sys.argv[2]
acc=json.load(open(f"{pq}/acc.json")); gene38=json.load(open(gene38f))
def fol(a): n=int(a[4:]); lo=(n-1)//1000*1000+1; return f"GCST{lo:08d}-GCST{lo+999:08d}"
B="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
for g,a in acc.items():
    if g not in gene38: continue
    out=f"{pq}/{g}.pqtl.tsv"
    if os.path.exists(out) and os.path.getsize(out)>1000 and sum(1 for _ in open(out))>100:
        print(g,"exists",flush=True); continue
    c,tss,strand=gene38[g]; lo=tss-500000; hi=tss+500000
    hdir=f"{B}/{fol(a)}/{a}/harmonised/"
    fn=None
    for _ in range(4):
        r=subprocess.run(["curl","-s","-m","90",hdir],capture_output=True,text=True)
        mm=re.findall(rf'35078996-{a}-[A-Z_0-9]+\.h\.tsv\.gz', r.stdout)
        if mm: fn=mm[0]; break
    if not fn: print(g,"NO FILENAME",flush=True); continue
    url=hdir+fn
    tmp=f"{pq}/_{g}.gz"
    if os.path.exists(tmp): os.remove(tmp)
    ok=False
    for _ in range(6):
        subprocess.run(["curl","-s","-m","500","-C","-","-o",tmp,url])
        if os.path.exists(tmp) and os.path.getsize(tmp)>1e8 and \
           subprocess.run(["bash","-c",f"zcat {tmp} 2>/dev/null | tail -1 >/dev/null"]).returncode==0:
            ok=True; break
    if not ok:
        if os.path.exists(tmp): os.remove(tmp)
        print(g,"DL FAIL",flush=True); continue
    cmd=f"zcat {tmp} | awk -F'\\t' 'NR==1 || ($3==\"{c}\" && $4>={lo} && $4<={hi})' > {out}"
    subprocess.run(["bash","-c",cmd]); os.remove(tmp)
    print(g,a,fn.split('-')[-1].split('.')[0],"cis rows",sum(1 for _ in open(out))-1,flush=True)
PY
