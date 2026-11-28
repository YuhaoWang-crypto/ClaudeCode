import json,re,csv
from scipy.stats import spearmanr
try:
    from Bio.Align import substitution_matrices; B=substitution_matrices.load("BLOSUM62")
    def bl(a,b):
        try: return B[a,b]
        except: return B[b,a]
except Exception:
    B=None
D="/home/user/ClaudeCode/protein_design_benchmark/data"
seqs=json.load(open(f"{D}/seqs.json"))
assays={"SPG1_STRSG_Olson_2014":"GB1 binding","BLAT_ECOLX_Stiffler_2015":"beta-lactamase","GFP_AEQVI_Sarkisyan_2016":"GFP"}
out={}
for aid in assays:
    seq=seqs[aid]["seq"]; xs=[]; ys=[]
    for r in csv.DictReader(open(f"{D}/{aid}.csv")):
        m=r["mutant"]
        if ":" in m or ";" in m: continue
        mm=re.match(r"^([A-Z])(\d+)([A-Z])$",m)
        if not mm: continue
        wt,pos,mut=mm.group(1),int(mm.group(2)),mm.group(3)
        if pos<1 or pos>len(seq) or seq[pos-1]!=wt: continue
        try: y=float(r["DMS_score"])
        except: continue
        xs.append(bl(wt,mut)); ys.append(y)
    rho=spearmanr(xs,ys).correlation
    out[aid]=round(rho,3); print(f"{aid:28s} BLOSUM62 rho={rho:.3f}  n={len(xs)}")
json.dump(out,open("/home/user/ClaudeCode/protein_design_benchmark/results/blosum62.json","w"))
