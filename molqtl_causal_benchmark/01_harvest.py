"""For each validated LDL locus, pull the GWAS region and eQTL (liver+blood) region via remote tabix."""
import os, json, subprocess, numpy as np, pandas as pd
OUT = "results"; W = os.environ["CW_DIR"]; os.makedirs(f"{W}/loci", exist_ok=True)
loci = json.load(open(f"{OUT}/causal_loci.json"))
GWAS = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90239001-GCST90240000/GCST90239658/harmonised/GCST90239658.h.tsv.gz"
EQTL = {"liver": "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000266/QTD000266.all.tsv.gz",
        "blood": "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000356/QTD000356.all.tsv.gz"}
WIN = 500000
def tabix(url, region):
    for _ in range(3):
        r = subprocess.run(["tabix", url, region], capture_output=True, text=True, timeout=180)
        if r.returncode == 0: return r.stdout
    return ""
for g, v in loci.items():
    c, tss = v["chr"], v["tss"]; region = f"{c}:{max(1,tss-WIN)}-{tss+WIN}"
    gw = tabix(GWAS, region)
    gwl = [l.split("\t") for l in gw.strip().split("\n") if l]
    gdf = pd.DataFrame(gwl, columns=["chr","pos","ea","oa","beta","se","eaf","p","vid","rsid","n","hm","code","rs2"][:len(gwl[0])] if gwl else None)
    if gdf.empty: print(g, "no GWAS"); continue
    gdf = gdf[["chr","pos","ea","oa","beta","se","eaf","p","rsid"]].copy()
    gdf.to_csv(f"{W}/loci/{g}.gwas.tsv", sep="\t", index=False)
    for tis, url in EQTL.items():
        eq = tabix(url, region)
        rows = [l.split("\t") for l in eq.strip().split("\n") if l]
        # eQTL all.tsv cols: molecular_trait_id,chr,pos,ref,alt,variant,ma_samples,maf,pvalue,beta,se,type,ac,an,r2,mto,gene_id,median_tpm,rsid
        ed = pd.DataFrame(rows)
        if ed.empty: continue
        ed = ed.iloc[:, [0,1,2,3,4,7,8,9,10,16,18]]
        ed.columns = ["gene","chr","pos","ref","alt","maf","p","beta","se","gene_id","rsid"]
        ed.to_csv(f"{W}/loci/{g}.{tis}.tsv", sep="\t", index=False)
    gsig = gdf.p.astype(float).min()
    print(g, "GWAS minP", f"{gsig:.1e}", "| genes(liver)", pd.read_csv(f"{W}/loci/{g}.liver.tsv",sep='\t').gene.nunique() if os.path.exists(f"{W}/loci/{g}.liver.tsv") else 0, flush=True)
