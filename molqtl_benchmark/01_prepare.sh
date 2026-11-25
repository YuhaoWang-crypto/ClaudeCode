#!/bin/bash
# GEUVADIS (Lappalainen 2013): 1000 Genomes LCLs, chr21 + chr22 genotypes and PEER-normalised gene expression.
set -e; cd $MQ_WORK
for c in 21 22; do plink1.9 --vcf $MQ_DATA/chr$c.vcf.gz --snps-only just-acgt --biallelic-only strict --geno 0.05 --maf 0.01 \
    --set-missing-var-ids @:#  --make-bed --out g$c --const-fid 0 > /dev/null; done
plink1.9 --bfile g21 --bmerge g22 --make-bed --out geno > /dev/null
python3 - <<'PY'
import os, numpy as np, pandas as pd
D, W = os.environ["MQ_DATA"], os.environ["MQ_WORK"]
s = pd.read_csv(f"{D}/E-GEUV-1.sdrf.txt", sep="\t"); pop = s.drop_duplicates("Source Name").set_index("Source Name")["Characteristics[ancestry category]"]
e = pd.read_csv(f"{D}/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz", sep="\t")
e = e[e.Chr.astype(str).isin(["21", "22"])]
fam = pd.read_csv(f"{W}/geno.fam", sep=r"\s+", header=None)[1]
ids = [i for i in fam if i in e.columns]
genes = e[["TargetID", "Gene_Symbol", "Chr", "Coord"]].set_index("TargetID"); genes.columns = ["symbol", "chr", "tss"]
e.set_index("TargetID")[ids].T.to_csv(f"{W}/expr.tsv", sep="\t"); genes.to_csv(f"{W}/genes.tsv", sep="\t")
rng = np.random.default_rng(0); sets = pd.DataFrame({"pop": pop.reindex(ids)}, index=ids)
eur = sets.index[~sets["pop"].str.contains("Yoruba")].values.copy(); rng.shuffle(eur); n = int(.7 * len(eur))
sets["set"] = "YRI"; sets.loc[eur[:n], "set"] = "EUR_train"; sets.loc[eur[n:], "set"] = "EUR_test"
sets.to_csv(f"{W}/sets.tsv", sep="\t"); print(sets.groupby(["set", "pop"]).size(), len(genes))
pd.Series(ids).to_frame().assign(f=0)[["f", 0]].to_csv(f"{W}/keep.txt", sep=" ", header=False, index=False)
PY
plink1.9 --bfile geno --keep keep.txt --make-bed --out geno > /dev/null; wc -l geno.bim geno.fam
