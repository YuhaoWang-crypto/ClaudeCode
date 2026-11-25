#!/bin/bash
# Natural-population GWAS on the Rice Diversity Panel (Zhao et al. 2011; 413 accessions, 44K SNP chip).
# Same SNPs, same traits, two models: naive linear regression vs GEMMA LMM with a kinship matrix.
set -e
R=$QTL_DATA/rice; W=$QTL_WORK/rice; mkdir -p $W; cd $W
plink1.9 --file $R/RiceDiversity_44K_Genotypes_PLINK/sativas413 --maf 0.05 --geno 0.1 --make-bed --out panel > /dev/null
awk '{$6=1; print}' panel.fam > fam.tmp && mv fam.tmp panel.fam   # GEMMA drops samples whose .fam phenotype is -9
plink1.9 --bfile panel --recode A --out panel > /dev/null          # 0/1/2 matrix for the GS script
plink1.9 --bfile panel --pca 4 --out panel > /dev/null
python3 - "$R/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt" <<'PY'
import sys, pandas as pd
p = pd.read_csv(sys.argv[1], sep="\t"); p.columns = [c.strip() for c in p.columns]
fam = pd.read_csv("panel.fam", sep=r"\s+", header=None)
p = fam[[0, 1]].merge(p, left_on=1, right_on="NSFTVID", how="left")
traits = ["Amylose content", "Seed length", "Seed width", "Pericarp color", "Plant height", "Flowering time at Arkansas"]
p[traits].to_csv("pheno.txt", sep="\t", header=False, index=False, na_rep="NA")
open("traits.txt", "w").write("\n".join(traits))
PY
gemma -bfile panel -gk 1 -o kin -outdir . > /dev/null
for i in 1 2 3 4 5 6; do
  gemma -bfile panel -p pheno.txt -n $i -lm 1 -o lm_$i -outdir . > /dev/null
  gemma -bfile panel -p pheno.txt -n $i -k kin.cXX.txt -lmm 1 -o lmm_$i -outdir . > /dev/null
done
