"""Partial-genome polygenic risk scores on 1000 Genomes from public PGS Catalog weights.
Shows the SAME weights give very different score distributions across ancestries - so a single
person's raw score is uninterpretable without an ancestry-matched reference, and population
calibration/AUC cannot be validated from one sample (the matrix row's stated boundary)."""
import os, glob, subprocess, json, numpy as np, pandas as pd
R = os.environ["RB_DATA"]; OUT = "results"
pgs = {"PGS000020": "T2D", "PGS002605": "LDL", "PGS000297": "Height", "PGS000337": "CAD"}
W = {k: pd.read_csv(f"{R}/weights/{k}.tsv", sep="\t", dtype={"chr": str}).set_index("id") for k in pgs}
W = {k: w[~w.index.duplicated()] for k, w in W.items()}
panel = pd.read_csv(f"{R}/panel.txt", sep="\t"); pop = panel.set_index("sample")["super_pop"]

done = sorted(int(open(f).read().split()[0][3:]) for f in glob.glob(f"{R}/chrvcf/chr*.done")
              if int(open(f).read().split()[1]) > 0)
vcf0 = f"{R}/chrvcf/chr{done[0]}.vcf.gz"
samples = subprocess.run(["bcftools", "query", "-l", vcf0], capture_output=True, text=True).stdout.split()
assert len(samples) > 100, f"no samples: {len(samples)}"
acc = {k: np.zeros(len(samples)) for k in pgs}; nvar = {k: 0 for k in pgs}
for c in done:
    out = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
                          f"{R}/chrvcf/chr{c}.vcf.gz"], capture_output=True, text=True).stdout
    for line in out.strip().split("\n"):
        if not line: continue
        f = line.split("\t"); ch, ps, ref, alt = f[0], f[1], f[2], f[3]
        dos = np.array([g[:3].count("1") if "." not in g[:3] else np.nan for g in f[4:]], float)
        key = f"{ch}:{ps}"
        for k in pgs:
            if key in W[k].index:
                row = W[k].loc[key]; ea, oa, wt = row["ea"], str(row["oa"]), float(row["w"])
                has_oa = oa not in ("", "nan", "None")
                if ea == alt and (oa == ref or not has_oa): d = dos.copy()
                elif ea == ref and (oa == alt or not has_oa): d = 2 - dos
                else: continue
                d[np.isnan(d)] = np.nanmean(d) if np.isfinite(np.nanmean(d)) else 0.0
                acc[k] += wt * d; nvar[k] += 1
df = pd.DataFrame(acc, index=samples); df["super_pop"] = pop.reindex(samples).values
df = df.dropna(subset=["super_pop"]); df.to_csv(f"{OUT}/prs_scores.csv")
summ = {"chromosomes_used": done, "n_samples": int(len(df)),
        "variants_scored": {pgs[k]: int(nvar[k]) for k in pgs},
        "variants_total_in_panel": {pgs[k]: int(len(W[k])) for k in pgs}, "by_ancestry": {}}
for k, name in pgs.items():
    eur = df[df.super_pop == "EUR"][k]; mu, sd = eur.mean(), eur.std()
    z = (df[k] - mu) / sd
    summ["by_ancestry"][name] = {sp: {"mean_PRS_z_vs_EUR": round(float(z[df.super_pop == sp].mean()), 2),
        "sd": round(float(z[df.super_pop == sp].std()), 2),
        "pct_above_EUR_95th_pctile": round(float((z[df.super_pop == sp] > 1.645).mean()) * 100, 1)}
        for sp in ["EUR", "SAS", "AMR", "EAS", "AFR"]}
json.dump(summ, open(f"{OUT}/prs_summary.json", "w"), indent=1); print(json.dumps(summ, indent=1))
