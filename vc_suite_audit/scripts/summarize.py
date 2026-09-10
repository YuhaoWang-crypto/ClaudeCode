import glob, os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
R = HERE + "/results"

frames = [pd.read_csv(f) for f in glob.glob(R + "/splits_*.csv")]
df = pd.concat(frames, ignore_index=True).drop_duplicates(
    subset=["dataset", "split", "qrasar"], keep="last")
order = ["random", "scaffold", "cluster0.4", "temporal"]
df["split"] = pd.Categorical(df.split, order, ordered=True)
df = df.sort_values(["dataset", "split", "qrasar"])
df.to_csv(R + "/ALL_splits.csv", index=False)

pd.set_option("display.width", 250)
print("=== Spearman by split severity (independent re-implementation) ===")
piv = df.pivot_table(index="dataset", columns=["split"], values="spearman", observed=True)
print(piv.round(3).to_string())

print("\n=== q-RASAR delta (with - without), by split ===")
d2 = df.pivot_table(index=["dataset", "split"], columns="qrasar", values="spearman", observed=True)
d2["delta"] = d2[True] - d2[False]
print(d2.round(3).to_string())

print("\n=== decision metrics, strictest split (cluster0.4, no qRASAR) ===")
sub = df[(df.split == "cluster0.4") & (~df.qrasar)]
print(sub[["dataset", "n", "spearman", "rmse", "r2", "ef5", "ef1",
           "within1log", "within2fold"]].round(3).to_string(index=False))

print("\n=== applicability-domain curves ===")
ad = pd.concat([pd.read_csv(f) for f in glob.glob(R + "/ad_*.csv")], ignore_index=True)
ad = ad.drop_duplicates(subset=["dataset", "lo"], keep="last")
ad.to_csv(R + "/ALL_ad.csv", index=False)
print(ad.pivot_table(index="dataset", columns="lo", values="spearman").round(3).to_string())
print("\n(n per bin)")
print(ad.pivot_table(index="dataset", columns="lo", values="n").to_string())
