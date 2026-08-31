"""Extract per-design experimental SHAPE + metadata from OpenKnotBench into a compact pickle."""
import csv, pickle, sys
import numpy as np
import pandas as pd

CSV = "/home/user/work/data-okb/Data/OpenKnotBench_data.v4.5.1.csv"
OUT = "/home/user/work/okb_designs.pkl"

hdr = next(csv.reader(open(CSV)))
react_cols = [h for h in hdr if h.startswith("reactivity_") and "error" not in h]
meta = ["id", "sequence", "SN_filter", "signal_to_noise", "reads", "round", "puzzle",
        "method", "target_openknot_score", "sub_start", "sub_end", "design_length",
        "design_sequence", "target_structure", "RNet_structure", "RNet_F1",
        "RNet_F1_crossed_pair"]

rows = []
for ch in pd.read_csv(CSV, usecols=meta + react_cols, chunksize=3000, low_memory=False):
    ch = ch[ch["design_sequence"].notna() & ch["target_structure"].notna()]
    for _, r in ch.iterrows():
        ds, ts = str(r["design_sequence"]), str(r["target_structure"])
        if len(ds) != len(ts):
            continue
        if set(ds) - set("ACGU"):
            continue
        s, e = int(r["sub_start"]) - 1, int(r["sub_end"])
        if e - s != len(ds):
            continue
        exp = r[react_cols].values.astype(np.float32)[s:e]
        rows.append(dict(
            id=r["id"], round=int(r["round"]), puzzle=r["puzzle"], method=r["method"],
            oks=float(r["target_openknot_score"]) if pd.notna(r["target_openknot_score"]) else np.nan,
            sn=float(r["signal_to_noise"]), sn_filter=int(r["SN_filter"]),
            reads=float(r["reads"]),
            rnet_f1=float(r["RNet_F1"]) if pd.notna(r["RNet_F1"]) else np.nan,
            rnet_f1_cp=float(r["RNet_F1_crossed_pair"]) if pd.notna(r["RNet_F1_crossed_pair"]) else np.nan,
            seq=ds, struct=ts, exp_react=exp,
        ))
    print(f"\rrows={len(rows)}", end="", file=sys.stderr, flush=True)

print(file=sys.stderr)
pickle.dump(rows, open(OUT, "wb"))
print(f"saved {len(rows)} designs -> {OUT}")
df = pd.DataFrame([{k: v for k, v in r.items() if k != "exp_react"} for r in rows])
print(df.groupby(["round"]).size())
print(df.groupby(["round", "method"]).size().to_string())
