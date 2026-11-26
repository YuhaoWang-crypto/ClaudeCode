"""Collapse the IEDB tcell_search export (iedb_tcell_I.jsonl) to peptide x allele pairs with
label = ever-Positive vs always-Negative, plus earliest reference year for the leakage split."""
import os, json, numpy as np, pandas as pd
W = os.environ["IB_WORK"]
df = pd.DataFrame(json.loads(l) for l in open(f"{W}/iedb_tcell_I.jsonl"))
df["seq"] = df.linear_sequence.str.upper().str.strip()
df = df[df.seq.str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")].copy()
df["len"] = df.seq.str.len(); df = df[(df.len >= 8) & (df.len <= 11)]
df["allele"] = df.mhc_allele_name.fillna(df.mhc_restriction)
df = df[df.allele.str.match(r"HLA-[ABC]\*\d+:\d+$", na=False)]
df["pos"] = (df.qualitative_measure.str.startswith("Positive")).astype(int)
df["year"] = df.reference_dates.apply(lambda x: int(x[0]) if isinstance(x, list) and x and str(x[0])[:4].isdigit() else np.nan)
g = df.groupby(["seq", "allele", "len"]).agg(pos=("pos", "max"), n=("pos", "size"),
    minyear=("year", "min"), maxyear=("year", "max")).reset_index()
g.to_parquet(f"{W}/tcell_pairs.parquet"); print("pairs", len(g), "pos", int(g.pos.sum()))
