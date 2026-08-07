"""Individual-level longitudinal labs: real per-patient series, not summary CVs.

MIMIC-IV Clinical Database Demo v2.2 (open access, ODC-BY, 100 patients).
Tests the personal-baseline argument on raw per-person trajectories:

  1. variance decomposition -> EMPIRICAL CV_I / CV_G / index of individuality
  2. how many KDIGO-defined AKI events occur while the absolute creatinine is
     still INSIDE the population reference interval (i.e. would be invisible to
     a population cutoff)

Caveat carried into the report: these are ICU patients, so within-person
variation reflects acute illness, not healthy biological variation. The CV_I
here is NOT comparable to a biological-variation study. But it is exactly the
setting where the personal-baseline claim has to pay off.
"""
import numpy as np, pandas as pd

BASE = "/tmp/mimicdemo/mimic-iv-clinical-database-demo-2.2"
CREAT = 50912          # Creatinine, Blood
REF_LO, REF_HI = 0.5, 1.2   # conventional adult serum creatinine reference (mg/dL)

lab = pd.read_csv(f"{BASE}/hosp/labevents.csv.gz",
                  usecols=["subject_id", "itemid", "charttime", "valuenum"],
                  parse_dates=["charttime"])
cr = lab[(lab.itemid == CREAT) & lab.valuenum.notna()].copy()
cr = cr.sort_values(["subject_id", "charttime"])
print(f"creatinine results: {len(cr)} from {cr.subject_id.nunique()} patients")

n_per = cr.groupby("subject_id").size()
print(f"measurements per patient: median {n_per.median():.0f}, "
      f"IQR {n_per.quantile(.25):.0f}-{n_per.quantile(.75):.0f}, max {n_per.max()}")

# ---- 1. variance decomposition on log scale ---------------------------------
keep = n_per[n_per >= 5].index
d = cr[cr.subject_id.isin(keep)].copy()
d["ln"] = np.log(d.valuenum.clip(lower=0.05))
grand = d.ln.mean()
pm = d.groupby("subject_id").ln.mean()
# within-person variance (pooled), between-person variance of the means
within = d.groupby("subject_id").ln.var(ddof=1).mean()
between = pm.var(ddof=1)
cv_i = 100 * np.sqrt(within)
cv_g = 100 * np.sqrt(max(between - within / n_per[keep].mean(), 0))
print(f"\npatients with >=5 results: {len(keep)}")
print(f"  empirical CV_I (within-person) = {cv_i:.1f}%")
print(f"  empirical CV_G (between-person) = {cv_g:.1f}%")
print(f"  index of individuality II = {cv_i/cv_g:.2f}")

# ---- 2. KDIGO AKI events hidden inside the population reference interval ----
# KDIGO: rise >= 0.3 mg/dL within 48 h, OR >= 1.5x lowest value in prior 7 d.
events, hidden, rows = 0, 0, []
for sid, g in cr.groupby("subject_id"):
    g = g.reset_index(drop=True)
    for i in range(len(g)):
        t, v = g.charttime[i], g.valuenum[i]
        prior = g[(g.charttime < t) & (g.charttime >= t - pd.Timedelta("7D"))]
        if prior.empty:
            continue
        p48 = g[(g.charttime < t) & (g.charttime >= t - pd.Timedelta("48h"))]
        crit_abs = (not p48.empty) and (v - p48.valuenum.min() >= 0.3)
        crit_rel = v >= 1.5 * prior.valuenum.min()
        if crit_abs or crit_rel:
            events += 1
            if REF_LO <= v <= REF_HI:          # still "normal" by population RI
                hidden += 1
                rows.append(dict(subject_id=sid, value=v,
                                 baseline=round(prior.valuenum.min(), 2),
                                 ratio=round(v / prior.valuenum.min(), 2)))

print(f"\nKDIGO-style AKI events detected: {events}")
print(f"  of which the ABSOLUTE value was still inside the population "
      f"reference interval ({REF_LO}-{REF_HI} mg/dL): {hidden} "
      f"({100*hidden/events:.1f}%)")
print("\n  examples (value normal by population RI, abnormal for the patient):")
ex = pd.DataFrame(rows).drop_duplicates("subject_id").head(8)
print(ex.to_string(index=False))

out = dict(n_results=int(len(cr)), n_patients=int(cr.subject_id.nunique()),
           n_ge5=int(len(keep)), cv_i=round(float(cv_i), 2),
           cv_g=round(float(cv_g), 2), ii=round(float(cv_i/cv_g), 3),
           kdigo_events=int(events), hidden_in_ri=int(hidden),
           hidden_frac=round(100*hidden/events, 1))
pd.Series(out).to_csv("/home/user/ClaudeCode/recorder_pipeline/data/mimic_personal.tsv",
                      sep="\t", header=False)
print("\nsaved ->", out)
