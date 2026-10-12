"""
PhysioNet "Autonomic Aging" (open access, ODC-By 1.0 per PhysioNet listing):
  Schumann A, Bär K. Autonomic Aging: A dataset to quantify changes of
  cardiovascular autonomic function during healthy aging (v1.0.0).
  PhysioNet 2021. doi:10.13026/2hsy-t491

1,121 healthy volunteers at rest, ECG (1000 Hz) + continuous non-invasive
finger blood pressure (vascular unloading). Age is given in 5-year groups
(1 = 18-19, 2 = 20-24, ..., 15 = 85-92). Device 0 = Task Force Monitor,
1 = CNAP 500 / MP150 (a potential confounder -> adjust for it).

We stream each WFDB record, detect R peaks, take systolic pressure per beat,
keep only the beat series (time, RR, SBP) and delete the raw waveform.
Cache: figures/_data_cache/autonomic_aging/ (git-ignored).
"""
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from netjac import DATA_CACHE

BASE = "https://physionet.org/files/autonomic-aging-cardiovascular/1.0.0"
DIR = os.path.join(DATA_CACHE, "autonomic_aging")
AGE_YEARS = {1: 18.5, 2: 22, 3: 27, 4: 32, 5: 37, 6: 42, 7: 47, 8: 52, 9: 57,
             10: 62, 11: 67, 12: 72, 13: 77, 14: 82, 15: 88}


def _curl(url, dest, tries=4):
    import time
    for k in range(tries):
        r = subprocess.run(["curl", "-sf", "--max-time", "300", url, "-o",
                            dest])
        if r.returncode == 0 and os.path.getsize(dest) > 0:
            return
        time.sleep(2 ** (k + 1))
    raise RuntimeError(f"download failed: {url}")


def subject_info():
    import pandas as pd
    os.makedirs(DIR, exist_ok=True)
    p = os.path.join(DIR, "subject-info.csv")
    if not os.path.exists(p):
        _curl(f"{BASE}/subject-info.csv", p)
    d = pd.read_csv(p, dtype={"ID": str})
    d["age"] = d["Age_group"].map(AGE_YEARS)
    return d


def _rpeaks(ecg, fs):
    from scipy.signal import butter, filtfilt, find_peaks
    b, a = butter(3, [5 / (fs / 2), 25 / (fs / 2)], "band")
    y = filtfilt(b, a, ecg)
    if np.percentile(-y, 99.5) > np.percentile(y, 99.5):
        y = -y
    thr = 0.35 * np.percentile(y, 99.5)
    pk, _ = find_peaks(y, height=thr, distance=int(0.33 * fs))
    return pk


def beats_from_record(path_noext):
    """Return dict(t, rr, sbp) per beat with artefact rejection."""
    import wfdb
    rec = wfdb.rdrecord(path_noext)
    fs = rec.fs
    names = [n.upper() for n in rec.sig_name]
    sig = rec.p_signal
    bp = sig[:, names.index("NIBP")]
    best = None
    for ch in [i for i, n in enumerate(names) if n.startswith("ECG")]:
        x = sig[:, ch]
        if not np.isfinite(x).all() or np.nanstd(x) == 0:
            continue
        pk = _rpeaks(np.nan_to_num(x), fs)
        rr = np.diff(pk) / fs
        ok = ((rr > 0.3) & (rr < 2.0)).mean() if len(rr) else 0
        if best is None or ok > best[0]:
            best = (ok, pk)
    pk = best[1]
    t = pk[1:] / fs
    rr = np.diff(pk) / fs * 1000.0                         # ms
    # systolic = max BP in the window [R_k, R_k+1) (pulse transit delay)
    sbp = np.array([np.nanmax(bp[a:b]) if b > a else np.nan
                    for a, b in zip(pk[:-1], pk[1:])])
    # artefact rejection: physiological ranges + ectopic/missed beats
    med = np.convolve(rr, np.ones(11) / 11, mode="same")
    good = ((rr > 300) & (rr < 2000) & (np.abs(rr - med) < 0.2 * med) &
            (sbp > 60) & (sbp < 250) & np.isfinite(sbp))
    return {"t": t, "rr": rr, "sbp": sbp, "good": good, "fs": fs}


def fetch_subject(sid, keep_raw=False):
    os.makedirs(DIR, exist_ok=True)
    out = os.path.join(DIR, f"{sid}.npz")
    if os.path.exists(out):
        return out
    with tempfile.TemporaryDirectory() as tmp:
        for ext in ("hea", "dat"):
            _curl(f"{BASE}/{sid}.{ext}", os.path.join(tmp, f"{sid}.{ext}"))
        b = beats_from_record(os.path.join(tmp, sid))
    np.savez_compressed(out, **{k: v for k, v in b.items()})
    return out


def fetch_cohort(ids, workers=6):
    def one(s):
        try:
            return s, fetch_subject(s)
        except Exception as e:                      # noqa: BLE001
            return s, f"ERR {type(e).__name__}: {e}"
    with ThreadPoolExecutor(workers) as ex:
        return dict(ex.map(one, ids))


def default_selection(seed=0):
    """All subjects aged 50+ (groups >= 8), 60 aged 30-49, 120 aged 18-29."""
    d = subject_info().dropna(subset=["Age_group"])
    rng = np.random.default_rng(seed)
    old = d[d.Age_group >= 8].ID.tolist()
    mid = d[(d.Age_group >= 4) & (d.Age_group <= 7)].ID.tolist()
    young = d[d.Age_group <= 3].ID.tolist()
    pick = old + list(rng.choice(mid, 60, replace=False)) + \
        list(rng.choice(young, 120, replace=False))
    return sorted(pick)


if __name__ == "__main__":
    ids = default_selection()
    res = fetch_cohort(ids)
    errs = {k: v for k, v in res.items() if str(v).startswith("ERR")}
    print(f"{len(ids)} selected, {len(ids) - len(errs)} ok, {len(errs)} errors")
    for k, v in list(errs.items())[:10]:
        print(k, v)
