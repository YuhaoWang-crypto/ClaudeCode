"""
PhysioNet Fantasia (open access, Open Data Commons Attribution v1.0):
  Iyengar N, Peng C-K, Morin R, Goldberger AL, Lipsitz LA. Age-related
  alterations in the fractal scaling of cardiac interbeat interval dynamics.
  Am J Physiol 1996;271:R1078-84. https://physionet.org/content/fantasia/1.0.0/

20 young (21-34 y) and 20 elderly (68-85 y) rigorously screened healthy
subjects, 120 min supine rest; ECG + respiration at 250 Hz, with verified beat
annotations (.ecg). We keep: beat times / RR intervals and the respiration
signal down-sampled to 4 Hz. Cache: figures/_data_cache/fantasia/.
"""
import os
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from netjac import DATA_CACHE

BASE = "https://physionet.org/files/fantasia/1.0.0"
DIR = os.path.join(DATA_CACHE, "fantasia")
RECORDS = ([f"f1y{i:02d}" for i in range(1, 11)] +
           [f"f2y{i:02d}" for i in range(1, 11)] +
           [f"f1o{i:02d}" for i in range(1, 11)] +
           [f"f2o{i:02d}" for i in range(1, 11)])


def _curl(url, dest, tries=4):
    import time
    for k in range(tries):
        r = subprocess.run(["curl", "-sf", "--max-time", "300", url, "-o",
                            dest])
        if r.returncode == 0 and os.path.getsize(dest) > 0:
            return
        time.sleep(2 ** (k + 1))
    raise RuntimeError(f"download failed: {url}")


def fetch_record(rec):
    import wfdb
    from scipy.signal import resample_poly
    os.makedirs(DIR, exist_ok=True)
    out = os.path.join(DIR, f"{rec}.npz")
    if os.path.exists(out):
        return out
    with tempfile.TemporaryDirectory() as tmp:
        for ext in ("hea", "dat", "ecg"):
            _curl(f"{BASE}/{rec}.{ext}", os.path.join(tmp, f"{rec}.{ext}"))
        hea = open(os.path.join(tmp, f"{rec}.hea")).read()
        m = re.search(r"Age:\s*(\d+)\s*Sex:\s*(\w)", hea)
        age, sex = (int(m.group(1)), m.group(2)) if m else (np.nan, "?")
        r = wfdb.rdrecord(os.path.join(tmp, rec))
        a = wfdb.rdann(os.path.join(tmp, rec), "ecg")
        fs = r.fs
        names = [n.upper() for n in r.sig_name]
        resp = r.p_signal[:, names.index("RESP")]
        beats = a.sample[np.isin(a.symbol, ["N"])]
        # exact rational resampling fs -> 4 Hz (anti-aliased)
        resp4 = resample_poly(np.nan_to_num(resp), 4, int(fs))
    np.savez_compressed(out, beats=beats / fs, resp=resp4, fs_resp=4.0,
                        age=age, sex=sex,
                        group="young" if "y" in rec[2] else "elderly")
    return out


def fetch_all(workers=6):
    def one(r):
        try:
            return r, fetch_record(r)
        except Exception as e:                       # noqa: BLE001
            return r, f"ERR {type(e).__name__}: {e}"
    with ThreadPoolExecutor(workers) as ex:
        return dict(ex.map(one, RECORDS))


if __name__ == "__main__":
    res = fetch_all()
    bad = {k: v for k, v in res.items() if str(v).startswith("ERR")}
    print(f"{len(res) - len(bad)}/{len(res)} records ok", bad)
