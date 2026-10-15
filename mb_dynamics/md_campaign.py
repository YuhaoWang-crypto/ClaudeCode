"""Run batches of CG Brownian-dynamics jobs in parallel (one OpenMM CPU thread per job).

    python md_campaign.py dt --dts 1.0,0.5,0.25     # time-step selection (+ production at the chosen dt)
    python md_campaign.py sensitivity --dt 0.5      # pre-declared linker geometry sensitivity
Trajectory descriptor files go to $MB_MD_DIR (default: ./md_runs, git-ignored); a JSON log of
wall time per job is appended to results/md_log.json (counts toward the 24 h compute cap).
"""
import argparse
import json
import os
import time
from multiprocessing import Pool

import cg_md

HERE = os.path.dirname(os.path.abspath(__file__))
MD_DIR = os.environ.get("MB_MD_DIR", os.path.join(HERE, "md_runs"))
LOG = os.path.join(HERE, "results", "md_log.json")


def _job(j):
    out = os.path.join(MD_DIR, j["name"] + ".npz")
    if os.path.exists(out):
        return dict(j, skipped=True)
    t = time.time()
    meta = cg_md.run(out, j["seed"], j["mode"], j["dt"], j["t_ns"], j["t_eq_ns"], j.get("save_ps", 10.0),
                     n_rep=j.get("n_rep", 64), n_link=j.get("n_link", 3), lp_link=j.get("lp_link", 1.0))
    return dict(j, wall_s=time.time() - t, sim_wall_s=meta["wall_s"])


# Every job = one OpenMM context with 64 independent replicas (alternating extended/collapsed starts).
# The runs at the selected time step double as production data: seeds 101, 102 build the R2 state
# model, 103 is the held-out validation seed, 104 an extra block for equilibrium agreement checks.

def jobs_dt(dts=(1.0, 0.5, 0.25), t_ns=1000.0, t_eq_ns=100.0):
    return [dict(name=f"dt{dt}_s{s}", seed=s, mode="mixed", dt=dt, t_ns=t_ns, t_eq_ns=t_eq_ns)
            for dt in dts for s in (101, 102, 103, 104)]


def jobs_sensitivity(dt, t_ns=500.0, t_eq_ns=100.0):
    js = []
    for n_link in (2, 6):            # total MB linker 0.6 / 1.8 nm (default 3 x 0.3 = 0.9 nm)
        js.append(dict(name=f"sens_link{n_link}_dt{dt}", seed=200 + n_link, mode="mixed", dt=dt,
                       t_ns=t_ns, t_eq_ns=t_eq_ns, n_link=n_link))
    for lp in (0.7, 2.0):           # linker persistence length (default 1 nm)
        js.append(dict(name=f"sens_lp{lp}_dt{dt}", seed=300 + int(10 * lp), mode="mixed", dt=dt,
                       t_ns=t_ns, t_eq_ns=t_eq_ns, lp_link=lp))
    return js


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["dt", "sensitivity"])
    ap.add_argument("--dts", type=str, default="1.0,0.5,0.25")
    ap.add_argument("--dt", type=float, default=0.5)
    ap.add_argument("--t_ns", type=float, default=None)
    ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(MD_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    kw = {} if a.t_ns is None else {"t_ns": a.t_ns}
    if a.what == "dt":
        js = jobs_dt(tuple(float(x) for x in a.dts.split(",")), **kw)
    else:
        js = jobs_sensitivity(a.dt, **kw)
    t0 = time.time()
    with Pool(a.procs) as p:
        res = p.map(_job, js, chunksize=1)
    log = json.load(open(LOG)) if os.path.exists(LOG) else []
    log.append(dict(batch=a.what, started=time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(t0)),
                    wall_s=time.time() - t0, cpu_s=sum(r.get("wall_s", 0) for r in res), jobs=res))
    json.dump(log, open(LOG, "w"), indent=1)
    print(json.dumps(log[-1]["wall_s"]), "s wall;", sum(r.get("wall_s", 0) for r in res), "s cpu")
