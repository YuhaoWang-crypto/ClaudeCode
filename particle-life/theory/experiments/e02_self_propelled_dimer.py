"""
E02 - The self-propelled dimer: where motion is born.
====================================================

Tests the exact two-body result of ``twobody.py`` against direct integration:

    bond length   r*  solves  f_S(r*) = 0      (reciprocal part)
    drift speed   V = mu_discrete * f_A(r*)    (non-reciprocal part)

and demonstrates the degeneracy that makes the *symmetric-radius* chase stall:
when ``alpha_ab = alpha_ba`` and ``beta_ab = beta_ba`` the two force curves share
their zeros, so ``f_A(r*) = 0`` and the pair comes to rest.  Sustained motion
needs either asymmetric radii or three-body frustration.
"""

import numpy as np
from common import ACC, dump, save
import matplotlib.pyplot as plt

from plife import matrices as M, model, twobody as TB

BOX = 9000.0
rows = []


def measure(p, types, x0, settle=80.0, window=40.0):
    sim = model.ParticleLife(p, N=len(types), seed=1)
    sim.types = np.asarray(types)
    sim.pos = np.asarray(x0, float).copy()
    sim.vel = np.zeros_like(sim.pos)
    sim.run(settle)
    X0, t0 = sim.pos.mean(0).copy(), sim.t
    sim.run(window)
    V = np.hypot(*(sim.pos.mean(0) - X0)) / (sim.t - t0)
    d = sim.pos[1] - sim.pos[0]
    return V, float(np.hypot(*d)), sim


start = np.array([[BOX / 2, BOX / 2], [BOX / 2 + 20.0, BOX / 2]])

# ---- (a) symmetric radii: the chase stalls -------------------------------- #
al, be = M.uniform_radii(2, 15.0, 45.0)
p = model.Params(A=M.chase_pair(0.9, -0.4), alpha=al, beta=be, L=(BOX, BOX))
pred = TB.dimer(p, 0, 1)
V, rm, _ = measure(p, [0, 1], start)
rows.append(dict(case="symmetric radii", alpha10=15.0, r_pred=pred["r_star"],
                 r_meas=rm, V_pred=abs(pred["V"]), V_meas=V))

# ---- (b) asymmetric min radii: a genuine self-propelled dimer -------------- #
for a10 in (8.0, 10.0, 12.0, 18.0, 20.0, 25.0, 30.0):
    al = np.array([[15.0, 15.0], [a10, 15.0]])
    be = np.full((2, 2), 45.0)
    p = model.Params(A=M.chase_pair(0.9, -0.4), alpha=al, beta=be, L=(BOX, BOX))
    pred = TB.dimer(p, 0, 1)
    V, rm, _ = measure(p, [0, 1], start)
    rows.append(dict(case="asymmetric radii", alpha10=a10, r_pred=pred["r_star"],
                     r_meas=rm, V_pred=abs(pred["V"]), V_meas=V))

print(f"{'case':18s} {'a10':>5s} {'r*_pred':>9s} {'r*_meas':>9s} "
      f"{'V_pred':>9s} {'V_meas':>9s} {'err %':>7s}")
for r in rows:
    e = 100 * abs(r["V_meas"] - r["V_pred"]) / max(r["V_pred"], 1e-9)
    print(f"{r['case']:18s} {r['alpha10']:5.0f} {r['r_pred']:9.4f} {r['r_meas']:9.4f} "
          f"{r['V_pred']:9.5f} {r['V_meas']:9.5f} {e:7.2f}")

# ---- (c) three-body frustration with *symmetric* radii -------------------- #
al3, be3 = M.uniform_radii(3, 15.0, 45.0)
p3 = model.Params(A=M.g_rps(3), alpha=al3, beta=be3, L=(BOX, BOX))
sim = model.ParticleLife(p3, N=3, seed=2)
sim.types = np.array([0, 1, 2])
sim.pos = np.array([[BOX / 2, BOX / 2], [BOX / 2 + 16, BOX / 2], [BOX / 2 + 8, BOX / 2 + 14]])
sim.vel = np.zeros((3, 2))
sim.run(60.0)
X0, t0 = sim.pos.mean(0).copy(), sim.t
ang0 = np.arctan2(*(sim.pos[1] - sim.pos[0])[::-1])
sim.run(30.0)
V3 = np.hypot(*(sim.pos.mean(0) - X0)) / (sim.t - t0)
ang1 = np.arctan2(*(sim.pos[1] - sim.pos[0])[::-1])
Fi = TB.cluster_internal_force(sim.pos, sim.types, p3)
Ti = TB.cluster_internal_torque(sim.pos, sim.types, p3)
V3_pred = p3.mobility_discrete * np.hypot(*Fi) / 3.0
omega3 = float(np.unwrap([ang0, ang1])[1] - ang0) / (sim.t - t0)
tri = dict(V_meas=V3, V_pred=V3_pred, F_int=float(np.hypot(*Fi)),
           tau_int=float(Ti), omega_meas=omega3)
print("\nRPS triangle, symmetric radii (three-body frustration):")
print(f"  |F_int| = {tri['F_int']:.5f}   V_pred = {V3_pred:.5f}   V_meas = {V3:.5f}")
print(f"  tau_int = {Ti:+.5f}            omega_meas = {omega3:+.5f} rad/s")

# ------------------------------------------------------------------- figures
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))

ax = axes[0]
rr = np.linspace(0.1, 50, 2000)
al = np.array([[15.0, 15.0], [25.0, 15.0]])
be = np.full((2, 2), 45.0)
pp = model.Params(A=M.chase_pair(0.9, -0.4), alpha=al, beta=be, L=(BOX, BOX))
ax.plot(rr, TB._pf(rr, pp, 0, 1), color=ACC[0], label=r"$f_{01}$  (0 chases 1)")
ax.plot(rr, TB._pf(rr, pp, 1, 0), color=ACC[1], label=r"$f_{10}$  (1 flees 0)")
ax.plot(rr, TB.f_sym(rr, pp, 0, 1), color=ACC[2], lw=2, label=r"$f_S$ (sets $r^*$)")
ax.plot(rr, TB.f_asym(rr, pp, 0, 1), color=ACC[3], lw=2, label=r"$f_A$ (sets $V$)")
rs = TB.bond_length_pair(pp, 0, 1)
ax.axvline(rs, color="w", ls="--", lw=0.9)
ax.axhline(0, color="#8b949e", lw=0.6)
ax.annotate(rf"$r^*$={rs:.2f},  $f_A(r^*)\neq0$", xy=(rs, TB.f_asym(np.array([rs]), pp, 0, 1)[0]),
            xytext=(rs + 4, 0.55), color="#e6edf3", fontsize=8,
            arrowprops=dict(arrowstyle="->", color="#8b949e", lw=0.8))
ax.set_xlabel("r"); ax.set_ylabel("force")
ax.set_title("decomposition of a chase pair", fontsize=10)
ax.legend(fontsize=7.5)

ax = axes[1]
asym = [r for r in rows if r["case"] == "asymmetric radii"]
ax.plot([r["V_pred"] for r in asym], [r["V_meas"] for r in asym], "o",
        color=ACC[0], ms=6, label="asymmetric radii")
sym = [r for r in rows if r["case"] == "symmetric radii"]
ax.plot([r["V_pred"] for r in sym], [r["V_meas"] for r in sym], "s",
        color=ACC[1], ms=6, label="symmetric radii (stalls)")
lim = max(0.05, max(r["V_meas"] for r in rows) * 1.15)
ax.plot([0, lim], [0, lim], color="#8b949e", ls="--", lw=1)
ax.set_xlim(0, lim); ax.set_ylim(0, lim)
ax.set_xlabel(r"predicted  $V=\mu\,f_A(r^*)$")
ax.set_ylabel("measured drift speed")
ax.set_title("self-propulsion: theory vs simulation", fontsize=10)
ax.legend(fontsize=8)

ax = axes[2]
a10 = np.linspace(6, 34, 160)
Vp = []
for a in a10:
    al = np.array([[15.0, 15.0], [a, 15.0]])
    pp = model.Params(A=M.chase_pair(0.9, -0.4), alpha=al, beta=np.full((2, 2), 45.0),
                      L=(BOX, BOX))
    Vp.append(abs(TB.dimer(pp, 0, 1)["V"]))
ax.plot(a10, Vp, color=ACC[2], lw=1.6, label="theory")
ax.plot([r["alpha10"] for r in asym], [r["V_meas"] for r in asym], "o",
        color=ACC[0], ms=5, label="simulation")
ax.axvline(15.0, color="#8b949e", ls="--", lw=0.9)
ax.text(15.2, max(Vp) * 0.9, r"$\alpha_{10}=\alpha_{01}$: $V=0$", color="#8b949e", fontsize=8)
ax.set_xlabel(r"$\alpha_{10}$   ($\alpha_{01}=15$ fixed)")
ax.set_ylabel("drift speed V")
ax.set_title("radius asymmetry drives the motor", fontsize=10)
ax.legend(fontsize=8)

fig.suptitle("E02  A non-reciprocal bound pair is a self-propelled particle", fontsize=12)
fig.tight_layout()
save(fig, "e02_self_propelled_dimer.png")
dump(dict(rows=rows, triangle=tri), "e02_dimer.json")
