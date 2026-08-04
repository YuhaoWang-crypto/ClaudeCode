"""
E01 - Anatomy of the interaction kernel.
========================================

Plots the force law, its potential, and the Fourier (Hankel) transform that the
continuum theory in ``stability.py`` is built on, and verifies the closed-form
results derived in ``kernels.py``:

  * ``U'(r) = f(r)``                              (potential is consistent)
  * ``U(alpha) = -A (beta-alpha)/2``              (well depth)
  * ``f(alpha) = 0``, stable for ``A > 0``        (bond length = min radius)
  * ``Uhat_core(0) = pi R alpha^3 / 12``          (core contribution at k=0)
"""

import numpy as np
from common import ACC, dump, save
import matplotlib.pyplot as plt

from plife import kernels as K

R = 1.0
checks = {}

# --------------------------------------------------------------------- checks
A, al, be = 0.7, 15.0, 45.0
r = np.linspace(1e-4, 50, 400001)
U = K.pair_potential(r, A, al, be, R)
f = K.pair_force(r, A, al, be, R)
m = (r > 0.3) & (r < be - 0.3) & (np.abs(r - al) > 0.2) & (np.abs(r - (al + be) / 2) > 0.2)
checks["max|dU/dr - f|"] = float(np.abs(np.gradient(U, r)[m] - f[m]).max())
checks["U(alpha)"] = float(K.pair_potential(al, A, al, be, R))
checks["-A(beta-alpha)/2"] = float(-A * (be - al) / 2)
checks["f(alpha)"] = float(K.pair_force(al, A, al, be, R))
checks["Uhat_core(0)"] = float(K.integrated_interaction(0.0, al, be, R))
checks["pi R alpha^3/12"] = float(np.pi * R * al**3 / 12)
for k, v in checks.items():
    print(f"  {k:24s} = {v:.6g}")

# ---------------------------------------------------------------------- plots
fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2))

ax = axes[0, 0]
rr = np.linspace(0, 60, 3000)
for A_, c in zip([1.0, 0.5, -0.5, -1.0], ACC):
    ax.plot(rr, K.pair_force(rr, A_, al, be, R), color=c, label=f"A = {A_:+.1f}")
ax.axvline(al, color="#8b949e", ls="--", lw=0.8)
ax.axvline(be, color="#8b949e", ls=":", lw=0.8)
ax.axhline(0, color="#8b949e", lw=0.6)
ax.text(al, 1.02, r" $\alpha$ = min radius", color="#8b949e", fontsize=8, va="bottom")
ax.text(be, 1.02, r" $\beta$ = max radius", color="#8b949e", fontsize=8, va="bottom")
ax.set_xlabel("r"); ax.set_ylabel("f(r)   (+ = attraction)")
ax.set_title("force law: repulsive core + attraction tent", fontsize=10)
ax.legend(fontsize=8)

ax = axes[0, 1]
for A_, c in zip([1.0, 0.5, -0.5, -1.0], ACC):
    ax.plot(rr, K.pair_potential(rr, A_, al, be, R), color=c)
ax.axvline(al, color="#8b949e", ls="--", lw=0.8)
ax.axhline(0, color="#8b949e", lw=0.6)
ax.plot([al], [K.pair_potential(al, 1.0, al, be, R)], "o", color="w", ms=4, zorder=5)
ax.annotate(r"bond: $r^*=\alpha$,  $U=-A(\beta-\alpha)/2$",
            xy=(al, K.pair_potential(al, 1.0, al, be, R)), xytext=(24, -12),
            color="#e6edf3", fontsize=8,
            arrowprops=dict(arrowstyle="->", color="#8b949e", lw=0.8))
ax.set_xlabel("r"); ax.set_ylabel("U(r)")
ax.set_title(r"pair potential  ($U'=f$, piecewise quadratic)", fontsize=10)

ax = axes[1, 0]
kk = np.linspace(1e-3, 0.6, 700)
for A_, c in zip([0.9, 0.5, 0.2, -0.4], ACC):
    ax.plot(kk, K.hankel_transform(kk, A_, al, be, R), color=c, label=f"A = {A_:+.1f}")
ax.axhline(0, color="#8b949e", lw=0.6)
ax.set_xlabel("k"); ax.set_ylabel(r"$\hat U(k)$")
ax.set_title(r"Hankel transform: $\hat U(k)<0 \Rightarrow$ that mode is unstable", fontsize=10)
ax.legend(fontsize=8)

ax = axes[1, 1]
alphas = np.linspace(4, 42, 60)
for A_, c in zip([0.9, 0.5, 0.3, 0.2], ACC):
    u0 = [K.integrated_interaction(A_, a, be, R) for a in alphas]
    ax.plot(alphas, u0, color=c, label=f"A = {A_:.1f}")
ax.plot(alphas, np.pi * R * alphas**3 / 12, color="#8b949e", ls="--", lw=1,
        label=r"core only $\pi R\alpha^3/12$")
ax.axhline(0, color="#8b949e", lw=0.6)
ax.set_xlabel(r"min radius $\alpha$   ($\beta$ = 45 fixed)")
ax.set_ylabel(r"$\hat U(0)$")
ax.set_title(r"$\hat U(0)$ changes sign: cohesion $\to$ dispersal", fontsize=10)
ax.legend(fontsize=8)

fig.suptitle("E01  Anatomy of the Particle Life interaction kernel", fontsize=12)
fig.tight_layout()
save(fig, "e01_force_law.png")
dump(checks, "e01_checks.json")
