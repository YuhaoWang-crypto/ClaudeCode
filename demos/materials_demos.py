"""
Two contrasting materials-science DEMOS runnable with only numpy/scipy/matplotlib
(no VASP/CP2K/phonopy/DeepMD, 4 cores).

DEMO A  -- RIGOROUS model reproduction.
    Reproduces the physics of PRL 106, 227401 (2011): exciton fine-structure
    splitting (FSS) and polarization angle of a quantum dot vs uniaxial stress.
    The bright-exciton 2x2 Hamiltonian is given explicitly in the paper (Eq. 2),
    so diagonalizing it is a faithful calculation, not a cartoon.

DEMO B  -- ILLUSTRATIVE workflow only (synthetic numbers).
    A single-atom-catalyst screening "volcano" of the kind used for
    Cr-Co3O4 / M-MnOC. The REAL version needs DFT adsorption energies (VASP/CP2K).
    Here the descriptors are invented to show the pipeline shape, NOT results.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# DEMO A : Exciton FSS & polarization angle under uniaxial stress
# Bright-state Hamiltonian  H = [[E+d+a3*p, k+b*p],[k+b*p, E-d+a4*p]]
# FSS(p)   = sqrt( (2d + (a3-a4) p)^2 + 4 (k + b p)^2 )
# angle(p): tan(2 theta) = 2(k+b p)/(2d + (a3-a4) p)
# Units: energies in ueV, stress p in MPa. Representative InAs/GaAs-scale params.
# ----------------------------------------------------------------------
d   = 15.0     # ueV   (delta: anisotropic exchange, zero-stress diagonal split)
k   = 8.0      # ueV   (kappa: off-diagonal coupling at zero stress -> C1 dot)
a3  = -0.55    # ueV/MPa
a4  =  0.10    # ueV/MPa
b   = -0.18    # ueV/MPa  (off-diagonal stress response)

p = np.linspace(0, 120, 600)                      # MPa
diag = 2*d + (a3 - a4)*p
offd = k + b*p
FSS  = np.sqrt(diag**2 + 4*offd**2)               # ueV
theta = 0.5*np.degrees(np.arctan2(2*offd, diag))  # deg

i_c = int(np.argmin(FSS))
p_c, FSS_min = p[i_c], FSS[i_c]

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(p, FSS, lw=2.2, color="#1f6feb")
ax[0].axvline(p_c, ls="--", color="#888")
ax[0].plot(p_c, FSS_min, "o", color="#d1242f", zorder=5)
ax[0].annotate(f"critical stress\n$p_c$={p_c:.0f} MPa\nFSS$_{{min}}$={FSS_min:.1f} $\\mu$eV",
               (p_c, FSS_min), textcoords="offset points", xytext=(20, 30),
               fontsize=9, arrowprops=dict(arrowstyle="->", color="#666"))
ax[0].set_xlabel("uniaxial stress $p$ (MPa)")
ax[0].set_ylabel("fine-structure splitting (μeV)")
ax[0].set_title("(A) Exciton FSS vs stress — faithful to PRL 106, 227401")

ax[1].plot(p, theta, lw=2.2, color="#8250df")
ax[1].axvline(p_c, ls="--", color="#888")
ax[1].set_xlabel("uniaxial stress $p$ (MPa)")
ax[1].set_ylabel("exciton polarization angle θ (deg)")
ax[1].set_title("(A) Polarization angle rotates through $p_c$")
fig.tight_layout()
fig.savefig("/home/user/ClaudeCode/demos/demoA_exciton_FSS.png", dpi=150)

# ----------------------------------------------------------------------
# DEMO B : ILLUSTRATIVE single-atom-catalyst volcano (synthetic descriptors)
# ----------------------------------------------------------------------
rng = np.random.default_rng(7)          # fixed seed (Date/random unavailable in workflows; fine in plain script)
metals = ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Ir"]
dGO = np.array([2.6,1.9,1.1,0.35,0.05,-0.3,-0.75,-1.2,-1.9,0.2])  # *O binding, fake eV
# toy 2-branch "volcano": activity limited by whichever step is uphill
overpotential = np.maximum(dGO - 0.0, -(dGO)) * 0.5 + 0.35 + rng.normal(0, 0.02, len(metals))

fig2, ax2 = plt.subplots(figsize=(6.4, 4.4))
ax2.scatter(dGO, overpotential, s=60, color="#2da44e", zorder=3)
for m, x, y in zip(metals, dGO, overpotential):
    ax2.annotate(m, (x, y), textcoords="offset points", xytext=(5, 4), fontsize=9)
xx = np.linspace(-2.1, 2.8, 200)
ax2.plot(xx, np.maximum(xx, -xx)*0.5 + 0.35, color="#999", ls="--", lw=1)
ax2.set_xlabel("ΔG$_{*O}$ descriptor (eV)  —  SYNTHETIC")
ax2.set_ylabel("theoretical overpotential η (V)  —  SYNTHETIC")
ax2.set_title("(B) SAC screening volcano — ILLUSTRATIVE ONLY\n(real ΔG needs DFT: VASP/CP2K, not run here)")
ax2.text(0.02, 0.02, "NOT PHYSICAL RESULTS — pipeline shape only",
         transform=ax2.transAxes, color="#d1242f", fontsize=8)
fig2.tight_layout()
fig2.savefig("/home/user/ClaudeCode/demos/demoB_SAC_volcano.png", dpi=150)

print(f"DEMO A  p_c = {p_c:.1f} MPa, FSS_min = {FSS_min:.2f} ueV, "
      f"FSS(0) = {FSS[0]:.2f} ueV")
print("Wrote demoA_exciton_FSS.png and demoB_SAC_volcano.png")
