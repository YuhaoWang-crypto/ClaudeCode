"""Visualise the DA/AA/UA sensor-descriptor fingerprint:
 (A) simulated DPV -- Gaussian peaks at each predicted oxidation potential,
     showing whether a sensor could resolve them;
 (B) normalised descriptor fingerprint heatmap across the panel."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RDIR = os.path.join(HERE, "..", "results")
FDIR = os.path.join(HERE, "..", "figures")
os.makedirs(FDIR, exist_ok=True)

d = json.load(open(os.path.join(RDIR, "sensor_descriptors.json")))["descriptors"]
names = list(d)
labels = {"dopamine": "Dopamine", "ascorbic_acid": "Ascorbic acid",
          "uric_acid": "Uric acid"}
colors = {"dopamine": "#1a73e8", "ascorbic_acid": "#e8710a", "uric_acid": "#188038"}

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5))

# (A) simulated DPV
V = np.linspace(0.5, 2.0, 600)
for n in names:
    e = d[n]["E_ox_vs_SHE_V"]
    peak = np.exp(-0.5 * ((V - e) / 0.05) ** 2)
    axA.plot(V, peak, color=colors[n], lw=2, label=f"{labels[n]}  ({e:+.2f} V)")
    axA.fill_between(V, peak, alpha=0.15, color=colors[n])
    axA.axvline(e, color=colors[n], ls=":", lw=1)
axA.set_xlabel("Predicted oxidation potential vs SHE (V)")
axA.set_ylabel("Response (a.u.)")
axA.set_title("Predicted voltammetric peaks\n(1e⁻ oxidation, ddCOSMO water)")
axA.legend(frameon=False, fontsize=9)

# (B) descriptor fingerprint heatmap (z-scored per descriptor)
keys = ["HOMO_eV", "LUMO_eV", "gap_eV", "eta_eV", "omega_eV", "dipole_D",
        "IP_gas_eV", "E_ox_vs_SHE_V"]
klabel = ["HOMO", "LUMO", "gap", "η", "ω", "μ", "IP", "E_ox"]
M = np.array([[d[n][k] for k in keys] for n in names], float)
Z = (M - M.mean(0)) / (M.std(0) + 1e-9)
im = axB.imshow(Z, aspect="auto", cmap="coolwarm", vmin=-1.5, vmax=1.5)
axB.set_xticks(range(len(keys))); axB.set_xticklabels(klabel)
axB.set_yticks(range(len(names))); axB.set_yticklabels([labels[n] for n in names])
for i in range(len(names)):
    for j in range(len(keys)):
        axB.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=8)
axB.set_title("Descriptor fingerprint (cell = raw value; colour = z-score)")
fig.colorbar(im, ax=axB, shrink=0.7, label="z-score across panel")

fig.tight_layout()
out = os.path.join(FDIR, "sensor_fingerprint.png")
fig.savefig(out, dpi=150)
print("wrote", out)
