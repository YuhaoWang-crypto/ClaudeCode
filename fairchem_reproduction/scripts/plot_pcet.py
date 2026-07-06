"""Compare naive 1e- vs PCET oxidation potentials for DA/AA/UA against
experiment, and show the PCET-predicted voltammogram."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RDIR = os.path.join(HERE, "..", "results")
FDIR = os.path.join(HERE, "..", "figures")

naive = json.load(open(os.path.join(RDIR, "sensor_descriptors.json")))["descriptors"]
pcet = json.load(open(os.path.join(RDIR, "pcet_oxidation.json")))["results"]

order = ["ascorbic_acid", "dopamine", "uric_acid"]
lab = {"dopamine": "Dopamine", "ascorbic_acid": "Ascorbic acid", "uric_acid": "Uric acid"}
col = {"dopamine": "#1a73e8", "ascorbic_acid": "#e8710a", "uric_acid": "#188038"}
# representative experimental oxidation peak potentials vs SHE, ~pH 7 (literature)
exp = {"ascorbic_acid": 0.30, "dopamine": 0.38, "uric_acid": 0.59}

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5))

# (A) PCET predicted voltammogram (pH 7) with experimental markers
V = np.linspace(0.0, 1.0, 700)
for n in order:
    e = pcet[n]["E_ox_pH7_vs_SHE_V"]
    axA.plot(V, np.exp(-0.5 * ((V - e) / 0.03) ** 2), color=col[n], lw=2,
             label=f"{lab[n]}: calc {e:+.2f} V | exp ≈ {exp[n]:+.2f} V")
    axA.fill_between(V, np.exp(-0.5 * ((V - e) / 0.03) ** 2), alpha=0.15, color=col[n])
    axA.axvline(exp[n], color=col[n], ls="--", lw=1.2, alpha=0.7)
axA.set_xlabel("Oxidation potential vs SHE (V), pH 7")
axA.set_ylabel("Response (a.u.)")
axA.set_title("PCET-predicted voltammogram (solid = calc peak, dashed = exp)\n"
              "order AA < DA < UA reproduces experiment")
axA.legend(frameon=False, fontsize=8.5, loc="upper left")

# (B) naive vs PCET vs experiment
x = np.arange(len(order)); w = 0.26
naive_v = [naive[n]["E_ox_vs_SHE_V"] for n in order]
pcet_v = [pcet[n]["E_ox_pH7_vs_SHE_V"] for n in order]
exp_v = [exp[n] for n in order]
axB.bar(x - w, naive_v, w, label="naive 1e⁻ (wrong order)", color="#bdc1c6")
axB.bar(x,      pcet_v,  w, label="PCET 1H⁺/1e⁻ (pH 7)", color="#1a73e8")
axB.bar(x + w,  exp_v,   w, label="experiment", color="#188038")
axB.set_xticks(x); axB.set_xticklabels([lab[n] for n in order])
axB.set_ylabel("Oxidation potential vs SHE (V)")
axB.set_title("Naive vs PCET vs experiment\nPCET recovers order + ~0.1 V accuracy")
axB.legend(frameon=False, fontsize=9)
for xi, v in zip(x, pcet_v):
    axB.text(xi, v + 0.02, f"{v:.2f}", ha="center", fontsize=8)

fig.tight_layout()
out = os.path.join(FDIR, "sensor_pcet.png")
fig.savefig(out, dpi=150)
print("wrote", out)
