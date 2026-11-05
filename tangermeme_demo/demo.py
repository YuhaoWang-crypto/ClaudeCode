"""tangermeme end-to-end demo on a synthetic task with a known answer.

We write the "regulatory grammar" ourselves, train a small CNN on it, and then
use tangermeme to recover that grammar from the trained model. Because the
ground truth is known, every analysis can be checked.

Ground-truth grammar (per 200 bp sequence):
    GATA  "GATAAG"   activator  +1.0 per copy
    EBOX  "CACGTG"   activator  +1.5 per copy
    AP1   "TGACTCA"  repressor  -1.0 per copy
    GATA -> EBOX cooperation: +2.0 * exp(-(gap - 10)^2 / 8)
        gap = bases between the end of GATA and the start of EBOX (GATA upstream)

Run:  python3 tangermeme_demo/demo.py
"""

import os
import re
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from torch import nn

from tangermeme.utils import random_one_hot, one_hot_encode, characters
from tangermeme.ersatz import substitute
from tangermeme.predict import predict
from tangermeme.marginalize import marginalize
from tangermeme.space import space
from tangermeme.deep_lift_shap import deep_lift_shap
from tangermeme.saturation_mutagenesis import saturation_mutagenesis
from tangermeme.seqlet import recursive_seqlets
from tangermeme.annotate import annotate_seqlets
from tangermeme.variant_effect import substitution_effect
from tangermeme.design import greedy_substitution
from tangermeme.plot import plot_logo

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

torch.manual_seed(0)
torch.set_num_threads(os.cpu_count() or 1)
rng = np.random.default_rng(0)

L = 200
MOTIFS = {"GATA": "GATAAG", "EBOX": "CACGTG", "AP1": "TGACTCA"}
WEIGHTS = {"GATA": 1.0, "EBOX": 1.5, "AP1": -1.0}
BG_PROBS = np.array([[0.3, 0.2, 0.2, 0.3]])  # AT-rich, genome-like background
OPT_GAP, COOP = 10, 2.0


def coop_bonus(gap):
    return COOP * np.exp(-((gap - OPT_GAP) ** 2) / 8.0)


def true_activity(seq):
    """Score a sequence with the ground-truth grammar (scans the actual string,
    so accidental motif occurrences in the background are counted too)."""
    hits = {k: [m.start() for m in re.finditer(f"(?={v})", seq)] for k, v in MOTIFS.items()}
    y = sum(WEIGHTS[k] * len(h) for k, h in hits.items())
    best = 0.0
    for g in hits["GATA"]:
        for e in hits["EBOX"]:
            gap = e - (g + len(MOTIFS["GATA"]))
            if gap >= 0:
                best = max(best, coop_bonus(gap))
    return y + best


def make_dataset(n):
    X = random_one_hot((n, 4, L), probs=BG_PROBS, random_state=int(rng.integers(1e9))).numpy()
    seqs = []
    for i in range(n):
        s = list(characters(torch.from_numpy(X[i])))
        used = np.zeros(L, bool)

        def place(motif, pos=None):
            if pos is None:
                for _ in range(50):
                    p = int(rng.integers(0, L - len(motif)))
                    if not used[p:p + len(motif)].any():
                        pos = p
                        break
                else:
                    return None
            s[pos:pos + len(motif)] = motif
            used[pos:pos + len(motif)] = True
            return pos

        if rng.random() < 0.5:  # a GATA->EBOX pair at a random gap (0..30)
            gap = int(rng.integers(0, 31))
            g = int(rng.integers(0, L - 6 - gap - 6))
            place(MOTIFS["GATA"], g)
            place(MOTIFS["EBOX"], g + 6 + gap)
        for k, v in MOTIFS.items():
            for _ in range(rng.binomial(2, 0.3)):
                place(v)
        seqs.append("".join(s))
    X = torch.stack([one_hot_encode(s) for s in seqs])
    y = torch.tensor([true_activity(s) for s in seqs], dtype=torch.float32)
    y = y + 0.1 * torch.randn_like(y)
    return X, y[:, None], seqs


class CNN(nn.Module):
    """Conv -> dilated convs -> mean over length -> linear. Only ops DeepLIFT
    has rules for (Conv1d, ReLU, mean, Linear)."""

    def __init__(self, h=64):
        super().__init__()
        self.c1 = nn.Conv1d(4, h, 11, padding=5)
        self.c2 = nn.Conv1d(h, h, 9, padding=4)
        self.c3 = nn.Conv1d(h, h, 9, padding=8, dilation=2)
        self.r1, self.r2, self.r3 = nn.ReLU(), nn.ReLU(), nn.ReLU()
        self.fc = nn.Linear(h, 1)

    def forward(self, X):
        X = X.float()
        h = self.r1(self.c1(X))
        h = self.r2(self.c2(h)) + h
        h = self.r3(self.c3(h)) + h
        return self.fc(h.mean(dim=-1)) * L / 50.0


def section(t):
    print("\n" + "=" * 70 + f"\n{t}\n" + "=" * 70)


t0 = time.time()

# ---------------------------------------------------------------- 1. data + train
section("1. Synthetic data + train a CNN")
X_tr, y_tr, _ = make_dataset(30000)
X_te, y_te, seqs_te = make_dataset(2000)
print(f"train {tuple(X_tr.shape)}, test {tuple(X_te.shape)}, y range [{y_tr.min():.2f}, {y_tr.max():.2f}]")

model = CNN()
EPOCHS = int(os.environ.get("EPOCHS", 16))
CKPT = os.environ.get("CKPT")
opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=1e-4)
sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=3e-3, total_steps=EPOCHS * (len(X_tr) // 128 + 1))
if CKPT and os.path.exists(CKPT):
    model.load_state_dict(torch.load(CKPT))
    EPOCHS = 0
for epoch in range(EPOCHS):
    model.train()
    perm = torch.randperm(len(X_tr))
    for b in range(0, len(X_tr), 128):
        idx = perm[b:b + 128]
        loss = nn.functional.mse_loss(model(X_tr[idx]), y_tr[idx])
        opt.zero_grad()
        loss.backward()
        opt.step()
        sched.step()
    model.eval()
    # tangermeme.predict: batched inference
    y_hat = predict(model, X_te, batch_size=512, device="cpu")
    r = np.corrcoef(y_hat[:, 0].numpy(), y_te[:, 0].numpy())[0, 1]
    print(f"epoch {epoch + 1:2d}  train mse {loss.item():.3f}  test pearson r {r:.3f}")
model.eval()
if CKPT and not os.path.exists(CKPT):
    torch.save(model.state_dict(), CKPT)
y_hat = predict(model, X_te, batch_size=512, device="cpu")
r = np.corrcoef(y_hat[:, 0].numpy(), y_te[:, 0].numpy())[0, 1]

fig, ax = plt.subplots(figsize=(4, 4))
ax.scatter(y_te[:, 0], y_hat[:, 0], s=4, alpha=0.4)
lim = [y_te.min().item() - 0.5, y_te.max().item() + 0.5]
ax.plot(lim, lim, "k--", lw=0.8)
ax.set(xlabel="true activity", ylabel="predicted (tangermeme.predict)", title=f"held-out test, r = {r:.3f}")
fig.tight_layout()
fig.savefig(f"{FIG}/1_predict.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------- 2. marginalize
section("2. marginalize: insert each motif into neutral backgrounds")
X_bg = random_one_hot((200, 4, L), probs=BG_PROBS, random_state=1)
rows = []
for name, motif in {**MOTIFS, "control(shuffled)": "AGTCAG"}.items():
    y0, y1 = marginalize(model, X_bg, motif, device="cpu")
    d = (y1 - y0)[:, 0]
    rows.append((name, motif, WEIGHTS.get(name, 0.0), d.mean().item(), d.std().item()))
    print(f"{name:18s} {motif:8s} true {WEIGHTS.get(name, 0.0):+.2f}   recovered {d.mean():+.3f} ± {d.std():.3f}")
marg = pd.DataFrame(rows, columns=["name", "motif", "true", "effect", "sd"])

fig, ax = plt.subplots(figsize=(5, 3.2))
xx = np.arange(len(marg))
ax.bar(xx - 0.2, marg["true"], 0.4, label="ground truth", color="#bbb")
ax.bar(xx + 0.2, marg["effect"], 0.4, yerr=marg["sd"], label="tangermeme.marginalize", color="#3b6fb6")
ax.set_xticks(xx, marg["name"], fontsize=8)
ax.axhline(0, c="k", lw=0.6)
ax.set_ylabel("Δ prediction")
ax.legend(fontsize=8)
ax.set_title("motif effect in isolation")
fig.tight_layout()
fig.savefig(f"{FIG}/2_marginalize.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------- 3. space
section("3. space: GATA -> EBOX spacing sweep (cooperativity)")
gaps = torch.arange(0, 41)[:, None]
y0, y_sp = space(model, X_bg[:100], [MOTIFS["GATA"], MOTIFS["EBOX"]], spacing=gaps, start=60, device="cpu")
# Subtract the two single-motif effects to isolate the interaction term.
_, yg = marginalize(model, X_bg[:100], MOTIFS["GATA"], start=60, device="cpu")
single_e = []
for gp in gaps[:, 0].tolist():
    _, ye = marginalize(model, X_bg[:100], MOTIFS["EBOX"], start=60 + 6 + gp, device="cpu")
    single_e.append((ye - y0)[:, 0].mean().item())
pair = (y_sp[:, :, 0] - y0[:, None, 0]).mean(0).numpy()
interaction = pair - (yg - y0)[:, 0].mean().item() - np.array(single_e)
truth = np.array([coop_bonus(g) for g in gaps[:, 0].numpy()])
best = int(gaps[np.argmax(interaction), 0])
print(f"recovered optimal gap = {best} bp (truth {OPT_GAP}); peak interaction {interaction.max():.2f} (truth {COOP})")

fig, ax = plt.subplots(figsize=(5.5, 3.2))
ax.plot(gaps[:, 0], truth, "k--", label="ground-truth bonus")
ax.plot(gaps[:, 0], interaction, "o-", ms=3, c="#c0392b", label="tangermeme.space (interaction)")
ax.set(xlabel="bases between GATA and EBOX", ylabel="Δ beyond additive", title=f"cooperative spacing; best gap = {best} bp")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{FIG}/3_space.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------- 4. DeepLIFT/SHAP + seqlets
section("4. deep_lift_shap attributions + recursive_seqlets + TOMTOM annotation")
X_a = X_te[:500]
X_attr = deep_lift_shap(model, X_a, n_shuffles=20, random_state=0, device="cpu", batch_size=256)
attr_sums = X_attr.sum(dim=1)

# Positive seqlets (activators) and negative seqlets (repressors) separately.
pos = recursive_seqlets(attr_sums, threshold=0.01)
neg = recursive_seqlets(-attr_sums, threshold=0.01)
neg["attribution"] = -neg["attribution"]
seqlets = pd.concat([pos, neg], ignore_index=True)
print(f"called {len(pos)} positive and {len(neg)} negative seqlets on {len(X_a)} sequences")

pwms = {k: one_hot_encode(v).float() * 0.97 + 0.0075 for k, v in MOTIFS.items()}
names = list(pwms)
motif_idxs, pvals = annotate_seqlets(X_a, seqlets, pwms, n_nearest=1, n_jobs=4)
seqlets.insert(0, "motif_name", [names[j] if j >= 0 else "none" for j in motif_idxs[:, 0].numpy()])
seqlets["tomtom_p"] = pvals[:, 0].numpy()
print(seqlets.groupby("motif_name").attribution.agg(["count", "mean"]).round(3))

# Recall against the planted / occurring motif instances.
truth_sites = []
for i, s in enumerate(seqs_te[:500]):
    for k, v in MOTIFS.items():
        for m in re.finditer(f"(?={v})", s):
            truth_sites.append((i, k, m.start(), m.start() + len(v)))
found = 0
for i, k, st, en in truth_sites:
    sub = seqlets[(seqlets.example_idx == i) & (seqlets.start < en) & (seqlets.end > st)]
    found += int((sub.motif_name == k).any())
precision = (seqlets.motif_name != "none").mean()
print(f"recall of true motif sites: {found}/{len(truth_sites)} = {found / len(truth_sites):.1%}")

# Pick an example with all three motifs for the logo figure.
def has_all(s):
    return all(v in s for v in MOTIFS.values())
ex = next(i for i, s in enumerate(seqs_te[:500]) if has_all(s))
fig, axes = plt.subplots(2, 1, figsize=(13, 4.8))
s_ex = seqlets[seqlets.example_idx == ex].sort_values("start").copy()
plot_logo(X_attr[ex], ax=axes[0], annotations=s_ex, score_key="attribution")
axes[0].set_title(f"DeepLIFT/SHAP attribution, test seq #{ex}  (true: " +
                  ", ".join(f"{k}@{m.start()}" for k, v in MOTIFS.items()
                            for m in re.finditer(f"(?={v})", seqs_te[ex])) + ")")

# ---------------------------------------------------------------- 5. ISM
section("5. saturation_mutagenesis (ISM) on the same sequence")
X_ism = saturation_mutagenesis(model, X_te[ex:ex + 1], device="cpu", target=0)
plot_logo(X_ism[0], ax=axes[1])
axes[1].set_title("in-silico saturation mutagenesis (forward passes only)")
fig.tight_layout()
fig.savefig(f"{FIG}/4_attribution_ism.png", dpi=150)
plt.close(fig)
a = X_attr[ex].sum(0).numpy()
b = X_ism[0].sum(0).numpy()
print(f"Pearson(DeepLIFT, ISM) on seq #{ex}: {np.corrcoef(a, b)[0, 1]:.3f}")

# ---------------------------------------------------------------- 6. variant effect
section("6. substitution_effect: point mutations inside the GATA site")
g = seqs_te[ex].find(MOTIFS["GATA"])
idx = {c: i for i, c in enumerate("ACGT")}
var_rows, labels = [], []
for off in range(6):
    ref = seqs_te[ex][g + off]
    for alt in "ACGT":
        if alt != ref:
            var_rows.append([len(var_rows), g + off, idx[alt]])
            labels.append(f"{g + off}:{ref}>{alt}")
subs = torch.tensor(var_rows)
y0, y1 = substitution_effect(model, X_te[ex:ex + 1].repeat(len(subs), 1, 1), subs, device="cpu")
veff = (y1 - y0)[:, 0].numpy()
order = np.argsort(veff)
for j in order[:5]:
    print(f"  {labels[j]:10s} Δ = {veff[j]:+.3f}")
print(f"mean Δ over all 18 SNVs in GATA: {veff.mean():+.3f}")

# ---------------------------------------------------------------- 7. design
section("7. greedy_substitution: design a high-activity sequence")
X0 = random_one_hot((1, 4, L), probs=BG_PROBS, random_state=7)
y_start = predict(model, X0, device="cpu")[0, 0].item()
target = torch.tensor([[12.0]])
X_bar = greedy_substitution(model, X0, y=target, motifs=list(MOTIFS.values()), max_iter=6,
                            reverse_complement=False, device="cpu")
y_end = predict(model, X_bar, device="cpu")[0, 0].item()
s_bar = characters(X_bar[0])
placed = sorted((m.start(), k) for k, v in MOTIFS.items() for m in re.finditer(f"(?={v})", s_bar))
print(f"start pred {y_start:.2f} -> designed pred {y_end:.2f}  (true grammar score {true_activity(s_bar):.2f})")
print("motifs in designed sequence:", placed)
gata = [p for p, k in placed if k == "GATA"]
ebox = [p for p, k in placed if k == "EBOX"]
gap_pairs = [e - (gg + 6) for gg in gata for e in ebox if e - (gg + 6) >= 0]
print("GATA->EBOX gaps in design:", gap_pairs)

X_bar_attr = deep_lift_shap(model, X_bar, n_shuffles=20, random_state=0, device="cpu")
fig, ax = plt.subplots(figsize=(13, 2.4))
plot_logo(X_bar_attr[0], ax=ax)
ax.set_title(f"designed sequence (pred {y_start:.1f} -> {y_end:.1f}); motifs: " +
             ", ".join(f"{k}@{p}" for p, k in placed))
fig.tight_layout()
fig.savefig(f"{FIG}/5_design.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------- summary
section("Summary")
summary = {
    "test_pearson_r": round(float(r), 3),
    "marginalize": {n: round(e, 3) for n, e in zip(marg["name"], marg["effect"])},
    "space_best_gap": best,
    "space_peak_interaction": round(float(interaction.max()), 3),
    "seqlet_recall": round(found / len(truth_sites), 3),
    "seqlet_labelled_frac": round(float(precision), 3),
    "deeplift_vs_ism_r": round(float(np.corrcoef(a, b)[0, 1]), 3),
    "design_pred_start_end": (round(y_start, 2), round(y_end, 2)),
    "design_true_score": round(true_activity(s_bar), 2),
    "design_gata_ebox_gaps": gap_pairs,
    "runtime_s": round(time.time() - t0, 1),
}
for k, v in summary.items():
    print(f"{k:24s} {v}")
pd.Series(summary).to_json(os.path.join(HERE, "results.json"), indent=2)
