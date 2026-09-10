"""Neural components of Gen-COMPAS.

1. `DiffusionModel` -- a denoising diffusion probabilistic model (Ho et al.,
   NeurIPS 2020) over rotation/translation-invariant heavy-atom coordinates,
   with deterministic DDIM (Song et al., ICLR 2021) inversion and sampling so
   that *intermediates* can be produced by interpolating two end-state
   structures in the latent space of the model.

2. `CommittorNet` -- a network that predicts the committor q(x) directly in
   conformational space, i.e. without any predefined collective variable.

Both act on the same feature vector, so no chemical intuition is injected.
"""

import numpy as np
import torch
import torch.nn as nn

torch.set_num_threads(1)


# ---------------------------------------------------------------------------
# Diffusion model
# ---------------------------------------------------------------------------


def timestep_embedding(t, dim):
    half = dim // 2
    freqs = torch.exp(
        -np.log(10000.0) * torch.arange(half, dtype=torch.float32) / half
    ).to(t.device)
    args = t.float()[:, None] * freqs[None]
    return torch.cat([torch.cos(args), torch.sin(args)], dim=-1)


class EpsNet(nn.Module):
    """Noise-prediction network eps_theta(x_t, t)."""

    def __init__(self, dim, hidden=512, tdim=128, n_blocks=4):
        super().__init__()
        self.tdim = tdim
        self.tmlp = nn.Sequential(
            nn.Linear(tdim, hidden), nn.SiLU(), nn.Linear(hidden, hidden)
        )
        self.inp = nn.Linear(dim, hidden)
        self.blocks = nn.ModuleList(
            [
                nn.Sequential(
                    nn.LayerNorm(hidden),
                    nn.Linear(hidden, hidden),
                    nn.SiLU(),
                    nn.Linear(hidden, hidden),
                )
                for _ in range(n_blocks)
            ]
        )
        self.out = nn.Sequential(nn.LayerNorm(hidden), nn.Linear(hidden, dim))

    def forward(self, x, t):
        h = self.inp(x) + self.tmlp(timestep_embedding(t, self.tdim))
        for blk in self.blocks:
            h = h + blk(h)
        return self.out(h)


class DiffusionModel:
    def __init__(self, dim, n_steps=400, hidden=512, seed=0):
        torch.manual_seed(seed)
        self.dim = dim
        self.T = n_steps
        self.net = EpsNet(dim, hidden=hidden)
        # cosine schedule
        s = 0.008
        ts = torch.arange(n_steps + 1, dtype=torch.float64) / n_steps
        f = torch.cos((ts + s) / (1 + s) * np.pi / 2) ** 2
        ab = (f / f[0]).float()
        self.alpha_bar = ab[1:].clamp(1e-5, 0.9999)
        self.mu = None
        self.sd = None

    # -- data scaling -------------------------------------------------------
    def set_scaler(self, x):
        self.mu = x.mean(0)
        self.sd = x.std(0) + 1e-6

    def _enc(self, x):
        return (x - self.mu) / self.sd

    def _dec(self, z):
        return z * self.sd + self.mu

    # -- training -----------------------------------------------------------
    def fit(self, x, steps=8000, batch=256, lr=2e-4, verbose=False,
            max_samples=6000, seed=0):
        """Train for a fixed budget of gradient steps, so that the cost per
        Gen-COMPAS iteration does not grow with the accumulated dataset."""
        if len(x) > max_samples:
            sel = np.random.default_rng(seed).choice(len(x), max_samples,
                                                     replace=False)
            x = x[sel]
        if self.mu is None:
            self.set_scaler(x)
        per_epoch = max(1, int(np.ceil(len(x) / batch)))
        epochs = max(50, int(np.ceil(steps / per_epoch)))
        z = torch.tensor(self._enc(x), dtype=torch.float32)
        opt = torch.optim.AdamW(self.net.parameters(), lr=lr, weight_decay=1e-6)
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)
        n = len(z)
        self.net.train()
        for ep in range(epochs):
            perm = torch.randperm(n)
            tot = 0.0
            for i in range(0, n, batch):
                idx = perm[i:i + batch]
                x0 = z[idx]
                t = torch.randint(0, self.T, (len(idx),))
                ab = self.alpha_bar[t][:, None]
                eps = torch.randn_like(x0)
                xt = ab.sqrt() * x0 + (1 - ab).sqrt() * eps
                loss = ((self.net(xt, t) - eps) ** 2).mean()
                opt.zero_grad()
                loss.backward()
                opt.step()
                tot += loss.item() * len(idx)
            sched.step()
            if verbose and (ep + 1) % max(1, epochs // 5) == 0:
                print(f"    [ddpm] epoch {ep + 1}/{epochs} "
                      f"({(ep + 1) * per_epoch} steps) loss={tot / n:.4f}",
                      flush=True)
        self.net.eval()
        return self

    # -- deterministic DDIM -------------------------------------------------
    def _ddim_steps(self, n_infer):
        return np.linspace(0, self.T - 1, n_infer).round().astype(int)

    @torch.no_grad()
    def invert(self, x, n_infer=50):
        """Map data -> latent (deterministic DDIM inversion)."""
        z = torch.tensor(self._enc(x), dtype=torch.float32)
        steps = self._ddim_steps(n_infer)
        for i in range(len(steps) - 1):
            t, t_next = steps[i], steps[i + 1]
            ab, ab_n = self.alpha_bar[t], self.alpha_bar[t_next]
            tt = torch.full((len(z),), int(t), dtype=torch.long)
            eps = self.net(z, tt)
            x0 = (z - (1 - ab).sqrt() * eps) / ab.sqrt()
            z = ab_n.sqrt() * x0 + (1 - ab_n).sqrt() * eps
        return z

    @torch.no_grad()
    def sample(self, z, n_infer=50):
        """Map latent -> data (deterministic DDIM sampling)."""
        steps = self._ddim_steps(n_infer)[::-1]
        for i in range(len(steps) - 1):
            t, t_prev = steps[i], steps[i + 1]
            ab, ab_p = self.alpha_bar[t], self.alpha_bar[t_prev]
            tt = torch.full((len(z),), int(t), dtype=torch.long)
            eps = self.net(z, tt)
            x0 = (z - (1 - ab).sqrt() * eps) / ab.sqrt()
            z = ab_p.sqrt() * x0 + (1 - ab_p).sqrt() * eps
        tt = torch.zeros(len(z), dtype=torch.long)
        eps = self.net(z, tt)
        ab0 = self.alpha_bar[0]
        x0 = (z - (1 - ab0).sqrt() * eps) / ab0.sqrt()
        return self._dec(x0.numpy())

    def interpolate(self, xa, xb, lambdas, n_infer=50):
        """Generate intermediates between paired end-state structures.

        xa, xb : (m, dim) arrays of end-state conformations
        lambdas: sequence of interpolation weights in (0, 1)
        returns (m * len(lambdas), dim)
        """
        za = self.invert(xa, n_infer)
        zb = self.invert(xb, n_infer)
        out = []
        for lam in lambdas:
            z = _slerp(za, zb, float(lam))
            out.append(self.sample(z, n_infer))
        return np.concatenate(out, axis=0)


def _slerp(a, b, lam):
    """Spherical interpolation, the correct geometry for gaussian latents."""
    an = a / a.norm(dim=-1, keepdim=True)
    bn = b / b.norm(dim=-1, keepdim=True)
    dot = (an * bn).sum(-1, keepdim=True).clamp(-0.9999, 0.9999)
    omega = torch.acos(dot)
    so = torch.sin(omega)
    w = (torch.sin((1 - lam) * omega) / so) * a + (torch.sin(lam * omega) / so) * b
    return w


# ---------------------------------------------------------------------------
# Committor network
# ---------------------------------------------------------------------------


class CommittorNet(nn.Module):
    def __init__(self, dim, hidden=256, n_blocks=3):
        super().__init__()
        layers = [nn.Linear(dim, hidden), nn.SiLU()]
        for _ in range(n_blocks):
            layers += [nn.Linear(hidden, hidden), nn.SiLU()]
        layers += [nn.Linear(hidden, 1)]
        self.net = nn.Sequential(*layers)
        self.mu = None
        self.sd = None

    def set_scaler(self, x):
        self.mu = torch.tensor(x.mean(0), dtype=torch.float32)
        self.sd = torch.tensor(x.std(0) + 1e-6, dtype=torch.float32)

    def logits(self, x):
        return self.net((x - self.mu) / self.sd).squeeze(-1)

    def forward(self, x):
        return torch.sigmoid(self.logits(x))

    @torch.no_grad()
    def predict(self, x):
        self.eval()
        t = torch.tensor(np.asarray(x), dtype=torch.float32)
        return torch.sigmoid(self.logits(t)).numpy()


def train_committor(x, y, w=None, dim=None, steps=4000, lr=1e-3, batch=512,
                    hidden=256, n_blocks=3, weight_decay=1e-5, seed=0,
                    verbose=False):
    """Fit q(x) = P(reach B before A | x) by weighted logistic regression on
    observed commitment outcomes.  This is the definition of the committor, so
    the fit is a Monte-Carlo estimator of it -- no CV, no biasing potential."""
    torch.manual_seed(seed)
    dim = dim or x.shape[1]
    model = CommittorNet(dim, hidden=hidden, n_blocks=n_blocks)
    model.set_scaler(x)
    xt = torch.tensor(x, dtype=torch.float32)
    yt = torch.tensor(y, dtype=torch.float32)
    wt = torch.ones(len(x)) if w is None else torch.tensor(w, dtype=torch.float32)
    per_epoch = max(1, int(np.ceil(len(xt) / batch)))
    epochs = max(50, int(np.ceil(steps / per_epoch)))
    opt = torch.optim.AdamW(model.parameters(), lr=lr,
                            weight_decay=weight_decay)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)
    bce = nn.BCEWithLogitsLoss(reduction="none")
    n = len(xt)
    model.train()
    for ep in range(epochs):
        perm = torch.randperm(n)
        tot = 0.0
        for i in range(0, n, batch):
            idx = perm[i:i + batch]
            loss = (bce(model.logits(xt[idx]), yt[idx]) * wt[idx]).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += loss.item() * len(idx)
        sched.step()
        if verbose and (ep + 1) % max(1, epochs // 5) == 0:
            print(f"    [committor] epoch {ep + 1}/{epochs} loss={tot / n:.4f}",
                  flush=True)
    model.eval()
    return model
