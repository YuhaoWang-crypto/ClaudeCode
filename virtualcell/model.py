"""A virtual cell model for zero-shot, cross-context knockdown prediction.

Task shape, following Virtual Cell Challenge 2: a cell line is presented with
*only* its non-targeting control profiles, together with a list of genes to
knock down.  The model must produce post-perturbation expression profiles
without ever having seen a perturbation measured in that cell line.

Design follows what actually won the 2025 challenge rather than what is
fashionable: every entry in the top three combined transferred summary-level
statistics with explicit context conditioning, and Arc's own write-up concluded
that purely learned approaches did not reliably beat statistical baselines.  So
this is a structured statistical model whose every component is switchable, with
the switches chosen by nested leave-one-cell-line-out cross-validation over the
*source* lines only.  The target line contributes nothing but its controls.

The prediction is assembled in seven stages:

1. **Context-weighted consensus.**  Average the effect of knockdown ``p`` over
   source lines, weighting each by how much its control transcriptome resembles
   the target's on the genes that actually vary between cell lines.
2. **Perturbation-similarity smoothing.**  Borrow strength from knockdowns with
   correlated effects, which denoises the weak ones.
3. **Reliability shrinkage.**  A James-Stein style factor per perturbation:
   effects that replicate across source lines survive, effects that do not are
   pulled towards the generic response.
4. **Program-basis denoising.**  Project onto the leading directions of the
   source response matrix, blended rather than applied outright.
5. **Context modulation.**  Rescale each gene by how its baseline expression in
   the target compares to the sources.  A gene that is silent in the target
   cannot be repressed there, however reliably it moves elsewhere.
6. **On-target knockdown.**  CRISPRi represses the targeted transcript itself,
   applied in count space where a knockdown is multiplicative, at the per-gene
   efficiency measured in the source lines.
7. **Global calibration.**  One scalar trading effect magnitude against error,
   fitted on source lines.

A knockdown absent from *every* source line takes a separate route
(:meth:`ContextTransferModel._from_neighbours`): the gene is still measured, so
its behaviour across the rest of the knockdown panel identifies functionally
related genes whose own knockdowns stand in for it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from sklearn.utils.extmath import randomized_svd

from .data import CellLine

EPS = 1e-6


@dataclass
class Hyper:
    """Model switches, all chosen by cross-validation over source lines."""

    temp: float = 0.25        # context-similarity softmax temperature (inf -> uniform)
    smooth: float = 0.30      # weight on perturbation-neighbour smoothing
    n_neighbors: int = 15     # neighbours used by that smoothing
    shrink: float = 1.0       # strength of the reliability shrinkage (0 -> off)
    gamma: float = 0.5        # context-modulation exponent (0 -> off)
    mod_clip: float = 3.0     # clip on the modulation ratio
    beta: float = 1.0         # global calibration scale
    use_global: float = 0.15  # weight on the generic across-perturbation response
    on_target: bool = True    # apply explicit on-target knockdown
    rank: int = 80            # program-basis rank used for denoising
    rank_mix: float = 0.0     # how much of the denoised effect to blend in (0 -> off)
    unseen_k: int = 25        # neighbours used for a knockdown absent from every source
    renorm: float = 0.0       # restore each effect's pre-denoising magnitude (0 -> off)


# --------------------------------------------------------------------------


@dataclass
class SourceBank:
    """Everything derivable from the source lines before hyperparameters apply.

    Cross-validation refits the model dozens of times per fold, but the source
    consensus, its reliability, the perturbation-similarity graph and the
    measured knockdown efficiency all stay put across most of the search.
    Computing them once per fold turns a ~7 s fit into a ~0.4 s one.
    """

    names: np.ndarray
    consensus: np.ndarray        # (n_pert, G) context-weighted mean effect
    global_effect: np.ndarray    # (G,) generic across-perturbation response
    reliability: np.ndarray      # (n_pert,) fraction of the effect that is signal
    sim: np.ndarray              # (n_pert, n_pert) cosine between effects
    kd: np.ndarray               # (n_pert,) residual on-target expression, NaN if unseen
    kd_fallback: float           # global median residual, used where kd is NaN
    basis: np.ndarray            # (max_rank, G) response-program basis
    center: np.ndarray           # (G,) mean effect the basis is taken around
    gene_sig: np.ndarray         # (G, n_pert) each gene's response across knockdowns
    pert_col: np.ndarray         # (n_pert,) gene column of each knocked-down gene, -1 if absent

    MAX_RANK = 400

    @classmethod
    def build(cls, sources: list[CellLine], weights: np.ndarray) -> "SourceBank":
        names = np.array(sorted({n for s in sources for n in s.names}))
        index = {n: i for i, n in enumerate(names)}
        n_pert, n_gene = names.size, sources[0].pert.shape[1]

        acc = np.zeros((n_pert, n_gene))
        wsum = np.zeros(n_pert)
        rowmaps = []
        for s, w in zip(sources, weights):
            rows = np.array([index[n] for n in s.names])
            acc[rows] += w * s.delta
            wsum[rows] += w
            rowmaps.append(rows)

        consensus = acc / np.maximum(wsum, EPS)[:, None]
        global_effect = consensus[wsum > 0].mean(axis=0)
        reliability = cls._reliability(consensus, sources, rowmaps, n_pert)

        z = consensus - consensus.mean(axis=1, keepdims=True)
        z /= np.linalg.norm(z, axis=1, keepdims=True) + EPS
        sim = z @ z.T
        np.fill_diagonal(sim, -np.inf)

        # Perturbation responses are not 6,000-dimensional: they collapse onto a
        # far smaller set of shared programs (translation stress, cell-cycle
        # arrest, the integrated stress response).  Keeping a basis lets the
        # model project a noisy transferred effect onto directions the source
        # data actually supports.  Taken once here, from the raw consensus, so
        # the basis is a property of the source lines rather than of the
        # smoothing hyperparameters.
        center = consensus.mean(axis=0)
        k = min(cls.MAX_RANK, min(consensus.shape) - 1)
        basis = randomized_svd(consensus - center, n_components=k,
                               random_state=0)[2]

        # How each *gene* behaves across the whole knockdown panel.  Genes in one
        # complex move together when their partners are silenced, so this is a
        # functional similarity that never touches a gene's own knockdown -- the
        # only handle available for a perturbation absent from every source line.
        gene_sig = (consensus - center).T.copy()
        gene_sig -= gene_sig.mean(axis=1, keepdims=True)
        gene_sig /= np.linalg.norm(gene_sig, axis=1, keepdims=True) + EPS

        symbols = sources[0].symbols
        col = {sym: i for i, sym in enumerate(symbols)}
        pert_col = np.array([col.get(n, -1) for n in names])

        kd, kd_fallback = cls._knockdown(sources, names)
        return cls(names=names, consensus=consensus, global_effect=global_effect,
                   reliability=reliability, sim=sim, kd=kd,
                   kd_fallback=kd_fallback, basis=basis, center=center,
                   gene_sig=gene_sig, pert_col=pert_col)

    @staticmethod
    def _reliability(consensus: np.ndarray, sources: list[CellLine],
                     rowmaps: list[np.ndarray], n_pert: int) -> np.ndarray:
        """How much of a perturbation's measured effect is real signal.

        With two or more source lines, reproducibility across lines is the
        cleanest available estimate: a knockdown whose effect looks the same in
        two unrelated cell lines is measuring biology, one that does not is
        mostly measuring sampling noise from a few dozen cells.  With one source
        line, fall back on a Wiener factor comparing effect energy to the
        typical energy of the weakest decile, which is dominated by noise.
        """
        energy = (consensus ** 2).mean(axis=1)
        floor = np.quantile(energy[energy > 0], 0.10) if np.any(energy > 0) else 1.0
        wiener = energy / (energy + floor)

        pairs = np.zeros(n_pert)
        counts = np.zeros(n_pert)
        for i in range(len(sources)):
            # local row of each global perturbation index, -1 when absent
            loc_i = np.full(n_pert, -1)
            loc_i[rowmaps[i]] = np.arange(rowmaps[i].size)
            for j in range(i + 1, len(sources)):
                loc_j = np.full(n_pert, -1)
                loc_j[rowmaps[j]] = np.arange(rowmaps[j].size)
                both = np.flatnonzero((loc_i >= 0) & (loc_j >= 0))
                if both.size == 0:
                    continue
                a = sources[i].delta[loc_i[both]]
                b = sources[j].delta[loc_j[both]]
                a = a - a.mean(axis=1, keepdims=True)
                b = b - b.mean(axis=1, keepdims=True)
                num = (a * b).sum(axis=1)
                den = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + EPS
                pairs[both] += np.clip(num / den, 0, 1)
                counts[both] += 1

        repro = np.where(counts > 0, pairs / np.maximum(counts, 1), wiener)
        return np.clip(0.5 * repro + 0.5 * wiener, 0.0, 1.0)

    @staticmethod
    def _knockdown(sources: list[CellLine], names: np.ndarray
                   ) -> tuple[np.ndarray, float]:
        """Residual expression of each targeted gene after CRISPRi.

        Knockdown efficiency is a property of the guide pair, not of the cell
        line: the same gene retains a similar residual fraction across contexts
        (Spearman 0.38-0.54 between these four lines).  So the per-gene value
        measured in the source lines transfers, and is much sharper than one
        global constant -- the residual ranges from under 5% to over 30%.

        Returns per-perturbation residuals (NaN where unmeasured) and the global
        median used as the fallback.
        """
        index = {n: i for i, n in enumerate(names)}
        acc = np.zeros(names.size)
        cnt = np.zeros(names.size)
        for s in sources:
            col = {sym: i for i, sym in enumerate(s.symbols)}
            for row, name in enumerate(s.names):
                g = col.get(name)
                if g is None or s.mu[g] < 0.1:
                    continue
                before = np.expm1(s.mu[g])
                if before <= 0:
                    continue
                acc[index[name]] += min(np.expm1(s.pert[row, g]) / before, 1.0)
                cnt[index[name]] += 1

        per_pert = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
        fallback = float(np.nanmedian(per_pert)) if np.any(cnt > 0) else 0.3
        return per_pert, fallback


@dataclass
class ContextTransferModel:
    hp: Hyper = field(default_factory=Hyper)

    # fitted state
    _genes: np.ndarray | None = None
    _symbols: np.ndarray | None = None
    _pert_names: np.ndarray | None = None
    _consensus: np.ndarray | None = None     # (n_pert, G) source consensus effect
    _global: np.ndarray | None = None        # (G,) generic response
    _reliability: np.ndarray | None = None   # (n_pert,)
    _modulation: np.ndarray | None = None    # (G,)
    _target_mu: np.ndarray | None = None     # (G,)
    _kd: np.ndarray | None = None            # (n_pert,) residual on-target expression
    _kd_fallback: float = 0.15               # used for perturbations never seen
    _gene_sig: np.ndarray | None = None      # (G, n_pert) per-gene response profile
    _pert_col: np.ndarray | None = None      # (n_pert,) gene column of each knockdown

    # ---------------------------------------------------------------- fitting

    def fit(self, sources: list[CellLine], target_mu: np.ndarray,
            symbols: np.ndarray, bank: SourceBank | None = None
            ) -> "ContextTransferModel":
        """Fit using source-line perturbations plus the target's controls only.

        ``bank`` lets a caller reuse the hyperparameter-independent precomputation
        across a search; it must have been built for these same source lines.
        """
        self._symbols = symbols
        self._target_mu = target_mu

        weights = self._context_weights(sources, target_mu)
        if bank is None:
            bank = SourceBank.build(sources, weights)
        self._build_consensus(bank)
        self._modulation = self._context_modulation(sources, target_mu)
        self._kd = bank.kd
        self._kd_fallback = bank.kd_fallback
        self._gene_sig = bank.gene_sig
        self._pert_col = bank.pert_col
        return self

    def context_weights(self, sources: list[CellLine],
                        target_mu: np.ndarray) -> np.ndarray:
        """Public accessor so a caller can build a matching :class:`SourceBank`."""
        return self._context_weights(sources, target_mu)

    def _context_weights(self, sources: list[CellLine],
                         target_mu: np.ndarray) -> np.ndarray:
        """Softmax over how much each source's baseline resembles the target's.

        Correlation is taken on the genes that vary *between* cell lines; on the
        full transcriptome any two human lines correlate at r > 0.9 and the
        weights carry no signal.
        """
        mus = np.vstack([s.mu for s in sources])
        if len(sources) == 1:
            return np.ones(1)

        spread = np.vstack([mus, target_mu[None, :]]).std(axis=0)
        variable = spread > np.quantile(spread, 0.5)

        corr = np.array([
            np.corrcoef(mu[variable], target_mu[variable])[0, 1] for mu in mus])
        if not np.isfinite(self.hp.temp) or self.hp.temp <= 0:
            return np.ones(len(sources)) / len(sources)
        w = np.exp((corr - corr.max()) / self.hp.temp)
        return w / w.sum()

    def _build_consensus(self, bank: SourceBank) -> None:
        """Apply the hyperparameter-dependent smoothing and shrinkage to a bank."""
        self._global = bank.global_effect

        consensus = self._smooth(bank.consensus, bank.reliability, bank.sim)

        if self.hp.shrink > 0:
            factor = bank.reliability ** self.hp.shrink
            consensus = (factor[:, None] * consensus
                         + (1 - factor)[:, None] * self._global[None, :])

        if self.hp.rank_mix > 0 and self.hp.rank > 0:
            # Denoising by projection sharpens agreement with the measured effect
            # but blurs perturbations together, which discrimination punishes.
            # Blending keeps both, at a ratio cross-validation picks.
            v = bank.basis[:min(self.hp.rank, bank.basis.shape[0])]
            centred = consensus - bank.center
            projected = (centred @ v.T) @ v + bank.center
            consensus = ((1 - self.hp.rank_mix) * consensus
                         + self.hp.rank_mix * projected)

        if self.hp.renorm > 0:
            # Smoothing, shrinkage and projection all pull effects toward each
            # other, which shrinks how much predicted magnitudes differ between
            # perturbations.  Discrimination is a retrieval task scored by L1
            # distance, so that spread is most of the signal: if every prediction
            # comes out the same size, the nearest measured effect is whichever
            # one happens to be that size, not the right one.  Restoring each
            # row's original norm keeps the denoised *shape* without paying for
            # it in magnitude.
            before = np.linalg.norm(bank.consensus, axis=1, keepdims=True)
            after = np.linalg.norm(consensus, axis=1, keepdims=True)
            consensus = consensus * (before / (after + EPS)) ** self.hp.renorm

        self._pert_names = bank.names
        self._consensus = consensus
        self._reliability = bank.reliability

    def _smooth(self, consensus: np.ndarray, reliability: np.ndarray,
                sim: np.ndarray) -> np.ndarray:
        """Blend each effect with its nearest neighbours in effect space."""
        if self.hp.smooth <= 0 or self.hp.n_neighbors <= 0:
            return consensus

        n = sim.shape[0]
        k = min(self.hp.n_neighbors, n - 1)
        nbr = np.argpartition(-sim, k, axis=1)[:, :k]
        rows = np.arange(n)[:, None]
        w = np.clip(sim[rows, nbr], 0, None)
        # weight neighbours by their own reliability too
        w = w * reliability[nbr]
        w = w / (w.sum(axis=1, keepdims=True) + EPS)

        # As a sparse mixing matrix: gathering consensus[nbr] densely would
        # allocate n * k * G floats, which is gigabytes at this problem size.
        weights = sparse.csr_matrix(
            (w.ravel(), (np.repeat(np.arange(n), k), nbr.ravel())), shape=(n, n))
        borrowed = weights @ consensus

        lam = self.hp.smooth * (1.0 - reliability)[:, None]
        return (1 - lam) * consensus + lam * borrowed

    def _context_modulation(self, sources: list[CellLine],
                            target_mu: np.ndarray) -> np.ndarray:
        """Per-gene factor for how expressible an effect is in the target line."""
        if self.hp.gamma == 0:
            return np.ones_like(target_mu)
        src_mu = np.vstack([s.mu for s in sources]).mean(axis=0)
        ratio = (target_mu + 0.05) / (src_mu + 0.05)
        ratio = np.clip(ratio, 1.0 / self.hp.mod_clip, self.hp.mod_clip)
        return ratio ** self.hp.gamma

    # ------------------------------------------------------------- prediction

    def _from_neighbours(self, gene: int | None) -> np.ndarray:
        """Predict a knockdown that appears in no source line at all.

        Nothing in the training data names this gene as a perturbation, but the
        gene is still *measured*, so we know how it responds to every other
        knockdown.  Genes acting together move together under perturbation --
        the observation the whole Perturb-seq literature is built on -- so the
        knockdowns whose targets behave most like this one are the best available
        stand-ins for its own.  This uses no measurement of the gene being
        silenced, only of the gene being watched.
        """
        if gene is None or self._gene_sig is None:
            return self._global

        seen = np.flatnonzero(self._pert_col >= 0)
        if seen.size == 0:
            return self._global

        sim = self._gene_sig[self._pert_col[seen]] @ self._gene_sig[gene]
        k = min(self.hp.unseen_k, seen.size)
        top = seen[np.argpartition(-sim, k - 1)[:k]]
        w = np.clip(self._gene_sig[self._pert_col[top]] @ self._gene_sig[gene], 0, None)
        if w.sum() <= EPS:
            return self._global
        return (w / w.sum()) @ self._consensus[top]

    def predict(self, pert_names: np.ndarray) -> np.ndarray:
        """Predicted effect (delta from the target's control mean) per perturbation."""
        assert self._consensus is not None, "call fit() first"
        index = {n: i for i, n in enumerate(self._pert_names)}
        col = {sym: i for i, sym in enumerate(self._symbols)}
        mu = self._target_mu

        out = np.zeros((pert_names.size, mu.size))
        for i, p in enumerate(pert_names):
            j = index.get(p)
            base = (self._consensus[j] if j is not None
                    else self._from_neighbours(col.get(p)))
            d = (1 - self.hp.use_global) * base + self.hp.use_global * self._global
            d = d * self._modulation
            out[i] = self.hp.beta * d

        if self.hp.on_target:
            for i, p in enumerate(pert_names):
                g = col.get(p)
                if g is None:
                    continue
                j = index.get(p)
                residual = self._kd[j] if j is not None else np.nan
                if not np.isfinite(residual):
                    residual = self._kd_fallback
                before = np.expm1(mu[g])
                out[i, g] = np.log1p(before * residual) - mu[g]
        return out

    def predict_expression(self, pert_names: np.ndarray) -> np.ndarray:
        """Absolute post-perturbation profile, which is what the challenge scores."""
        return self._target_mu[None, :] + self.predict(pert_names)


# --------------------------------------------------------------------------
# baselines
# --------------------------------------------------------------------------


@dataclass
class ControlBaseline:
    """Predict the unperturbed control profile: delta is identically zero.

    Deceptively strong on mean absolute error, and at chance on discrimination.
    """

    _mu: np.ndarray | None = None

    def fit(self, sources: list[CellLine], target_mu: np.ndarray,
            symbols: np.ndarray) -> "ControlBaseline":
        self._mu = target_mu
        return self

    def predict(self, pert_names: np.ndarray) -> np.ndarray:
        return np.zeros((pert_names.size, self._mu.size))

    def predict_expression(self, pert_names: np.ndarray) -> np.ndarray:
        return np.repeat(self._mu[None, :], pert_names.size, axis=0)


@dataclass
class GlobalMeanBaseline:
    """Predict the same average perturbation response for every knockdown.

    This is the reference model ``cell-eval`` builds with ``build_base_mean_adata``
    and the one leaderboard scores are normalised against.
    """

    _mu: np.ndarray | None = None
    _delta: np.ndarray | None = None

    def fit(self, sources: list[CellLine], target_mu: np.ndarray,
            symbols: np.ndarray) -> "GlobalMeanBaseline":
        self._mu = target_mu
        self._delta = np.vstack([s.delta.mean(axis=0) for s in sources]).mean(axis=0)
        return self

    def predict(self, pert_names: np.ndarray) -> np.ndarray:
        return np.repeat(self._delta[None, :], pert_names.size, axis=0)

    def predict_expression(self, pert_names: np.ndarray) -> np.ndarray:
        return self._mu[None, :] + self.predict(pert_names)


@dataclass
class NaiveTransferBaseline:
    """Copy the effect measured in the source lines, unweighted and unmodulated.

    The honest thing to beat: it is what anyone would try first, and it already
    uses cross-context information.
    """

    _mu: np.ndarray | None = None
    _names: np.ndarray | None = None
    _delta: np.ndarray | None = None
    _global: np.ndarray | None = None

    def fit(self, sources: list[CellLine], target_mu: np.ndarray,
            symbols: np.ndarray) -> "NaiveTransferBaseline":
        self._mu = target_mu
        names = sorted({n for s in sources for n in s.names})
        index = {n: i for i, n in enumerate(names)}
        acc = np.zeros((len(names), target_mu.size))
        cnt = np.zeros(len(names))
        for s in sources:
            rows = np.array([index[n] for n in s.names])
            acc[rows] += s.delta
            cnt[rows] += 1
        self._delta = acc / np.maximum(cnt, 1)[:, None]
        self._global = self._delta[cnt > 0].mean(axis=0)
        self._names = np.array(names)
        return self

    def predict(self, pert_names: np.ndarray) -> np.ndarray:
        index = {n: i for i, n in enumerate(self._names)}
        out = np.zeros((pert_names.size, self._mu.size))
        for i, p in enumerate(pert_names):
            j = index.get(p)
            out[i] = self._delta[j] if j is not None else self._global
        return out

    def predict_expression(self, pert_names: np.ndarray) -> np.ndarray:
        return self._mu[None, :] + self.predict(pert_names)
