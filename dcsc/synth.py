"""
Synthetic single-cell world with KNOWN ground truth, used to score every DCSC
claim against truth.

Each cell is generated from four separate factor blocks plus a batch effect:

    log mu = base[leaf] + (W_mod * gain[leaf]) @ a      (identity, programs)
                        + W_mem @ m + W_env @ e         (memory, environment)
                        + batch_shift[b]
    counts ~ NegBinomial(libsize * softmax(log mu), theta)

Deliberately hard parts (so the prototype is not a toy that wins by construction):
  * program gene sets overlap each other (random sampling from one pool);
  * Dissociation_stress shares genes with TNF_NFkB / Heat_shock and Culture_2D
    shares genes with EMT / Proliferation, so a naive gene-set score leaks
    artifacts into programs;
  * program response gain differs by cell type (cross-identity transfer is
    imperfect);
  * counts are overdispersed and sparse, and the encoder only ever sees
    log1p-normalised counts, so its linear model is misspecified.

⚠️ This is synthetic: it validates the algorithms and the bookkeeping, not
biology. Real loadings would come from stimulation atlases (cytokine screens,
Tahoe-100M drug contrasts) and dissociation-protocol experiments.
"""
import numpy as np

PROGRAMS = [
    "IFNa", "IFNg", "TNF_NFkB", "IL6_JAK_STAT3", "TGFb", "Hypoxia",
    "Prolif_G2M", "Prolif_S", "ER_stress_UPR", "Oxphos", "Glycolysis", "mTORC1",
    "p53_DNAdamage", "Wnt", "Notch", "Hedgehog", "Androgen", "Estrogen",
    "IL4_STAT6", "IL10_STAT3", "Cholesterol", "Heat_shock", "Apoptosis", "EMT",
    "Antigen_pres", "Complement", "Inflammasome", "AhR", "NRF2_oxidative",
    "Autophagy", "SASP", "Type2_cytokine",
]
MEMORY = ["Trained_immunity", "TOX_exhaustion", "LPS_tolerance", "Hypoxic_memory"]
ENV = ["Dissociation_stress", "Culture_2D", "Mito_damage", "Ambient_RNA"]
TREE = {
    "Myeloid": ["CD14_Mono", "CD16_Mono", "Macrophage"],
    "Lymphoid": ["CD4_T", "CD8_T", "NK"],
    "Stromal": ["Fibroblast", "Endothelial", "Epithelial"],
}
LINEAGES = list(TREE)
LEAVES = [(lin, sub) for lin, subs in TREE.items() for sub in subs]
P, NM, NE = len(PROGRAMS), len(MEMORY), len(ENV)
K_LATENT = P + NM + NE          # 40 latent coordinates after identity
N_ID_GENES = 600


class World:
    def __init__(self, n_genes=2000, seed=0, n_batches=6, theta=8.0):
        rng = np.random.default_rng(seed)
        G = n_genes
        self.G, self.theta, self.n_batches = G, theta, n_batches

        # ---- identity: lineage + subtype markers on genes [0, 600) ----------
        base = rng.normal(0, 1.2, G)
        lin_eff = rng.normal(0, 0.25, (len(LINEAGES), G))
        lin_markers = np.array_split(np.arange(0, 240), len(LINEAGES))
        sub_markers = np.array_split(np.arange(240, N_ID_GENES), len(LEAVES))
        self.leaf_base = np.zeros((len(LEAVES), G))
        for i, (lin, _) in enumerate(LEAVES):
            li = LINEAGES.index(lin)
            b = base + lin_eff[li] + rng.normal(0, 0.1, G)
            b[lin_markers[li]] += 2.5
            b[sub_markers[i]] += 1.8
            self.leaf_base[i] = b
        self.leaf_lineage = np.array([LINEAGES.index(l) for l, _ in LEAVES])

        # ---- programs: 40-gene sets drawn from one shared pool (overlaps) ---
        pool = np.arange(N_ID_GENES, G - 200)
        self.W_mod = np.zeros((G, P))
        self.prog_sets = []
        for p in range(P):
            s = rng.choice(pool, 40, replace=False)
            sign = np.where(rng.random(40) < 0.8, 1.0, -1.0)
            self.W_mod[s, p] = sign * rng.uniform(0.6, 1.4, 40)
            self.prog_sets.append(s)

        # ---- memory + environment: own pool + deliberate overlap ------------
        own = np.arange(G - 200, G)
        self.W_mem = np.zeros((G, NM))
        self.W_env = np.zeros((G, NE))
        own_split = np.array_split(rng.permutation(own), NM + NE)
        for j in range(NM):
            s = own_split[j]
            self.W_mem[s, j] = rng.uniform(0.6, 1.4, len(s))
            self.W_mem[rng.choice(pool, 10, replace=False), j] = 0.6
        P_IDX = {n: i for i, n in enumerate(PROGRAMS)}
        overlaps = {  # env factor -> programs it shares genes with
            "Dissociation_stress": ["TNF_NFkB", "Heat_shock"],
            "Culture_2D": ["EMT", "Prolif_G2M"],
            "Mito_damage": ["Oxphos"],
            "Ambient_RNA": [],
        }
        for j, name in enumerate(ENV):
            s = own_split[NM + j]
            self.W_env[s, j] = rng.uniform(0.6, 1.4, len(s))
            for prog in overlaps[name]:
                shared = self.prog_sets[P_IDX[prog]][:12]
                self.W_env[shared, j] = np.abs(self.W_mod[shared, P_IDX[prog]])
        # ambient RNA = contamination by epithelial markers
        self.W_env[sub_markers[-1], ENV.index("Ambient_RNA")] = 1.0

        # ---- identity-specific response gain + batch shifts ------------------
        self.gain = np.exp(rng.normal(0, 0.2, (len(LEAVES), P)))
        self.batch_shift = rng.normal(0, 0.15, (n_batches, G))

    # ----------------------------------------------------------------------
    def log_mu(self, leaf, a, m, e, batch):
        """leaf/batch: int arrays (n,); a (n,P); m (n,NM); e (n,NE)."""
        lm = self.leaf_base[leaf].copy()
        lm += np.einsum("np,gp->ng", a * self.gain[leaf], self.W_mod)
        lm += m @ self.W_mem.T + e @ self.W_env.T
        lm += self.batch_shift[batch]
        return lm

    def sample_counts(self, leaf, a, m, e, batch, rng):
        lm = self.log_mu(leaf, a, m, e, batch)
        lm -= lm.max(1, keepdims=True)
        p = np.exp(lm)
        p /= p.sum(1, keepdims=True)
        lib = np.exp(rng.normal(np.log(2500), 0.4, len(leaf)))[:, None]
        mu = lib * p
        lam = rng.gamma(self.theta, mu / self.theta)
        return rng.poisson(lam).astype(np.float32)

    def sample_norm(self, leaf, a, m, e, batch, rng, chunk=4000):
        """Generate cells and return log1p(CP10k) matrix (float32)."""
        out = []
        for s in range(0, len(leaf), chunk):
            sl = slice(s, s + chunk)
            c = self.sample_counts(leaf[sl], a[sl], m[sl], e[sl], batch[sl], rng)
            out.append(normalize(c))
        return np.vstack(out)


def normalize(counts):
    return np.log1p(counts / np.maximum(counts.sum(1, keepdims=True), 1) * 1e4).astype(np.float32)


def random_factors(n, rng, p_prog=0.15, p_mem=0.2):
    """Independent, sparse ground-truth factors for test cells."""
    a = rng.normal(0, 0.08, (n, P))
    act = rng.random((n, P)) < p_prog
    a[act] = rng.choice([-1, 1], act.sum()) * rng.uniform(0.6, 1.5, act.sum())
    m = np.where(rng.random((n, NM)) < p_mem, rng.uniform(0.5, 1.5, (n, NM)), 0.0)
    e = np.zeros((n, NE))
    e[:, 0] = rng.uniform(0, 1.5, n)                                # dissociation
    e[:, 1] = (rng.random(n) < 0.5).astype(float)                   # 2D culture
    e[:, 2] = np.where(rng.random(n) < 0.3, rng.uniform(0, 1, n), 0)  # mito
    e[:, 3] = rng.uniform(0, 0.5, n)                                # ambient
    return a, m, e
