"""Per-cell-type identifiability of IFNa vs IFNg reference responses.

Question asked: if a target sample of this cell type carried a full-strength
reference response plus measurement noise, could the score for each stimulus be
told apart from zero?

Model (unscaled log2FC): b = x_a * a + x_g * g + e, e ~ N(0, sigma^2 I).
x is in units of "one full reference response". Its covariance is
sigma^2 (A^T A)^-1, so z_k = 1 / SE(x_k) is the z-score a full-strength response
of stimulus k would reach. Collinearity and weak responses both shrink z.

Caveats: sigma is an assumption (no donor-level data in the aggregate source),
and genes are treated as independent, which overstates information when genes
co-vary. The z threshold is therefore set conservatively; it is an engineering
gate, not a validated detection limit.
"""
import numpy as np

SIGMA_REF = 0.5      # assumed per-gene noise of a target log2FC (log2 units)
Z_MIN = 5.0          # full-strength response must reach this z to count as identifiable
MIN_GENES = 30

# Copied from the author's IFN_Score_Generate.R. Note "Memory_CD8" does not match
# the matrix label "CD8_Memory", so in the author code the rule never fires for it.
AUTHOR_IFNG_NONRESPONDING = {"NK.CD56hi", "NK.CD56dim", "gdT", "MAIT", "Memory_CD8", "Plasma"}


def check(a, g, sigma=SIGMA_REF, z_min=Z_MIN, min_genes=MIN_GENES):
    a = np.asarray(a, float)
    g = np.asarray(g, float)
    n = len(a)
    A = np.column_stack([a, g])
    gram = A.T @ A
    cos = float(a @ g / np.linalg.norm(a) / np.linalg.norm(g))
    cond = float(np.sqrt((1 + abs(cos)) / (1 - abs(cos))))
    cov = sigma ** 2 * np.linalg.inv(gram)
    se = np.sqrt(np.diag(cov))
    z = 1.0 / se
    ok_a = n >= min_genes and z[0] >= z_min
    ok_g = n >= min_genes and z[1] >= z_min
    status = {(True, True): "both", (True, False): "ifna_only",
              (False, True): "ifng_only", (False, False): "none"}[(ok_a, ok_g)]
    return {
        "n_genes": n,
        "rms_ifna": float(np.sqrt(np.mean(a ** 2))),
        "rms_ifng": float(np.sqrt(np.mean(g ** 2))),
        "cosine": cos,
        "condition_number": cond,
        "se_ifna": float(se[0]), "se_ifng": float(se[1]),
        "z_ifna": float(z[0]), "z_ifng": float(z[1]),
        "status": status,
    }


def identifiable_stimuli(status):
    return {"both": ["IFNA", "IFNG"], "ifna_only": ["IFNA"],
            "ifng_only": ["IFNG"], "none": []}[status]
