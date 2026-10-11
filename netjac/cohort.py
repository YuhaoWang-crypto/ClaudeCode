"""
cohort — virtual cohorts on an exact kinetic model.

An individual = the curated model with every enzyme amount multiplied by an
independent log-normal factor (inter-individual expression variation, CV
`cv_enzyme`) and optionally a SHARED factor (global "basal metabolic rate",
CV `cv_global`) that moves all enzymes together. A phenotype is a dict of
parameter overrides applied on top (e.g. loss of an allosteric feedback).

Assays are simulated with multiplicative log-normal noise (CV `cv_assay`,
~5-15 % for MS / ELISA / quantitative LFA).

⚠️ The cohort is synthetic: it tests whether a READOUT can in principle
discriminate a mechanistically defined phenotype on real kinetics. It is not
evidence that the phenotype exists in patients.
"""
import numpy as np


def lognormal(rng, cv, size=None):
    s = np.sqrt(np.log(1 + cv ** 2))
    return np.exp(rng.normal(-s ** 2 / 2, s, size))


def make_individuals(model, n, rng, cv_enzyme=0.25, cv_global=0.0,
                     extra_params=None):
    """Return list of dicts {rid: factor} (+ '__global__')."""
    ez = list(model.enzyme_params())
    out = []
    for _ in range(n):
        g = lognormal(rng, cv_global) if cv_global > 0 else 1.0
        ind = {rid: g * lognormal(rng, cv_enzyme) for rid in ez}
        ind["__global__"] = g
        if extra_params:
            for p, cv in extra_params.items():
                ind[f"param:{p}"] = lognormal(rng, cv)
        out.append(ind)
    return out


class applied:
    """Context manager: apply phenotype + individual to the model, restore on
    exit. Usage:  with applied(model, base, pheno, ind): x = model.steady_state()"""

    def __init__(self, model, base_params, phenotype, individual):
        self.m, self.base = model, base_params
        self.pheno, self.ind = phenotype or {}, individual or {}

    def __enter__(self):
        self.m.set_params(self.base)
        self.m.set_params(self.pheno)
        for k, fac in self.ind.items():
            if k == "__global__":
                continue
            if k.startswith("param:"):
                p = k[6:]
                self.m.set_params({p: self.m.get(p) * fac})
            else:
                self.m.scale_enzyme(k, fac)
        return self.m

    def __exit__(self, *a):
        self.m.set_params(self.base)


def assay(values, rng, cv):
    return np.asarray(values) * lognormal(rng, cv, np.shape(values))
