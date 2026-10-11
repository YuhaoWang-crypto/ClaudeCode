"""
core — KineticModel: any SBML kinetic model -> steady state, state Jacobian,
elasticities, control matrix (dFlux/dEnzyme), with the MCA theorems checked.

Conventions
  x   : independent state = floating species NOT defined by assignment rules
  v   : reaction rates
  N   : stoichiometry restricted to x, already divided by compartment volume,
        so that dx/dt = N v  (verified numerically at load time)
  B   : orthonormal basis of range(N). Conservation laws (NAD+NADH, ...) make
        the full Jacobian singular for STRUCTURAL reasons; the dynamics live on
        x = x* + B z, so the dynamical Jacobian is  Jz = B^T (N eps) B.
        Its eigenvalues are the physically meaningful ones (no fake zeros).
  eps : unscaled elasticities dv/dx  (rows = reactions, cols = state species);
        a nonzero eps[r, k] is an EDGE: reaction r senses species k
        (substrate, product or allosteric effector).
  C^S : unscaled concentration control  dx*/dv  = -B Jz^{-1} B^T N
  C^J : unscaled flux control           dJ*/dv  = I + eps C^S
        scaled C^J_ij = dln J_i / dln E_j  ( = the "dFlux_i/dEnzyme_j" matrix )
"""
import os
import subprocess
import numpy as np
from scipy.linalg import orth

from netjac import CACHE

MIRROR = "https://raw.githubusercontent.com/biomodels/{id}/master/{id}/{id}.xml"


def fetch_biomodel(biomd_id):
    """Download a curated BioModels SBML via the GitHub mirror (cached)."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, f"{biomd_id}.xml")
    if not os.path.exists(path) or os.path.getsize(path) < 500:
        subprocess.run(["curl", "-sf", "--max-time", "90",
                        MIRROR.format(id=biomd_id), "-o", path], check=True)
    return path


class quiet:
    """Silence C-level stdout/stderr (CVODE warnings on stiff/failing runs)."""

    def __enter__(self):
        import sys
        sys.stdout.flush(); sys.stderr.flush()
        self.null = os.open(os.devnull, os.O_WRONLY)
        self.saved = [os.dup(1), os.dup(2)]
        os.dup2(self.null, 1); os.dup2(self.null, 2)

    def __exit__(self, *a):
        os.dup2(self.saved[0], 1); os.dup2(self.saved[1], 2)
        for fd in self.saved + [self.null]:
            os.close(fd)


def _load_sbml_string(path):
    """Promote reaction-local parameters to globals (e.g. vPFK.VmPFK ->
    vPFK_VmPFK) so every kinetic constant is addressable by name."""
    import libsbml
    doc = libsbml.readSBML(path)
    props = libsbml.ConversionProperties()
    props.addOption("promoteLocalParameters", True)
    doc.convert(props)
    return libsbml.writeSBMLToString(doc), doc


class KineticModel:
    def __init__(self, source, name=None):
        import roadrunner
        roadrunner.Logger.setLevel(roadrunner.Logger.LOG_FATAL)
        roadrunner.Logger.disableLogging()
        path = fetch_biomodel(source) if source.startswith("BIOMD") else source
        sbml, doc = _load_sbml_string(path)
        self.name = name or source
        self.rr = roadrunner.RoadRunner(sbml)
        m = doc.getModel()
        rule_vars = {r.getVariable() for r in m.getListOfRules()
                     if r.isAssignment()}
        self.all_species = list(self.rr.model.getFloatingSpeciesIds())
        self.state_ids = [s for s in self.all_species if s not in rule_vars]
        self.rxn_ids = list(self.rr.model.getReactionIds())
        self.param_ids = list(self.rr.model.getGlobalParameterIds())
        self._sidx = [self.all_species.index(s) for s in self.state_ids]
        self._default_params = {p: self.rr[p] for p in self.param_ids}
        self._default_x = self.x()
        Nfull = self.rr.getFullStoichiometryMatrix()
        rows = list(Nfull.rownames)
        N = np.array(Nfull)[[rows.index(s) for s in self.state_ids], :]
        vol = np.array([m.getCompartment(m.getSpecies(s).getCompartment())
                        .getSize() for s in self.state_ids])
        self.N = N / vol[:, None]
        self.B = orth(self.N)
        self.n_cons = len(self.state_ids) - self.B.shape[1]
        self._check_rhs()
        self._enzyme_params = None

    # ---------------------------------------------------------------- state
    def x(self):
        return np.array([self.rr[s] for s in self.state_ids])

    def set_x(self, x):
        for s, val in zip(self.state_ids, x):
            self.rr[s] = float(val)

    def f(self, x=None):
        if x is not None:
            self.set_x(x)
        return np.asarray(self.rr.model.getFloatingSpeciesConcentrationRates()
                          )[self._sidx]

    def v(self, x=None):
        if x is not None:
            self.set_x(x)
        return np.asarray(self.rr.model.getReactionRates()).copy()

    def _check_rhs(self):
        x = self.x()
        f, v = self.f(x), self.v(x)
        err = np.abs(f - self.N @ v).max() / (np.abs(f).max() + 1e-12)
        if err > 1e-6:
            raise ValueError(f"dx/dt != N v (rel err {err:.2e}); model has "
                             "rate rules or non-unit compartments not handled")

    # ----------------------------------------------------------- parameters
    def get(self, p):
        return self.rr[p]

    def set_params(self, d):
        for k, val in d.items():
            self.rr[k] = float(val)

    def reset(self):
        """Restore SBML default parameters and initial state."""
        self.set_params(self._default_params)
        self.set_x(self._default_x)

    def snapshot(self):
        return {p: self.rr[p] for p in self.param_ids}

    # --------------------------------------------------------- steady state
    def steady_state(self, t_end=2000.0, x0=None, tol=1e-7, max_rounds=6):
        """Integrate to steady state (robust to conservation-law singularity
        that defeats Newton solvers). Returns x*; raises if it never settles
        (e.g. the system is on a limit cycle -> use relax/spectrum tools)."""
        if x0 is not None:
            self.set_x(x0)
        for _ in range(max_rounds):
            try:
                with quiet():
                    self.rr.simulate(0, t_end, 3)
            except RuntimeError:
                raise RuntimeError(f"{self.name}: integration failed "
                                   "(runaway accumulation / no steady state)")
            x = self.x()
            res = np.abs(self.f(x)).max() / (np.abs(self.N).sum(1).max()
                                             * np.abs(self.v(x)).max() + 1e-12)
            if res < tol:
                return x
            t_end *= 4
        raise RuntimeError(f"{self.name}: no steady state (residual {res:.1e})"
                           " - system may be oscillating")

    def fixed_point(self, x0=None, tol=1e-9):
        """Newton solve of the fixed point WITHOUT integrating, so it also
        finds UNSTABLE steady states (inside an oscillatory regime). The
        unknown is z with x = x0 + B z, which keeps every conservation law at
        the value it has in x0."""
        from scipy.optimize import root
        x0 = self.x() if x0 is None else np.asarray(x0, float)
        B = self.B

        def g(z):
            return B.T @ self.f(np.maximum(x0 + B @ z, 1e-14))
        sol = root(g, np.zeros(B.shape[1]), method="hybr", tol=1e-12)
        x = np.maximum(x0 + B @ sol.x, 1e-14)
        res = np.abs(self.f(x)).max() / (np.abs(self.v(x)).max() + 1e-12)
        if res > tol * 1e3:
            raise RuntimeError(f"{self.name}: fixed point not found ({res:.1e})")
        return x

    # ----------------------------------------------------------- jacobians
    def _fd(self, fun, x, h=1e-6):
        x = np.asarray(x, float)
        f0 = fun(x)
        out = np.zeros((len(f0), len(x)))
        for k in range(len(x)):
            dx = h * max(abs(x[k]), 1e-6)
            xp, xm = x.copy(), x.copy()
            xp[k] += dx
            xm[k] -= dx
            out[:, k] = (fun(xp) - fun(xm)) / (2 * dx)
        self.set_x(x)
        return out

    def elasticities(self, x):
        """Unscaled eps = dv/dx (reactions x state species)."""
        return self._fd(self.v, x)

    def jacobian(self, x):
        """Full state Jacobian d(dx/dt)/dx (has n_cons structural zeros)."""
        return self._fd(self.f, x)

    def reduced_jacobian(self, x, eps=None):
        eps = self.elasticities(x) if eps is None else eps
        return self.B.T @ self.N @ eps @ self.B

    # ------------------------------------------------------------------ MCA
    def control(self, x):
        """Unscaled + scaled control matrices and theorem residuals."""
        v = self.v(x)
        eps = self.elasticities(x)
        Jz = self.B.T @ self.N @ eps @ self.B
        CS = -self.B @ np.linalg.solve(Jz, self.B.T @ self.N)
        CJ = np.eye(len(v)) + eps @ CS
        with np.errstate(divide="ignore", invalid="ignore"):
            CJs = CJ * v[None, :] / v[:, None]
            CSs = CS * v[None, :] / x[:, None]
        CJs[~np.isfinite(CJs)] = 0.0
        CSs[~np.isfinite(CSs)] = 0.0
        # theorem checks (unscaled forms, valid with conservation laws)
        K = _null(self.N)
        summ = np.abs(CJ @ K - K).max() if K.size else 0.0
        conn = np.abs(CJ @ eps @ self.B).max()
        return {"x": x, "v": v, "eps": eps, "Jz": Jz, "CS": CS, "CJ": CJ,
                "CJ_scaled": CJs, "CS_scaled": CSs,
                "summation_residual": summ, "connectivity_residual": conn}

    # -------------------------------------------------- enzyme (Vmax) knobs
    def enzyme_params(self):
        """Auto-detect, per reaction, the parameter set that scales its rate
        linearly (the 'enzyme amount'): doubling it doubles v_r and leaves the
        other rates unchanged. Detected numerically -> works for any SBML."""
        if self._enzyme_params is not None:
            return self._enzyme_params
        x_now = self.x()
        # generic test state: every species strictly positive (at SBML initial
        # conditions products are often 0, which hides reverse Vmax terms)
        pos = x_now[x_now > 0]
        fill = 0.1 * np.median(pos) if len(pos) else 1e-3
        x = np.where(x_now > 0, x_now, fill) * np.exp(
            0.05 * np.sin(np.arange(len(x_now)) + 1.0))
        v0 = self.v(x)
        affects = {}
        for p in self.param_ids:
            p0 = self.rr[p]
            if p0 == 0:
                continue
            self.rr[p] = 2 * p0
            dv = self.v(x) - v0
            self.rr[p] = p0
            hit = np.where(np.abs(dv) > 1e-9 * (np.abs(v0) + 1e-12))[0]
            if len(hit) == 1:
                affects.setdefault(int(hit[0]), []).append(p)
        out = {}
        for r, cands in affects.items():
            rid = self.rxn_ids[r]
            if abs(v0[r]) < 1e-12:
                continue
            def linear(ps):
                old = {p: self.rr[p] for p in ps}
                self.set_params({p: 2 * q for p, q in old.items()})
                ratio = self.v(x)[r] / v0[r]
                self.set_params(old)
                return abs(ratio - 2) < 1e-6
            # prefer the full Vmax/kcat group (e.g. GAPDH forward+reverse):
            # a single one can LOOK linear at states where the other
            # direction carries no flux (e.g. SBML initial conditions)
            vm = [p for p in cands if "vm" in p.lower() or "kcat" in p.lower()]
            if vm and linear(vm):
                out[rid] = vm
                continue
            singles = [p for p in cands if linear([p])]
            if singles:
                out[rid] = [singles[0]]
        self.set_x(x_now)
        self._enzyme_params = out
        return out

    def scale_enzyme(self, rid, factor):
        for p in self.enzyme_params()[rid]:
            self.rr[p] = self.rr[p] * factor

    def control_fd(self, x, rel=1e-3, t_end=500.0):
        """Scaled flux-control matrix by brute force: perturb each enzyme by
        +-rel and re-solve the steady state. Independent of the analytic
        formula -> used to validate it."""
        J0 = self.v(x)
        ez = self.enzyme_params()
        cols = {}
        for rid in ez:
            out = []
            for s in (+1, -1):
                self.scale_enzyme(rid, 1 + s * rel)
                xs = self.steady_state(t_end=t_end, x0=x)
                out.append(self.v(xs))
                self.scale_enzyme(rid, 1 / (1 + s * rel))
            with np.errstate(divide="ignore", invalid="ignore"):
                cols[rid] = (np.log(np.abs(out[0])) - np.log(np.abs(out[1]))
                             ) / (np.log(1 + rel) - np.log(1 - rel))
        self.set_x(x)
        return cols, J0


def _null(A, tol=1e-10):
    u, s, vh = np.linalg.svd(A)
    rank = (s > tol * s.max()).sum()
    return vh[rank:].T
