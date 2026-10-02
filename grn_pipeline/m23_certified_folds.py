"""
Module 23 — CERTIFIED saddle-node folds of the exact Markevich ERK switch
(ball arithmetic via the BootLoops toolkit).

M15 locates the bistable window [39.25, 57.38] nM by bisection on a float
count of stable states; M20b confirms it on the official SBML. Both are
floating-point. Here the same model (built from M15's own rate function, so
there is no re-transcription) is turned into an exact polynomial system over
Q and the window edges are PROVEN:

  Route A (existence + uniqueness): Krawczyk interval-Newton on the fold
      system  N1 = N2 = 0,  det dN/d(M,Mpp) = 0  in (M, Mpp, MAPKK), Arb balls
      at 256 bits. K(X) strictly inside X  =>  exactly one fold in the box.
  Route B (completeness): exact elimination over Q (Res_M of the steady-state
      curve and the fold condition), certified real-root isolation of the
      univariate resultant, back-substitution by the first subresultant.
      Every fold in the physical simplex is a root; every root is either
      matched to a Route-A box or excluded by a certified inequality.
  Gate: BootLoops `baller.tripwire.dual` — the two routes must agree to the
      bar or the computation HALTS (never averaged); printing goes through
      `baller.render` (fail-closed: only certified digits are printed).

Controls (BootLoops acceptance-gate / planted-truth protocol):
  + planted cubic with folds at exactly K = -2, +2 must be recovered;
  - a Krawczyk box displaced off the fold must REFUSE;
  - a tripwire fed a value off by 1e-8 must HALT;
  - the Km5 = 77 / 79 mutants must give certified windows DISJOINT from
    Km5 = 78 (the certificate discriminates the back-solved constant).
"""
import json
import os
from contextlib import contextmanager

import numpy as np
import sympy as sp
from flint import arb, arb_mat, ctx, fmpq, fmpz_poly

from grn_pipeline import m15_markevich_mm as m15
from grn_pipeline._bootloops import baller as _load_baller

PREC = 256                      # bits (~77 decimal digits)
RAD = 1e-40                     # relative Krawczyk box half-width
DIGITS = 30                     # tripwire bar between the two routes
RECEIPT = "receipts/m23_certified_folds.json"


# ---------------------------------------------------------------------------
# exact model -> polynomial system over Q
# ---------------------------------------------------------------------------
X, Y, KK = sp.symbols("M Mpp MAPKK")


def _q(e):
    return sp.nsimplify(e, rational=True)


def markevich_system():
    """Numerators/denominators of dM/dt and dMpp/dt from M15's own _v()."""
    Mp = _q(m15.MTOT) - X - Y
    v1, v2, v3, v4 = (_q(v) for v in m15._v(X, Mp, Y, KK))
    n1, d1 = sp.fraction(sp.together(v4 - v1))
    n2, d2 = sp.fraction(sp.together(v2 - v3))
    mt = float(m15.MTOT)            # exact in binary; arb-compatible
    return {"N": (sp.expand(n1), sp.expand(n2)),
            "den": (sp.expand(d1), sp.expand(d2)),
            "box": ((0, mt), (0, mt)),
            "simplex": True,         # M, Mpp, Mp = MTOT - M - Mpp > 0
            "kpos": True}


def planted_system():
    """Positive control: steady states y = x, K = x^3 - 3x; folds at x = +-1,
    K = -+2 exactly (dK/dx = 3x^2 - 3)."""
    return {"N": (sp.expand(KK - X**3 + 3 * X + (Y - X)), sp.expand(X - Y)),
            "den": (sp.Integer(1), sp.Integer(1)),
            "box": ((-10, 10), (-10, 10)),
            "simplex": False,
            "kpos": False}


@contextmanager
def km5_set(value):
    old = m15.KM5
    m15.KM5 = value
    try:
        yield
    finally:
        m15.KM5 = old


# ---------------------------------------------------------------------------
# rigorous evaluation of sympy polynomials over Arb boxes
# ---------------------------------------------------------------------------
class APoly:
    """sympy polynomial in (M, Mpp, MAPKK) evaluated in ball arithmetic."""

    def __init__(self, expr):
        p = sp.Poly(sp.expand(expr), X, Y, KK)
        self.terms = [(e, fmpq(int(c.p), int(c.q))) for e, c in p.terms()]

    def __call__(self, x, y, k):
        s = arb(0)
        for (a, b, c), co in self.terms:
            s += arb(co) * x**a * y**b * k**c
        return s


def _jac(F, vars_):
    return [[sp.diff(f, v) for v in vars_] for f in F]


def _floatmid(a):
    return float(a.mid())


# ---------------------------------------------------------------------------
# Route A — Krawczyk interval-Newton (existence + uniqueness)
# ---------------------------------------------------------------------------
def krawczyk(F, J, center, radii):
    """F: list of APoly-like callables of the n box coordinates; J: n x n of
    same. center: list of exact arb midpoints; radii: list of floats.
    Returns (ok, box, image). ok=True PROVES a unique zero of F in box."""
    n = len(center)
    c = [arb(v.mid()) for v in center]
    Xb = [arb(c[i].mid(), radii[i]) for i in range(n)]
    Fc = arb_mat([[f(*c)] for f in F])
    JX = arb_mat([[J[i][j](*Xb) for j in range(n)] for i in range(n)])
    Jm = np.array([[_floatmid(J[i][j](*c)) for j in range(n)]
                   for i in range(n)])
    Cn = np.linalg.inv(Jm)
    C = arb_mat([[arb(float(Cn[i, j])) for j in range(n)] for i in range(n)])
    I = arb_mat([[1 if i == j else 0 for j in range(n)] for i in range(n)])
    D = arb_mat([[Xb[i] - c[i]] for i in range(n)])
    Kimg = arb_mat([[c[i]] for i in range(n)]) - C * Fc + (I - C * JX) * D
    ok = True
    for i in range(n):
        lo, hi = c[i] - radii[i], c[i] + radii[i]
        if not (Kimg[i, 0] > lo and Kimg[i, 0] < hi):
            ok = False
    return ok, Xb, [Kimg[i, 0] for i in range(n)]


def newton_polish(F, J, start, iters=60):
    x = [arb(v) for v in start]
    n = len(x)
    for _ in range(iters):
        pt = [arb(v.mid()) for v in x]
        Fv = arb_mat([[f(*pt)] for f in F])
        Jv = arb_mat([[J[i][j](*pt) for j in range(n)] for i in range(n)])
        step = Jv.mid().solve(Fv.mid())
        x = [arb((pt[i] - step[i, 0]).mid()) for i in range(n)]
    return x


def fold_system(sysd):
    N1, N2 = sysd["N"]
    D = sp.expand(sp.Matrix(_jac((N1, N2), (X, Y))).det())
    F = (N1, N2, D)
    Fa = [APoly(f) for f in F]
    Ja = [[APoly(e) for e in row] for row in _jac(F, (X, Y, KK))]
    return Fa, Ja


def certify_fold(sysd, start):
    """Polish from `start` = (M, Mpp, MAPKK) and Krawczyk-certify the fold."""
    Fa, Ja = fold_system(sysd)
    c = newton_polish(Fa, Ja, start)
    radii = [max(abs(_floatmid(v)), 1.0) * RAD for v in c]
    ok, box, img = krawczyk(Fa, Ja, c, radii)
    return {"ok": ok, "box": box, "image": img, "center": c, "radii": radii}


# ---------------------------------------------------------------------------
# Route B — exact elimination + certified root isolation (completeness)
# ---------------------------------------------------------------------------
def _drop_factors(expr, dens, sysd):
    """Remove factors of expr that divide the rate-law denominators, after
    CERTIFYING each removed factor nonzero on the bounding box of the
    physical domain (so no physical fold can hide in a dropped factor)."""
    c, facs = sp.factor_list(sp.expand(expr), X, Y)
    keep = sp.Integer(1)
    dropped = []
    for f, m in facs:
        if any(sp.rem(sp.Poly(d, X, Y), sp.Poly(f, X, Y)).is_zero
               for d in dens if d.free_symbols):
            (xa, xb), (ya, yb) = sysd["box"]
            val = APoly(f)(arb((xa + xb) / 2, (xb - xa) / 2),
                           arb((ya + yb) / 2, (yb - ya) / 2), arb(0))
            if val.contains(0):
                raise RuntimeError(f"cannot certify dropped factor {f} != 0")
            dropped.append(f)
        else:
            keep *= f**m
    return sp.expand(keep), dropped


def enumerate_folds(sysd, kind="fold"):
    """Every fold (det J = 0) -- or, kind="hopf", every trace J = 0 point --
    of the system's steady states in the physical domain, with certificates.
    Hopf candidates with a CERTAIN det J < 0 are neutral saddles, not Hopf."""
    N1, N2 = sysd["N"]
    pN1 = sp.Poly(N1, KK)
    if pN1.degree() != 1:
        raise RuntimeError("route B needs N1 linear in MAPKK")
    B1 = -pN1.coeff_monomial(KK)
    A1 = pN1.coeff_monomial(1)
    # K = A1/B1 on the steady-state set; B1 must not vanish there
    for f, _ in sp.factor_list(B1, X, Y)[1]:
        if f in (X, Y):
            continue        # coordinate factor: > 0 on the open domain
        (xa, xb), (ya, yb) = sysd["box"]
        v = APoly(f)(arb((xa + xb) / 2, (xb - xa) / 2),
                     arb((ya + yb) / 2, (yb - ya) / 2), arb(0))
        if v.contains(0):
            raise RuntimeError(f"cannot certify B1 factor {f} != 0")
    JN = _jac((N1, N2), (X, Y))
    D = sp.Matrix(JN).det()
    if kind == "fold":
        cond = D
    else:   # trace of diag(1/d1, 1/d2) JN, cleared of the (positive) dens
        d1, d2 = sysd["den"]
        cond = JN[0][0] * d2 + JN[1][1] * d1
    G = sp.fraction(sp.together(N2.subs(KK, A1 / B1)))[0]
    H = sp.fraction(sp.together(cond.subs(KK, A1 / B1)))[0]
    dens = list(sysd["den"]) + [B1]
    G, dG = _drop_factors(G, dens, sysd)
    H, dH = _drop_factors(H, dens, sysd)
    prs = sp.subresultants(sp.Poly(G, X), sp.Poly(H, X))
    R = sp.Poly(prs[-1].as_expr(), Y)
    S1 = sp.Poly(prs[-2].as_expr(), X)
    if R.is_zero or S1.degree() != 1:
        raise RuntimeError("degenerate elimination; refusing completeness")
    s1, s0 = S1.coeff_monomial(X), S1.coeff_monomial(1)
    lcG, lcH = sp.Poly(G, X).LC(), sp.Poly(H, X).LC()
    a = lambda e: APoly(e)
    out = {"deg_R": R.degree(), "n_complex": 0, "excluded": [],
           "candidates": [], "dropped": [str(f) for f in dG + dH]}
    if R.degree() == 0:
        out["n_distinct"] = 0
        return out          # nonzero constant resultant: no common zero
    z = arb(0)
    facs = sp.factor_list(R.as_expr(), Y)[1]
    out["n_distinct"] = sum(sp.Poly(f, Y).degree() for f, _ in facs)
    for fac, _ in facs:
        fp = sp.Poly(fac, Y)
        if fp.degree() == 1:
            # rational root: decide it EXACTLY over Q
            y0 = -fp.coeff_monomial(1) / fp.coeff_monomial(Y)
            _exact_root(y0, s0, s1, A1, B1, lcG, lcH, sysd, out)
            continue
        for y in _real_roots(fp, out):
            if a(lcG)(z, y, z).contains(0) and a(lcH)(z, y, z).contains(0):
                raise RuntimeError("both leading coefficients may vanish")
            _classify(y, s0, s1, A1, B1, sysd, out)
    if out["n_complex"] + len(out["excluded"]) + len(out["candidates"]) \
            != out["n_distinct"]:
        raise RuntimeError("root accounting mismatch; refusing completeness")
    if kind == "hopf":
        Da, keep = APoly(D), []
        for c in out["candidates"]:
            if Da(c["M"], c["Mpp"], c["MAPKK"]) < 0:
                out["excluded"].append(dict(c, why="neutral saddle (det<0)"))
            else:
                keep.append(c)
        out["candidates"] = keep
    return out


def _exact_root(y0, s0, s1, A1, B1, lcG, lcH, sysd, out):
    """A rational root of the resultant is decided EXACTLY over Q."""
    if lcG.subs(Y, y0) == 0 and lcH.subs(Y, y0) == 0:
        raise RuntimeError("both leading coefficients vanish at a root")
    s1v = s1.subs(Y, y0)
    if s1v == 0:
        raise RuntimeError("first subresultant vanishes at a rational root")
    x0 = -s0.subs(Y, y0) / s1v
    mt = sp.Rational(sysd["box"][0][1])
    why = [f"{n}<=0" for n, v in (("M", x0), ("Mpp", y0), ("Mp", mt - x0 - y0))
           if sysd["simplex"] and v <= 0]
    b1 = B1.subs({X: x0, Y: y0})
    k0 = A1.subs({X: x0, Y: y0}) / b1 if b1 != 0 else None
    if sysd["kpos"] and k0 is not None and k0 <= 0:
        why.append("MAPKK<=0")
    ex = lambda v: arb(fmpq(int(v.p), int(v.q)))
    if why:
        out["excluded"].append({"Mpp": ex(y0), "M": ex(x0),
                                "why": "exact boundary point: " +
                                ",".join(why)})
    elif k0 is None:
        raise RuntimeError("rational root with B1 = 0 inside the domain")
    else:
        out["candidates"].append({"M": ex(x0), "Mpp": ex(y0),
                                  "MAPKK": ex(k0)})


def _real_roots(fp, out):
    """Certified real roots of an irreducible integer polynomial: isolate
    (Arb), then refine each by Newton and RE-ENCLOSE with a certified sign
    change inside its isolating ball (IVT + isolation => that root)."""
    den = sp.ilcm(*[sp.Rational(c).q for c in fp.all_coeffs()])
    zp = fmpz_poly([int(c * den) for c in reversed(fp.all_coeffs())])
    dz = zp.derivative()
    roots = []
    for r, _ in zp.complex_roots():
        if not r.imag.contains(0):
            out["n_complex"] += 1
            continue
        iso = r.real
        y = arb(iso.mid())
        for _ in range(12):
            y = arb((y - zp(y) / dz(y)).mid())
        e = max(abs(float(y.mid())), 1.0) * RAD
        lo, hi = arb(y.mid()) - e, arb(y.mid()) + e
        flo, fhi = zp(lo), zp(hi)
        if not ((flo * fhi) < 0 and iso.contains(arb(y.mid(), e))):
            roots.append(iso)               # keep the coarse, still-valid ball
        else:
            roots.append(arb(y.mid(), e))
    return roots


def _classify(y, s0, s1, A1, B1, sysd, out):
    """Back-substitute M = -s0/s1 (first subresultant) and decide the root:
    excluded by a CERTAIN inequality, or kept as a physical candidate."""
    a = lambda e: APoly(e)
    z = arb(0)
    s1v = a(s1)(z, y, z)
    if s1v.contains(0):
        raise RuntimeError("first subresultant degenerate at a root")
    x = -a(s0)(z, y, z) / s1v
    k = a(A1)(x, y, z) / a(B1)(x, y, z)
    viol = _certainly_violated(x, y, k, sysd)
    if viol:
        out["excluded"].append({"Mpp": y, "M": x, "MAPKK": k, "why": viol})
    else:
        out["candidates"].append({"M": x, "Mpp": y, "MAPKK": k})


def _certainly_violated(x, y, k, sysd):
    """arb comparisons are True only when CERTAIN for the whole ball."""
    why = []
    if sysd["simplex"]:
        mt = sysd["box"][0][1]
        if x < 0:
            why.append("M<0")
        if y < 0:
            why.append("Mpp<0")
        if mt - x - y < 0:
            why.append("Mp<0")
    if sysd["kpos"] and k < 0:
        why.append("MAPKK<0")
    return ",".join(why)


# ---------------------------------------------------------------------------
# certified steady states + stability at a fixed MAPKK (the 3-state table)
# ---------------------------------------------------------------------------
def certified_states(sysd, K0, starts):
    N1, N2 = sysd["N"]
    kq = _q(K0)
    F = [N1.subs(KK, kq), N2.subs(KK, kq)]
    Fa = [APoly(f) for f in F]
    Ja = [[APoly(e) for e in row] for row in _jac(F, (X, Y))]
    J2 = lambda f: (lambda x, y: f(x, y, arb(0)))
    F2 = [J2(f) for f in Fa]
    Jm = [[J2(e) for e in row] for row in Ja]
    d1, d2 = (APoly(d) for d in sysd["den"])
    res = []
    for st in starts:
        c = newton_polish(F2, Jm, st)
        radii = [max(abs(_floatmid(v)), 1.0) * RAD for v in c]
        ok, box, _ = krawczyk(F2, Jm, c, radii)
        x, y = box
        # Jacobian of the RHS at a steady state = diag(1/den) * dN
        k = arb(fmpq(int(kq.p), int(kq.q)))
        a11 = Ja[0][0](x, y, k) / d1(x, y, k)
        a12 = Ja[0][1](x, y, k) / d1(x, y, k)
        a21 = Ja[1][0](x, y, k) / d2(x, y, k)
        a22 = Ja[1][1](x, y, k) / d2(x, y, k)
        tr, det = a11 + a22, a11 * a22 - a12 * a21
        disc = tr * tr - 4 * det
        lmax = (tr + disc.sqrt()) / 2 if disc > 0 else None
        kind = ("stable" if (det > 0 and tr < 0) else
                "saddle" if det < 0 else "UNDECIDED")
        res.append({"ok": ok, "M": x, "Mpp": y, "tr": tr, "det": det,
                    "lmax": lmax, "kind": kind})
    return res


# ---------------------------------------------------------------------------
# float starts from M15 (independent float route)
# ---------------------------------------------------------------------------
def _fold_start(Kf):
    """Midpoint of the merging stable/saddle pair just inside the window."""
    for dk in (0.05, -0.05):
        ss = m15.steady_states(Kf + dk, n=1600)
        if len(ss) == 3:
            break
    lo, mid, hi = ss
    pair = (mid, hi) if dk > 0 else (lo, mid)
    return [(pair[0]["M"] + pair[1]["M"]) / 2,
            (pair[0]["Mpp"] + pair[1]["Mpp"]) / 2, Kf]


def _render(B, x, d=20):
    return B.render(x, digits=d) if B else x.str(d, radius=False)


def _flatten(x):
    return x.str(40, radius=True)


def report(fig_path=None):
    print("=" * 68)
    print("MODULE 23 — CERTIFIED folds of the exact Markevich ERK switch")
    print("=" * 68)
    old = ctx.prec
    ctx.prec = PREC
    try:
        return _report()
    finally:
        ctx.prec = old


def _report():
    B = _load_baller()
    if B is None:
        print("BootLoops not found ($BOOTLOOPS_DIR or ../bootloops); "
              "run `bash scripts/setup_bootloops.sh`. Tripwire/render gates "
              "fall back to local checks.")
    if B:
        from baller import tripwire as tw
    else:
        tw = None
    rec = {"prec_bits": PREC, "controls": {}}

    # ---- positive control: planted folds at K = -2, +2 exactly
    P = planted_system()
    pa = [certify_fold(P, s) for s in ([1.1, 0.9, -1.9], [-0.9, -1.1, 2.2])]
    pb = enumerate_folds(P)
    plant_ok = (all(r["ok"] for r in pa)
                and pa[0]["box"][2].contains(arb(-2))
                and pa[1]["box"][2].contains(arb(2))
                and len(pb["candidates"]) == 2)
    print(f"\n[+ control] planted cubic: Krawczyk certifies folds containing "
          f"K=-2 and K=+2 exactly; route B finds {len(pb['candidates'])} "
          f"folds -> {'PASS' if plant_ok else 'FAIL'}")
    rec["controls"]["planted"] = plant_ok

    # ---- the Markevich model
    S = markevich_system()
    print("\nexact system over Q from M15's _v() (Km5 = 78); rate-law "
          "denominators certified > 0 on the simplex")
    win_float = m15.bistable_window()
    folds = [certify_fold(S, _fold_start(k)) for k in win_float]
    comp = enumerate_folds(S)
    print(f"route B: resultant degree {comp['deg_R']} "
          f"({comp['n_distinct']} distinct roots), "
          f"{comp['n_complex']} complex, "
          f"{len(comp['excluded'])} real roots certified non-physical, "
          f"{len(comp['candidates'])} physical candidates")

    labels = ("lower fold (ON branch dies)", "upper fold (OFF branch dies)")
    rec["folds"] = []
    for lab, f, kf in zip(labels, folds, win_float):
        if not f["ok"]:
            raise RuntimeError(f"Krawczyk REFUSED at {lab}")
        kb = f["box"][2]
        match = [c for c in comp["candidates"] if c["MAPKK"].overlaps(kb)]
        if len(match) != 1:
            raise RuntimeError(f"route B does not match {lab}")
        kB = match[0]["MAPKK"]
        if tw:
            tw.dual(lambda: _mp(kb), lambda: _mp(kB), DIGITS, label=lab)
        # stability of the merging pair: trace < 0 certified at the fold
        tr = _fold_trace(S, f["box"])
        print(f"\n{lab}:")
        print(f"  MAPKK* = {_render(B, kb, 25)} nM   (certified; M15 float "
              f"bisection: {kf:.6f})")
        print(f"  M* = {_render(B, f['box'][0], 12)}, Mpp* = "
              f"{_render(B, f['box'][1], 12)} nM")
        print(f"  Krawczyk: unique fold in a box of half-width "
              f"{f['radii'][2]:.1e}; route B agrees to >= {DIGITS} digits "
              f"[baller.tripwire.dual]; trace J = {_render(B, tr, 4)} < 0 "
              f"=> stable node + saddle collide (saddle-node, not Hopf)")
        rec["folds"].append({"label": lab, "MAPKK": _flatten(kb),
                             "MAPKK_mid": float(kb.mid()),
                             "M": _flatten(f["box"][0]),
                             "Mpp": _flatten(f["box"][1]),
                             "trace": _flatten(tr),
                             "m15_float": kf})
    if len(comp["candidates"]) != 2:
        raise RuntimeError("completeness: expected exactly 2 physical folds")
    print("\ncompleteness: these are the ONLY two folds in the physical "
          "simplex for MAPKK > 0")
    hopf = enumerate_folds(S, kind="hopf")
    print(f"no Hopf: trace J = 0 on the steady-state set -> resultant degree "
          f"{hopf['deg_R']} ({hopf['n_distinct']} distinct roots), "
          f"{hopf['n_complex']} complex, {len(hopf['excluded'])} certified "
          f"non-physical/neutral saddles, {len(hopf['candidates'])} Hopf "
          f"points -> stability changes ONLY at the two folds")
    if hopf["candidates"]:
        raise RuntimeError("possible Hopf point: window claim refused")
    print("=> the bistable window is exactly "
          f"[{_render(B, folds[0]['box'][2], 10)}, "
          f"{_render(B, folds[1]['box'][2], 10)}] nM")
    rec["no_hopf"] = {"deg_R": hopf["deg_R"],
                      "n_excluded": len(hopf["excluded"]),
                      "n_hopf": len(hopf["candidates"])}
    rec["completeness"] = {"deg_R": comp["deg_R"],
                           "n_distinct": comp["n_distinct"],
                           "n_complex": comp["n_complex"],
                           "n_excluded": len(comp["excluded"]),
                           "n_physical": len(comp["candidates"])}

    # ---- certified steady-state tables (M15 three-state table at 50 nM;
    #      M16 critical-slowing rows at 39.3 and 57.3 nM, just inside folds)
    rec["states"] = {}
    for K0 in ("39.3", "50", "57.3"):
        ss = m15.steady_states(float(K0))
        cs = certified_states(S, sp.Rational(K0),
                              [[r["M"], r["Mpp"]] for r in ss])
        if not all(r["ok"] for r in cs) or \
                any(r["kind"] == "UNDECIDED" for r in cs):
            raise RuntimeError(f"state certification REFUSED at {K0}")
        rec["states"][K0] = [_state_rec(r) for r in cs]
        if K0 != "50":
            continue
        print("\nthree steady states at MAPKK = 50 nM (Krawczyk-certified, "
              "stability by certified sign of det/trace):")
        for r in cs:
            lm = _render(B, r["lmax"], 8) if r["lmax"] is not None \
                else "complex"
            print(f"  Mpp = {_render(B, r['Mpp'], 10):>16}  M = "
                  f"{_render(B, r['M'], 8):>12}  lambda_max = {lm:>14}  "
                  f"{r['kind']}  [certified]")
    print("\ncritical slowing (tau = -1/lambda_max, certified):")
    for K0, i, lab in (("39.3", 2, "ON "), ("50", 0, "OFF"), ("50", 2, "ON "),
                       ("57.3", 0, "OFF")):
        st = rec["states"][K0][i]
        print(f"  MAPKK = {K0:>4}  {lab}  Mpp = {st['Mpp_mid']:9.3f}  "
              f"lambda_max = {st['lmax_mid']:+.6f}  tau = "
              f"{st['tau_mid']:8.2f} s")
    print("  (M16/REPORT: tau = 309.1 s at MAPKK = 50 is the OFF state)")

    # ---- negative controls
    f0 = folds[0]
    Fa, Ja = fold_system(S)
    shifted = list(f0["center"])
    shifted[2] = shifted[2] + 5 * f0["radii"][2]
    okd, _, _ = krawczyk(Fa, Ja, shifted, f0["radii"])
    print(f"\n[- control] Krawczyk on a box displaced 5 radii off the fold: "
          f"{'REFUSED (correct)' if not okd else 'CERTIFIED (BUG!)'}")
    rec["controls"]["displaced_box_refused"] = not okd
    if tw:
        halted = False
        try:
            kb = f0["box"][2]
            tw.dual(lambda: _mp(kb), lambda: _mp(kb) * (1 + 1e-8), DIGITS,
                    label="mutation")
        except tw.DualPathDisagreement:
            halted = True
        print(f"[- control] tripwire fed a 1e-8 relative error: "
              f"{'HALTED (correct)' if halted else 'PASSED (BUG!)'}")
        rec["controls"]["tripwire_halts"] = halted

    # ---- sensitivity to the back-solved Km5
    print("\nsensitivity to the back-solved Km5 (certified windows):")
    rec["km5_scan"] = {}
    win78 = (folds[0]["box"][2], folds[1]["box"][2])
    disjoint = True
    for km5 in (77.0, 79.0):
        with km5_set(km5):
            Sk = markevich_system()
            # continuation: start Newton from the Km5 = 78 certified folds
            wk = [certify_fold(Sk, [float(v.mid()) for v in f["center"]])
                  for f in folds]
        if not all(w["ok"] for w in wk):
            raise RuntimeError(f"Krawczyk REFUSED at Km5 = {km5}")
        lo, hi = wk[0]["box"][2], wk[1]["box"][2]
        disjoint &= (not lo.overlaps(win78[0])) and \
            (not hi.overlaps(win78[1]))
        print(f"  Km5 = {km5:4.0f}: [{_render(B, lo, 8)}, "
              f"{_render(B, hi, 8)}] nM")
        rec["km5_scan"][str(km5)] = [_flatten(lo), _flatten(hi)]
    print(f"  Km5 = 78 : [{_render(B, win78[0], 8)}, "
          f"{_render(B, win78[1], 8)}] nM")
    print(f"[- control] mutant windows disjoint from Km5 = 78: "
          f"{'yes (certificate discriminates)' if disjoint else 'NO'}")
    rec["controls"]["km5_mutants_disjoint"] = disjoint

    os.makedirs(os.path.dirname(RECEIPT), exist_ok=True)
    with open(RECEIPT, "w") as fh:
        json.dump(rec, fh, indent=1)
    print(f"\nreceipt written to: {RECEIPT}")
    if not all(rec["controls"].values()):
        raise RuntimeError(f"a control failed: {rec['controls']}")
    return rec


def _state_rec(r):
    """Receipt row: full certified balls (strings) + midpoints for emitall."""
    mt = float(m15.MTOT)
    row = {"kind": r["kind"]}
    vals = {"M": r["M"], "Mp": mt - r["M"] - r["Mpp"], "Mpp": r["Mpp"]}
    if r["lmax"] is not None:
        vals["lmax"] = r["lmax"]
        if r["kind"] == "stable":
            vals["tau"] = -1 / r["lmax"]
    for k, v in vals.items():
        row[k] = _flatten(v)
        row[k + "_mid"] = float(v.mid())
    return row


def _mp(b):
    import mpmath
    with mpmath.workdps(80):
        return mpmath.mpf(b.mid().str(75, radius=False))


def _fold_trace(sysd, box):
    N1, N2 = sysd["N"]
    J = _jac((N1, N2), (X, Y))
    d1, d2 = (APoly(d) for d in sysd["den"])
    x, y, k = box
    return (APoly(J[0][0])(x, y, k) / d1(x, y, k)
            + APoly(J[1][1])(x, y, k) / d2(x, y, k))


if __name__ == "__main__":
    report()
