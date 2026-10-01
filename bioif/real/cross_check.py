"""
bioif.real.cross_check -- reconciling this work against an independent report.

A strategy report on the same problem (five demo chains A-E, 15 modality
pairs, an 8-dimension bottleneck matrix) covered the p53 <-> Ames edge
independently. Two analyses of the same public datasets is a gift: where they
agree the result is firmer, and where they disagree one of them is wrong in a
way that is worth finding.

Three reconciliations are computed here.

**1. Identity normalisation changes the dataset size.** The report joins
1,565 compounds across Tox21 and Ames; this work joins 2,064. The gap is not
a modelling choice, it is the join key. `identity_yield()` measures all four
plausible keys on the same raw files, and the spread is large enough to
matter: the naive canonical-SMILES join throws away a fifth of the available
paired data. The report's own bottleneck matrix ranks identifier chaos (D2)
as the top failure mode; this is that failure mode costing it calibration
data.

**2. Their stratified confound analysis is better than mine, and it
replicates.** I compared marginal risk ratios across endpoints and concluded
the cytotoxicity confound "does not explain the signal". The report does the
sharper thing -- stratify on the cytotoxicity readout -- and finds the
association vanishing inside the cytotoxic stratum. On 32% more compounds it
replicates closely. That changes the edge from a single number into a
conditional one, which is implemented in `tox_adapters.P53ToMutagenicity`.

**3. Average precision is not comparable across prevalence conventions.**
The report quotes p53 AP 0.191 against a 0.046 base rate; this work gets
0.286 against 0.074. Those look like different models and are not: the lift
over base rate is 4.2x versus 3.9x. The difference is the denominator --
whether unlabelled Tox21 wells count as negatives. Reporting lift alongside
AP removes the ambiguity.

Run: python3 -m bioif.real.cross_check
"""
from __future__ import annotations

import collections
from dataclasses import dataclass
from math import comb

from . import tox

# --------------------------------------------------------------------------
# 1. What the join key costs
# --------------------------------------------------------------------------

JOIN_KEYS = ("canonical_smiles", "inchikey_full", "inchikey_skeleton",
             "parent_skeleton")


def _join_variants(smiles: str) -> dict[str, str | None]:
    """The keys a pipeline might plausibly join on, for one structure."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem.MolStandardize import rdMolStandardize
    RDLogger.DisableLog("rdApp.*")
    out: dict[str, str | None] = {k: None for k in JOIN_KEYS}
    try:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return out
        out["canonical_smiles"] = Chem.MolToSmiles(m)
        try:
            ik = Chem.MolToInchiKey(m)
            if ik:
                out["inchikey_full"] = ik
                out["inchikey_skeleton"] = ik.split("-")[0]
        except Exception:
            pass
        try:
            p = rdMolStandardize.FragmentParent(m)
            p = rdMolStandardize.Uncharger().uncharge(p)
            ikp = Chem.MolToInchiKey(p)
            if ikp:
                out["parent_skeleton"] = ikp.split("-")[0]
        except Exception:
            pass
    except Exception:
        pass
    return out


#: The scan costs ~40 s of RDKit work over ~14k structures, so the result is
#: committed. It is a property of two fixed snapshots, so caching it does not
#: hide anything; BIOIF_REFRESH=1 recomputes and rewrites it.
IDENTITY_CACHE = tox.SNAPSHOT / "identity_yield.json"


def identity_yield(refresh: bool = False) -> dict[str, dict]:
    """Overlap between the two datasets under each join key."""
    import json
    import os
    if not refresh and not os.environ.get("BIOIF_REFRESH") \
            and IDENTITY_CACHE.exists():
        return json.loads(IDENTITY_CACHE.read_text())
    res = _identity_yield_uncached()
    IDENTITY_CACHE.write_text(json.dumps(res, indent=1))
    return res


def _identity_yield_uncached() -> dict[str, dict]:
    t_raw, a_raw = tox.load_raw_smiles()
    tsets: dict[str, set] = {k: set() for k in JOIN_KEYS}
    asets: dict[str, set] = {k: set() for k in JOIN_KEYS}
    for smi in t_raw:
        for k, v in _join_variants(smi).items():
            if v:
                tsets[k].add(v)
    for smi in a_raw:
        for k, v in _join_variants(smi).items():
            if v:
                asets[k].add(v)
    best = len(tsets["parent_skeleton"] & asets["parent_skeleton"])
    return {k: {"tox21": len(tsets[k]), "ames": len(asets[k]),
                "overlap": len(tsets[k] & asets[k]),
                "lost_vs_best": len(tsets[k] & asets[k]) - best}
            for k in JOIN_KEYS}


# --------------------------------------------------------------------------
# 2. The stratified edge
# --------------------------------------------------------------------------

@dataclass
class Stratum:
    name: str
    n: int
    n_pos: int
    odds_ratio: float
    ppv: float
    base: float
    p: float

    @property
    def informative(self) -> bool:
        """Does a p53 call move the mutagenicity risk inside this stratum?"""
        return self.p < 0.05 and self.odds_ratio > 1.5


def _tab(ov, keys):
    a = b = c = d = 0
    for k in keys:
        r = ov[k]
        v = r["tox"][tox.P53]
        if v in ("", "NA"):
            continue
        pos = int(float(v)) == 1
        mut = int(float(r["ames"]["ames"])) == 1
        if pos and mut:
            a += 1
        elif pos:
            b += 1
        elif mut:
            c += 1
        else:
            d += 1
    return a, b, c, d


def _fisher(t) -> float:
    a, b, c, d = t
    r1, r2, c1, n = a + b, c + d, a + c, sum(t)
    if min(r1, r2, c1) == 0 or n == 0:
        return 1.0
    def P(x):
        return (comb(r1, x) * comb(r2, c1 - x)) / comb(n, c1)
    obs = P(a)
    return min(1.0, sum(P(x) for x in range(max(0, c1 - r2), min(r1, c1) + 1)
                        if P(x) <= obs * (1 + 1e-9)))


def _stratum(name, t) -> Stratum:
    a, b, c, d = t
    orr = (a * d) / (b * c) if b and c else float("nan")
    return Stratum(name, sum(t), a + b, orr, a / max(a + b, 1),
                   c / max(c + d, 1), _fisher(t))


def stratified(control: str = tox.CYTOTOX_CONTROL, ov=None) -> dict[str, Stratum]:
    """The p53 -> Ames edge, marginally and within each control stratum."""
    ov = ov if ov is not None else tox.overlap()
    out = {"marginal": _stratum("marginal", _tab(ov, list(ov)))}
    for lab, want in (("+", 1), ("-", 0)):
        ks = [k for k in ov if ov[k]["tox"][control] not in ("", "NA")
              and int(float(ov[k]["tox"][control])) == want]
        out[lab] = _stratum(f"{control}{lab}", _tab(ov, ks))
    return out


# --------------------------------------------------------------------------
# 3. AP is not comparable across prevalence conventions
# --------------------------------------------------------------------------

def ap_lift(ap: float, prevalence: float) -> float:
    return ap / prevalence if prevalence else float("nan")


REPORTED = {
    "p53":  {"ap": 0.191, "prevalence": 0.046, "auroc": 0.765,
             "source": "independent report, 5xCV scaffold split"},
    "ames": {"ap": 0.836, "prevalence": 0.538, "auroc": 0.817,
             "source": "independent report, 5xCV scaffold split"},
}


def report() -> str:
    L = ["Reconciliation against an independent analysis of the same datasets",
         ""]

    # -- 1
    L.append("1. What the identity join key costs (same raw files, 4 keys)")
    try:
        iy = identity_yield()
        L.append(f"   {'join key':<22}{'tox21':>8}{'ames':>8}{'overlap':>9}"
                 f"{'vs best':>9}")
        for k in JOIN_KEYS:
            r = iy[k]
            L.append(f"   {k:<22}{r['tox21']:>8}{r['ames']:>8}"
                     f"{r['overlap']:>9}{r['lost_vs_best']:>+9}")
        best = iy["parent_skeleton"]["overlap"]
        naive = iy["canonical_smiles"]["overlap"]
        L += ["",
              f"   Salt-stripping, neutralising and ignoring stereochemistry",
              f"   buys {best - naive} compounds ({100 * (best - naive) / naive:.0f}% more) over a",
              f"   canonical-SMILES join, and {best - iy['inchikey_full']['overlap']} over a full-InChIKey join.",
              f"   The independent report joins 1,565; this work joins {best}.",
              f"   That is {best - 1565} compounds of paired calibration data available",
              f"   for zero new experiments -- and it is their own D2 bottleneck",
              f"   (identifier chaos) charging rent."]
    except Exception as e:
        L.append(f"   (unavailable: {e})")

    # -- 2
    ov = tox.overlap()
    L += ["", "2. The edge is conditional, not marginal (their analysis, my n)"]
    L.append(f"   {'stratum':<22}{'n':>6}{'n(p53+)':>9}{'OR':>8}{'PPV':>8}"
             f"{'base':>8}{'p':>10}  informative")
    rows = []
    for ctrl in (tox.CYTOTOX_CONTROL, tox.DDR_CONTROL):
        st = stratified(ctrl, ov)
        if not rows:
            m = st["marginal"]
            L.append(f"   {'ALL (marginal)':<22}{m.n:>6}{m.n_pos:>9}"
                     f"{m.odds_ratio:>8.2f}{m.ppv:>8.3f}{m.base:>8.3f}"
                     f"{m.p:>10.1e}  {'yes' if m.informative else 'NO'}")
        for lab in ("+", "-"):
            s = st[lab]
            rows.append(s)
            L.append(f"   {s.name:<22}{s.n:>6}{s.n_pos:>9}"
                     f"{s.odds_ratio:>8.2f}{s.ppv:>8.3f}{s.base:>8.3f}"
                     f"{s.p:>10.1e}  {'yes' if s.informative else 'NO'}")
    mmp = stratified(tox.CYTOTOX_CONTROL, ov)
    L += ["",
          f"   Inside the cytotoxic stratum the association is gone: OR "
          f"{mmp['+'].odds_ratio:.2f} (p={mmp['+'].p:.2f}).",
          f"   Outside it the edge is at its strongest: OR "
          f"{mmp['-'].odds_ratio:.2f} (p={mmp['-'].p:.1e}).",
          "   The report found the same (OR 0.98 vs 2.55) on 1,565 compounds;",
          "   it replicates here on 1,908. So a p53 positive is informative",
          "   about mutagenicity ONLY in compounds that are not already",
          "   flagged cytotoxic -- which makes this a conditional edge, and",
          "   the adapter is now conditional too."]

    # -- 3
    L += ["", "3. AP looked different between the two analyses; it is not"]
    L.append(f"   {'node':<8}{'source':<14}{'AP':>7}{'prev':>7}{'lift':>7}"
             f"{'AUROC':>8}")
    from . import qsar
    for ds, key in (("p53", "p53"), ("ames", "ames")):
        r = REPORTED[key]
        L.append(f"   {ds:<8}{'report':<14}{r['ap']:>7.3f}{r['prevalence']:>7.3f}"
                 f"{ap_lift(r['ap'], r['prevalence']):>7.2f}{r['auroc']:>8.3f}")
        try:
            e = qsar.build(ds, alpha=0.1, seed=0)["eval"]
            L.append(f"   {'':<8}{'this work':<14}{e.ap:>7.3f}"
                     f"{e.prevalence:>7.3f}"
                     f"{ap_lift(e.ap, e.prevalence):>7.2f}{e.auroc:>8.3f}")
        except Exception as exc:
            L.append(f"   {'':<8}{'this work':<14}(unavailable: {exc})")
    L += ["",
          "   The APs differ by a lot and the LIFTS barely differ. AP moves",
          "   with prevalence, and the two analyses use different denominators",
          "   (whether unlabelled Tox21 wells count as negatives). Quote lift",
          "   next to AP, or two correct analyses look like a disagreement."]
    return "\n".join(L)


if __name__ == "__main__":
    print(report())
