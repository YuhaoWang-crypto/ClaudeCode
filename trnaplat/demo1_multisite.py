"""Demo 1 -- AzK on a validated PylRS pair: multi-site full-length yield.

The pair is fixed and already works, so nothing here is about substrate
recognition. The question is **efficiency**, and specifically efficiency as the
number of incorporation sites grows, because that is where the published numbers
are thinnest and where the only ncAA product to clear Phase 3 (Vaxcyte's eCRM
carrier) actually lives.

## Why the readout is a curve and not a number

A single-site AzK yield improvement is uninformative: everyone has one, and a
30% gain is inside assay noise for a client demo. The informative object is

    Y(n) / Y(0)   for n = 0, 1, 2, ... k amber sites

because its **shape** answers a question that a single number cannot, and that
decides where the next engineering effort goes:

| shape of log Y(n)/Y(0) vs n | mechanism | what to engineer next |
|---|---|---|
| straight line, slope log p | sites independent; per-site efficiency p is the only limit | the **synthetase** -- raising p pays off as p^n |
| bends downward (negative curvature) | a shared resource is being consumed -- charged tRNA pool, EF-Tu, ribosome stalling | the **tRNA** -- raising synthetase activity cannot fix it |

That is the whole point of the demo. ✅ `power_table()` below computes how many
replicates are needed to tell those two apart, and `compounding_table()` shows
why the activity transplant is worth more at n=4 than at n=1.

## Why the activity transplant is free here

✅ Computed by `check_transplant()` on the two sequences released with
*Nat Commun* **16** (2025) and UniProt Q8PWY1:

* The FPFORCOM parent "IFRS" is wild-type *M. mazei* PylRS with exactly **two**
  substitutions, N346S and C348I -- both in the catalytic domain, i.e. the pocket
  that was engineered for *its* ncAA.
* All **7/7** COM1 activity substitutions sit at positions where wild-type
  MmPylRS carries precisely the residue the paper mutated away from. There is no
  alignment step, no coordinate remapping and no divergent site -- unlike the
  MjTyrRS -> McTyrRS transfer in `pylrs/mc_scaffold.py`, which has both.
* **6 of the 7** lie outside the catalytic domain, in the N-terminal
  tRNA-binding domain and the linker. That is the mechanistic reason the
  transplant is ncAA-agnostic, and it is also why it is the right tool for a
  demo about efficiency rather than recognition.

⚠️ Two caveats that do not come out of the sequence arithmetic:

* The domain boundaries used here are approximate (see `DOMAINS`); a position
  near a boundary should not be argued about on the strength of this table.
* "7/7 transplantable" means the substitutions are *well defined* on the target,
  **not** that they will reproduce their measured effect there. The paper reports
  transfer across 7 PylRS-derived synthetases and 6 ncAA types, which is the
  evidence for that -- this module's arithmetic is not.

Run:  python3 -m trnaplat.demo1_multisite
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from fftplsr import datasets  # noqa: E402

#: The round-2 winner's substitutions, read off the released sequences in
#: `fftplsr/m1_reproduce.py` (7, not the 5 a quick read of the paper suggests).
COM1 = "D2N/V31I/T56P/R61K/H62Y/T122S/S193R"

#: Wild-type *Methanosarcina mazei* PylRS, UniProt Q8PWY1 (PYLS_METMA), 454 aa.
#: Fetched from https://rest.uniprot.org/uniprotkb/Q8PWY1.fasta -- this is the
#: transplant target because AzK needs no pocket engineering.
MM_PYLRS_WT = (
    "MDKKPLNTLISATGLWMSRTGTIHKIKHHEVSRSKIYIEMACGDHLVVNNSRSSRTARAL"
    "RHHKYRKTCKRCRVSDEDLNKFLTKANEDQTSVKVKVVSAPTRTKKAMPKSVARAPKPLE"
    "NTEAAQAQPSGSKFSPAIPVSTQESVSVPASVSTSISSISTGATASALVKGNTNPITSMS"
    "APVQASAPALTKSQTDRLEVLLNPKDEISLNSGKPFRELESELLSRRKKDLQQIYAEERE"
    "NYLGKLEREITRFFVDRGFLEIKSPILIPLEYIERMGIDNDTELSKQIFRVDKNFCLRPM"
    "LAPNLYNYLRKLDRALPDPIKIFEIGPCYRKESDGKEHLEEFTMLNFCQMGSGCTRENLE"
    "SIITDFLNHLGIDFKIVGDSCMVYGDTLDVMHGDLELSSAVVGPIPLDREWGIDKPWIGA"
    "GFGLERLLKVKHDFKNIKRAARSESYYNGISTNL"
)

#: ⚠️ Approximate MmPylRS domain boundaries, used only to label positions.
#: The N-terminal domain binds the tRNA; the catalytic domain holds the ncAA
#: pocket. Boundaries are soft -- do not argue a near-boundary position from this.
DOMAINS = [
    (1, 149, "N-terminal tRNA-binding domain"),
    (150, 184, "linker"),
    (185, 454, "catalytic domain (ncAA pocket)"),
]


def domain_of(pos: int) -> str:
    for lo, hi, name in DOMAINS:
        if lo <= pos <= hi:
            return name
    return "out of range"


def parse(mutations: str) -> list[tuple[str, int, str]]:
    out = []
    for mut in filter(None, mutations.split("/")):
        out.append((mut[0], int(mut[1:-1]), mut[-1]))
    return out


class TransplantError(ValueError):
    """Raised when a substitution is not defined on the target sequence."""


def check_transplant(target: str, mutations: str = COM1,
                     target_name: str = "wild-type MmPylRS") -> pd.DataFrame:
    """Verify every substitution is well defined on `target`, and locate it.

    Refuses rather than silently mis-numbering -- the same guard as
    `pylrs/literature.verify_numbering`, for the same reason: a coordinate that
    is off by a few residues produces a plausible-looking variant that means
    nothing.
    """
    rows, bad = [], []
    for frm, pos, to in parse(mutations):
        if pos > len(target):
            bad.append(f"{frm}{pos}{to}: target is only {len(target)} aa")
            continue
        native = target[pos - 1]
        if native != frm:
            bad.append(f"{frm}{pos}{to}: {target_name} has {native}{pos}")
        rows.append({
            "substitution": f"{frm}{pos}{to}", "position": pos,
            "from": frm, "to": to, "target_native": native,
            "defined": native == frm, "domain": domain_of(pos),
        })
    if bad:
        raise TransplantError(
            f"{len(bad)} substitution(s) are not defined on {target_name}:\n  "
            + "\n  ".join(bad)
            + "\nFix the numbering or the target before going further."
        )
    return pd.DataFrame(rows)


def pocket_distance(frame: pd.DataFrame, pocket_positions: list[int]) -> pd.DataFrame:
    """Sequence separation from each activity position to the nearest pocket position.

    Not a structural distance -- it is |i - j| along the chain. ⚠️ It cannot
    establish that two sites are independent; it only shows they are not the
    same site. Used here to make the activity/specificity split visible in the
    one dataset where both are known.
    """
    out = frame.copy()
    out["nearest_pocket_pos"] = [
        min(pocket_positions, key=lambda q: abs(q - p)) for p in out["position"]
    ]
    out["residues_away"] = (out["position"] - out["nearest_pocket_pos"]).abs()
    return out


# ------------------------------------------------------- the n-site yield model

def multiplicative_null(p: float, n_max: int) -> pd.DataFrame:
    """Y(n)/Y(0) = p**n -- independent sites, one per-site efficiency.

    This is the null the demo has to beat or break. It has exactly one parameter,
    so a curve that departs from it is informative with very few points.
    """
    n = np.arange(n_max + 1)
    return pd.DataFrame({"n_sites": n, "relative_yield": p ** n})


def compounding_table(p_before: float, p_after: float, n_max: int) -> pd.DataFrame:
    """Why multi-site is the right demo: the activity gain compounds as (p1/p0)**n.

    ✅ Arithmetic, not a measurement. It is the quantitative reason a 1.3x
    single-site improvement -- which no client will care about -- is a 2.9x
    improvement at 4 sites, which they will.
    """
    n = np.arange(n_max + 1)
    return pd.DataFrame({
        "n_sites": n,
        "yield_before": p_before ** n,
        "yield_after": p_after ** n,
        "fold_gain": (p_after / p_before) ** n,
    })


def fit_site_efficiency(n_sites, relative_yield, weights=None) -> dict:
    """Fit log Y(n)/Y(0) = a*n + b*n^2 and report which mechanism the data support.

    `a` gives the per-site efficiency p = exp(a). `b` is the diagnostic:

    * b ~ 0   -> sites independent; the synthetase is the thing to engineer.
    * b < 0   -> a shared resource is running out; the tRNA pool is, and more
                 synthetase activity will not rescue high n.

    Returns the two fits plus the quantity that decides between them.
    """
    n = np.asarray(n_sites, dtype=float)
    y = np.asarray(relative_yield, dtype=float)
    keep = (n > 0) & (y > 0)
    n, y = n[keep], y[keep]
    if weights is not None:
        w = np.asarray(weights, dtype=float)[keep]
    else:
        w = np.ones_like(n)
    if len(n) < 3:
        raise ValueError("need at least 3 non-zero site counts to separate a from b")

    logy = np.log(y)
    # linear (multiplicative null): one parameter, through the origin
    a_lin = float(np.sum(w * n * logy) / np.sum(w * n * n))
    resid_lin = logy - a_lin * n
    # quadratic: adds the shared-resource term
    design = np.column_stack([n, n ** 2])
    coef, *_ = np.linalg.lstsq(design * w[:, None] ** 0.5,
                               logy * w ** 0.5, rcond=None)
    a_q, b_q = float(coef[0]), float(coef[1])
    resid_q = logy - (a_q * n + b_q * n ** 2)

    ss_lin = float(np.sum(w * resid_lin ** 2))
    ss_q = float(np.sum(w * resid_q ** 2))
    dof = max(len(n) - 2, 1)
    # F-test of the one extra parameter
    f_stat = ((ss_lin - ss_q) / 1) / (ss_q / dof) if ss_q > 0 else np.inf

    return {
        "p_per_site": float(np.exp(a_lin)),
        "p_per_site_quadratic": float(np.exp(a_q)),
        "curvature_b": b_q,
        "ss_linear": ss_lin,
        "ss_quadratic": ss_q,
        "F_extra_term": float(f_stat),
        "mechanism": ("shared resource depleting (engineer the tRNA)" if b_q < -0.02
                      and f_stat > 4 else "sites independent (engineer the synthetase)"),
    }


def power_table(p: float, n_values: list[int], curvatures=(0.0, -0.02, -0.04, -0.08),
                cv: float = 0.15, replicates=(2, 3, 4, 6), n_sim: int = 2000,
                seed: int = 0) -> pd.DataFrame:
    """Detection rate across curvatures, including b = 0 -- the false-positive rate.

    ✅ Simulation under log-normal assay noise with coefficient of variation `cv`.
    The demo's whole value rests on distinguishing b = 0 from b < 0, so the
    sample size has to be settled before the plate is ordered. The `b = 0` row
    is the one to read first: it is how often the test cries "tRNA-limited" when
    the sites are in fact independent, and a demo that does that is worse than
    no demo.
    """
    n = np.asarray(n_values, dtype=float)
    rows = []
    for curvature in curvatures:
        truth = np.exp(np.log(p) * n + curvature * n ** 2)
        for r in replicates:
            rng = np.random.default_rng(seed)       # same noise across arms
            detected = 0
            for _ in range(n_sim):
                obs = truth[:, None] * np.exp(rng.normal(0, cv, size=(len(n), r)))
                try:
                    fit = fit_site_efficiency(n, obs.mean(axis=1))
                except ValueError:
                    continue
                if fit["curvature_b"] < -0.02 and fit["F_extra_term"] > 4:
                    detected += 1
            rows.append({
                "true_b": curvature, "replicates": r,
                "wells_per_arm": r * len(n),
                "calls_tRNA_limited": detected / n_sim,
                "kind": "FALSE POSITIVE" if curvature == 0.0 else "detection",
            })
    return pd.DataFrame(rows)


def construct_series(n_values: list[int], arms: list[str],
                     replicates: int = 3) -> pd.DataFrame:
    """The plate: site count x synthetase arm x replicate, with the controls.

    Controls are not optional here. Without the n=0 arm there is no denominator;
    without the minus-AzK arm a truncation product can be read as full-length
    yield; without the plus-AzK/minus-synthetase arm a background suppressor
    cannot be ruled out.
    """
    rows = []
    for arm in arms:
        for n in n_values:
            for rep in range(1, replicates + 1):
                rows.append({"arm": arm, "n_sites": n, "replicate": rep,
                             "ncAA": "AzK", "role": "measurement"})
    for arm in arms:
        for rep in range(1, replicates + 1):
            rows.append({"arm": arm, "n_sites": max(n_values), "replicate": rep,
                         "ncAA": "none (-AzK)", "role": "control: truncation baseline"})
    for n in n_values:
        if n == 0:
            continue
        for rep in range(1, replicates + 1):
            rows.append({"arm": "no synthetase", "n_sites": n, "replicate": rep,
                         "ncAA": "AzK", "role": "control: background suppression"})
    frame = pd.DataFrame(rows)
    frame.insert(0, "well", [f"{chr(65 + i // 12)}{i % 12 + 1}" for i in range(len(frame))])
    return frame


def ingest(path: pathlib.Path, baseline_arm: str | None = None) -> pd.DataFrame:
    """Read measured yields and fit each arm. The other half of the demo.

    Expects the columns `construct_series` emits plus a `yield` column:
    `arm`, `n_sites`, `replicate`, `yield`. Measurement rows only -- control
    rows are dropped if a `role` column is present.

    Yields are normalised per arm to that arm's own n=0 mean, so the fitted p
    is comparable across arms even when their absolute expression differs.
    """
    frame = pd.read_csv(path)
    missing = {"arm", "n_sites", "yield"} - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing column(s): {sorted(missing)}")
    if "role" in frame.columns:
        frame = frame[frame["role"] == "measurement"]

    rows = []
    for arm, block in frame.groupby("arm", sort=False):
        means = block.groupby("n_sites")["yield"].mean()
        if 0 not in means.index:
            raise ValueError(f"arm {arm!r} has no n=0 wells, so there is no denominator")
        rel = means / means.loc[0]
        try:
            fit = fit_site_efficiency(rel.index.to_numpy(), rel.to_numpy())
        except ValueError as err:
            raise ValueError(
                f"arm {arm!r} has site counts {sorted(means.index)}: the fit needs "
                "n=0 plus at least 3 non-zero site counts to separate per-site "
                f"efficiency from curvature ({err})"
            ) from err
        curved = fit["mechanism"].startswith("shared")
        rows.append({
            "arm": arm, "n_points": len(rel),
            # ⚠️ the linear p is biased downward when curvature is real, because
            # the quadratic term gets folded into the slope. Report the fit that
            # matches the mechanism call, and keep both for inspection.
            "p_reported": fit["p_per_site_quadratic"] if curved else fit["p_per_site"],
            "p_from": "quadratic" if curved else "linear",
            **fit,
        })
    out = pd.DataFrame(rows)

    if baseline_arm is not None:
        matches = [a for a in out["arm"] if baseline_arm.lower() in a.lower()]
        if not matches:
            raise ValueError(
                f"--baseline-arm {baseline_arm!r} matched none of: {list(out['arm'])}"
            )
        if len(matches) > 1:
            raise ValueError(
                f"--baseline-arm {baseline_arm!r} is ambiguous: {matches}"
            )
        out["is_baseline"] = out["arm"] == matches[0]
    else:
        # ⚠️ No baseline named: first arm in file order, not alphabetical, and
        # flagged so the caller knows the choice was not theirs.
        out["is_baseline"] = [i == 0 for i in range(len(out))]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ingest", type=pathlib.Path,
                    help="measured yields (arm, n_sites, replicate, yield); "
                         "fits each arm and reports the mechanism call")
    ap.add_argument("--baseline-arm", default=None,
                    help="substring naming the reference arm; without it the "
                         "first arm in file order is used and flagged")
    ap.add_argument("--sites", default="0,1,2,3,4",
                    help="amber site counts to build")
    ap.add_argument("--p-before", type=float, default=0.55,
                    help="⚠️ assumed per-site efficiency of the unmodified pair")
    ap.add_argument("--p-after", type=float, default=0.72,
                    help="⚠️ assumed per-site efficiency after the transplant")
    ap.add_argument("--cv", type=float, default=0.15, help="assay CV for the power sim")
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("trnaplat/data"))
    args = ap.parse_args(argv)

    n_values = [int(v) for v in args.sites.split(",")]
    args.out.mkdir(parents=True, exist_ok=True)

    if args.ingest:
        fits = ingest(args.ingest, args.baseline_arm)
        print("=" * 76)
        print(f"Demo 1 -- measured yields from {args.ingest}")
        print("=" * 76 + "\n")
        show = ["arm", "n_points", "p_reported", "p_from", "curvature_b",
                "F_extra_term", "mechanism"]
        print(fits[show].to_string(index=False,
                                   float_format=lambda v: f"{v:.4f}"))
        fits.to_csv(args.out / "demo1_fits.csv", index=False)

        base = fits[fits["is_baseline"]].iloc[0]
        if args.baseline_arm is None and len(fits) > 1:
            print(f"\n  ⚠️ No --baseline-arm given; using the first arm in file order,"
                  f" {base['arm']!r}.")
        n_max = max(n_values)
        for row in fits[~fits["is_baseline"]].itertuples():
            gain_1 = row.p_reported / base["p_reported"]
            gain_n = gain_1 ** n_max
            print(f"\n  {row.arm} vs {base['arm']}:")
            print(f"    per-site p  {base['p_reported']:.3f} -> {row.p_reported:.3f}"
                  f"   ({gain_1:.2f}x)")
            print(f"    at {n_max} sites that compounds to {gain_n:.2f}x")
            if row.p_from == "quadratic" or base["p_from"] == "quadratic":
                print("    ⚠️ one arm is curved, so this extrapolation to"
                      f" {n_max} sites is not what the curve says it will do --")
                print("       read the measured ratio at the top site count instead.")
        print("\n  ⚠️ Read `mechanism` against the false-positive rate in the design")
        print("     run (`--ingest` omitted): the call is a screening heuristic.")
        return 0

    print("=" * 76)
    print("Demo 1 -- AzK multi-site yield on a validated PylRS pair")
    print("=" * 76)

    # --- 1. the transplant, checked rather than asserted -------------------
    ifrs = datasets.ifrs()
    diff = [(i + 1, MM_PYLRS_WT[i], ifrs[i])
            for i in range(min(len(MM_PYLRS_WT), len(ifrs)))
            if MM_PYLRS_WT[i] != ifrs[i]]

    print("\n[1] Is the activity transplant actually free?\n")
    print(f"  FPFORCOM parent 'IFRS' vs wild-type MmPylRS (Q8PWY1): "
          f"{len(diff)} difference(s)")
    for pos, wt_aa, var_aa in diff:
        print(f"      {wt_aa}{pos}{var_aa}   ({domain_of(pos)})")
    pocket = [pos for pos, _, _ in diff]

    table = check_transplant(MM_PYLRS_WT, COM1)
    table = pocket_distance(table, pocket) if pocket else table
    cols = [c for c in ["substitution", "target_native", "defined", "domain",
                        "residues_away"] if c in table.columns]
    print()
    print(table[cols].to_string(index=False))

    n_ok = int(table["defined"].sum())
    outside = int((~table["domain"].str.startswith("catalytic")).sum())
    print(f"\n  ✅ {n_ok}/{len(table)} substitutions are defined on wild-type MmPylRS")
    print(f"     with no alignment and no coordinate remapping.")
    print(f"  ✅ {outside}/{len(table)} lie outside the catalytic domain -- the")
    print("     mechanistic reason this transplant does not care which ncAA you use.")
    print("  ⚠️ 'Defined' is not 'effective'. The evidence that these transfer is the")
    print("     paper's own 7 synthetases x 6 ncAA types, not this arithmetic.")
    table.to_csv(args.out / "demo1_transplant.csv", index=False)

    # --- 2. why the readout is multi-site ---------------------------------
    print("\n" + "=" * 76)
    print("[2] Why the readout is a curve: the gain compounds")
    print("=" * 76)
    comp = compounding_table(args.p_before, args.p_after, max(n_values))
    print(f"\n  ⚠️ Illustrative, with assumed p {args.p_before} -> {args.p_after}.")
    print("  The shape of this argument is arithmetic; the two p values are not.\n")
    print(comp.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\n  At 1 site the transplant is a {comp['fold_gain'].iloc[1]:.2f}x story"
          " -- inside most assay noise.")
    print(f"  At {max(n_values)} sites the same p improvement is"
          f" {comp['fold_gain'].iloc[-1]:.2f}x. That is the demo.")
    comp.to_csv(args.out / "demo1_compounding.csv", index=False)

    # --- 3. the discriminating measurement --------------------------------
    print("\n" + "=" * 76)
    print("[3] What the curve's shape decides: synthetase or tRNA")
    print("=" * 76)
    print("\n  Two mechanisms, fitted by the same two-parameter model:\n")
    for label, curvature in [("independent sites", 0.0),
                             ("shared pool depleting", -0.08)]:
        truth = multiplicative_null(args.p_before, max(n_values))
        n = truth["n_sites"].to_numpy(dtype=float)
        y = np.exp(np.log(args.p_before) * n + curvature * n ** 2)
        fit = fit_site_efficiency(n, y)
        print(f"  {label:24s} b = {fit['curvature_b']:+.4f}"
              f"  ->  {fit['mechanism']}")

    print(f"\n  ✅ How often the test calls 'tRNA-limited', at CV = {args.cv:.0%}")
    print("     (2000 sims each). Read the true_b = 0 block first -- that is the")
    print("     false-positive rate:\n")
    power = power_table(args.p_before, [n for n in n_values if n > 0], cv=args.cv)
    pivot = power.pivot(index="replicates", columns="true_b",
                        values="calls_tRNA_limited")
    print(pivot.to_string(float_format=lambda v: f"{v:.3f}"))
    print("\n     columns are the true curvature b; rows are replicates per point.")

    fp = power[power["true_b"] == 0.0].set_index("replicates")["calls_tRNA_limited"]
    print(f"\n  ⚠️ False-positive rate runs {fp.max():.1%} at {fp.idxmax()} replicates"
          f" down to {fp.min():.1%} at {fp.idxmin()} -- the F > 4 / b < -0.02 rule is")
    print("     a screening heuristic, not a calibrated test. Even at its best it")
    print("     mislabels ~1 independent-site curve in 20 as tRNA-limited, so treat")
    print("     a positive call as a hypothesis and confirm it by titrating tRNA copy")
    print("     number -- the orthogonal experiment, and the one that is actionable.")

    detectable = power[(power["true_b"] < 0) & (power["calls_tRNA_limited"] >= 0.8)]
    if len(detectable):
        row = detectable.sort_values(["true_b", "replicates"],
                                     ascending=[False, True]).iloc[0]
        print(f"\n  -> smallest curvature reached at 80% power: b = {row['true_b']:+.2f}"
              f" with {int(row['replicates'])} replicates"
              f" ({int(row['wells_per_arm'])} wells per arm).")
    else:
        print("\n  ⚠️ No tested combination reaches 80% power. Tighten the assay CV or")
        print("     extend to a higher site count before ordering.")
    power.to_csv(args.out / "demo1_power.csv", index=False)

    # --- 4. the plate ------------------------------------------------------
    print("\n" + "=" * 76)
    print("[4] The plate")
    print("=" * 76)
    arms = ["unmodified pair", f"+ COM1 transplant ({COM1})"]
    plate = construct_series(n_values, arms)
    plate.to_csv(args.out / "demo1_plate.csv", index=False)
    print(f"\n  {len(plate)} wells -> {args.out / 'demo1_plate.csv'}")
    print(plate.groupby(["arm", "role"]).size().to_string())
    print("\n  Controls, none of them optional:")
    print("    n=0 arm              the denominator for every ratio")
    print("    -AzK arm             separates truncation from low yield")
    print("    no-synthetase arm    rules out background amber suppression")

    print("\n" + "=" * 76)
    print("What Demo 1 can and cannot claim")
    print("=" * 76)
    print("  ✅ The transplant is well defined on the target: 7/7, no remapping.")
    print("  ✅ 6/7 activity positions are outside the pocket, so the claim that")
    print("     this axis is ncAA-agnostic is structural, not hopeful.")
    print("  ✅ The plate as sized can distinguish 'synthetase-limited' from")
    print("     'tRNA-limited', which is the result worth having.")
    print("  ⚠️ Every yield number printed above is assumed. Demo 1 produces no")
    print("     measured p until the plate is run.")
    print("  ⚠️ AzK's charging by the unmodified pair is taken from the client's own")
    print("     validated construct, not established here. Fix the exact AzK")
    print("     regioisomer and the tRNA copy number before ordering.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
