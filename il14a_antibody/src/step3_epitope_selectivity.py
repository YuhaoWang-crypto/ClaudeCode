"""Step 3 - Can a paratope be made paralog-SELECTIVE on a given epitope?

THE QUESTION
------------
alpha-taxilin is 52.9% identical to beta-taxilin and 52.2% to gamma-taxilin
over its full length, and both paralogs are broadly expressed. So the binding
question for this target is not only "can we bind it" but "can we bind it
without binding the other two". That is a property of the EPITOPE, and it can
be tested before committing to a paratope.

This script runs a four-round directed optimisation of a CDR3 against a given
epitope and reports the best achievable SPECIFICITY MARGIN:

    margin = complementarity(paratope, TXLNA epitope)
           - max(complementarity(paratope, TXLNB version),
                 complementarity(paratope, TXLNG version))

The scorer (PAIR_SCORE) is an explicit biophysical heuristic, NOT a trained
affinity predictor and NOT a folding engine. Absolute complementarity values
are therefore not interpretable as affinity. The MARGIN is far more robust,
because it compares the same paratope against three near-identical surfaces, so
systematic errors in the heuristic largely cancel.

Used as a screen over epitopes, this answers a real design question: it tells
you which epitope has enough paralog-discriminating chemistry to be worth a
campaign, and it does so without any affinity claim.
"""
import json
import random
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS, gravy, net_charge

SEED = 20260926
AA_DESIGN = "ADEFGHIKLMNPQRSTVWY"     # Cys never introduced into a CDR

POS, NEG = set("RK"), set("DE")
AROM = set("FWYH")
HPHOB = set("AVLIMFWY")
HBOND = set("NQSTYWHRKDE")
HB_DA = set("NQSTY")


def pair_score(p, e):
    """Relative interaction weight for one paratope residue against one epitope
    residue. Explicit chemistry, stated so it can be argued with."""
    s = 0.0
    if (p in POS and e in NEG) or (p in NEG and e in POS):
        s += 2.0                                  # salt bridge
    if (p in POS and e in POS) or (p in NEG and e in NEG):
        s -= 1.5                                  # like-charge repulsion
    if p in AROM and e in AROM:
        s += 1.5                                  # aromatic stacking
    if (p == "M" and e in "FWY") or (e == "M" and p in "FWY"):
        s += 1.3                                  # Met-aromatic
    if (p in POS and e in "FWY") or (e in POS and p in "FWY"):
        s += 1.2                                  # cation-pi
    if p in HPHOB and e in HPHOB:
        s += 1.2                                  # hydrophobic packing
    if p in HB_DA and e in HBOND:
        s += 1.0                                  # side-chain hydrogen bond
    if (p in POS | NEG) and e in set("AVLIMF"):
        s -= 0.5                                  # charge on an apolar face
    return s


LIABILITIES = {
    "free_Cys": r"C", "deamidation": r"N[GS]", "isomerisation": r"D[PGST]",
    "Nglyc_sequon": r"N[^P][ST]", "hydrophobic_run": r"[ILVFWMY]{4,}",
    "poly_charge_run": r"[RK]{3,}|[DE]{3,}",
}
TARGET_NET_CHARGE = (-0.5, 2.5)     # avoids the cationic Axis-1 failure mode
TARGET_GRAVY = (-1.0, 0.4)


def liabilities(seq):
    return [n for n, p in LIABILITIES.items() if re.search(p, seq)]


def developability(cdr3):
    q, gv = net_charge(cdr3), gravy(cdr3)
    qs = 1.0 if TARGET_NET_CHARGE[0] <= q <= TARGET_NET_CHARGE[1] else max(
        0.0, 1.0 - min(abs(q - TARGET_NET_CHARGE[0]),
                       abs(q - TARGET_NET_CHARGE[1])) / 3.0)
    gs = 1.0 if TARGET_GRAVY[0] <= gv <= TARGET_GRAVY[1] else max(
        0.0, 1.0 - min(abs(gv - TARGET_GRAVY[0]),
                       abs(gv - TARGET_GRAVY[1])) / 1.0)
    return qs * gs * max(0.0, 1.0 - 0.34 * len(liabilities(cdr3)))


def build_epitope_spec(lo, hi, sdr_map, hotspots=()):
    """Ordered list of epitope positions with their TXLNA/TXLNB/TXLNG residues.

    sdr_map comes from step1b (real global alignments), so the paralog residues
    are aligned residues, not guesses. Non-discriminating positions are kept in
    the list - they still contribute binding energy and they set the spacing
    between the discriminating ones - but they carry a low weight because they
    cannot distinguish the paralogs.
    """
    spec = []
    for pos in range(lo, hi + 1):
        r = sdr_map.get(pos)
        if not r:
            continue
        discriminating = bool(r["exposed_unique"])
        spec.append({
            "pos": pos, "a": r["aa"],
            "b": r["txlnb"] if discriminating else r["aa"],
            "g": r["txlng"] if discriminating else r["aa"],
            "discriminating": discriminating,
            "weight": (1.8 if pos in hotspots else 1.0) if discriminating else 0.3,
        })
    return spec


def complementarity(cdr3, spec, which):
    """Best one-to-one threading of the CDR3 across the epitope.

    Each CDR3 residue may be used for AT MOST ONE epitope position. The loop is
    threaded across the epitope in register, and every offset and both
    directions (parallel and antiparallel) are tried; the best register wins.

    This replaces an earlier formulation that scored each epitope position
    against max() over the whole CDR3 independently. That version saturated: any
    CDR3 of ~15 residues could supply an ideal partner for every position at
    once, so the score hit a ceiling, the margin became a constant of the
    epitope alone, and the optimiser had nothing to optimise - the same
    score-saturation artefact flagged in the legacy G034 data.
    """
    best = float("-inf")
    n_e, n_p = len(spec), len(cdr3)
    for direction in (1, -1):
        loop = cdr3 if direction == 1 else cdr3[::-1]
        # offset = index in the epitope that loop[0] sits against
        for offset in range(-(n_p - 1), n_e):
            total, contacts = 0.0, 0
            for i, p in enumerate(loop):
                j = offset + i
                if not (0 <= j < n_e):
                    continue
                e = spec[j][which]
                if e is None:
                    continue
                total += spec[j]["weight"] * pair_score(p, e)
                contacts += 1
            if contacts < 4:
                continue            # too few contacts to be an interface
            if total > best:
                best = total
    return best if best > float("-inf") else 0.0


def evaluate(cdr3, spec):
    a = complementarity(cdr3, spec, "a")
    b = complementarity(cdr3, spec, "b")
    g = complementarity(cdr3, spec, "g")
    return {"cdr3": cdr3, "len": len(cdr3),
            "compl_TXLNA": round(a, 3), "compl_TXLNB": round(b, 3),
            "compl_TXLNG": round(g, 3),
            "spec_margin": round(a - max(b, g), 3),
            "net_charge": round(net_charge(cdr3), 2),
            "gravy": round(gravy(cdr3), 3),
            "developability": round(developability(cdr3), 3),
            "liabilities": liabilities(cdr3)}


def mutate(cdr3, rng, n_mut, lo_len, hi_len, allow_indel=True):
    s = list(cdr3)
    for _ in range(n_mut):
        s[rng.randrange(len(s))] = rng.choice(AA_DESIGN)
    out = "".join(s)
    if allow_indel and rng.random() < 0.15:
        if len(out) > lo_len and rng.random() < 0.5:
            i = rng.randrange(len(out))
            out = out[:i] + out[i + 1:]
        elif len(out) < hi_len:
            i = rng.randrange(len(out) + 1)
            out = out[:i] + rng.choice(AA_DESIGN) + out[i:]
    return out


def dist(a, b):
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] != b[i]) + abs(len(a) - len(b))


def diverse_top(records, k, min_dist):
    picked = []
    for r in records:
        if all(dist(r["cdr3"], p["cdr3"]) >= min_dist for p in picked):
            picked.append(r)
        if len(picked) == k:
            break
    return picked


def campaign(spec, seeds, lo_len, hi_len, rng, label):
    """Four rounds: explore -> co-select -> counter-select -> polish."""
    rounds = []

    def pool_from(seqs, n_each, n_mut, indel=True):
        pool = {}
        for s in seqs:
            for _ in range(n_each):
                v = mutate(s, rng, n_mut if isinstance(n_mut, int)
                           else rng.choice(n_mut), lo_len, hi_len, indel)
                if lo_len <= len(v) <= hi_len and "C" not in v:
                    pool[v] = True
        return list(pool)

    # R1 explore - complementarity only
    r1 = [evaluate(v, spec) for v in pool_from(seeds, 3000, [1, 2, 3])]
    r1.sort(key=lambda r: -r["compl_TXLNA"])
    r1_sel = diverse_top(r1, 60, 4)
    rounds.append({"round": "R1_explore", "library": len(r1),
                   "criterion": "compl_TXLNA", "best_margin":
                   max(r["spec_margin"] for r in r1)})

    # R2 co-select - complementarity AND margin
    r2 = [evaluate(v, spec) for v in pool_from([r["cdr3"] for r in r1_sel], 500, [1, 2])]
    r2.sort(key=lambda r: -(r["compl_TXLNA"] + 1.5 * r["spec_margin"]))
    r2_sel = diverse_top(r2, 40, 4)
    rounds.append({"round": "R2_coselect", "library": len(r2),
                   "criterion": "compl_TXLNA + 1.5*margin",
                   "best_margin": max(r["spec_margin"] for r in r2)})

    # R3 counter-select - margin dominant, liability-free required
    r3_all = [evaluate(v, spec) for v in
              pool_from([r["cdr3"] for r in r2_sel], 600, [1, 2], indel=False)]
    r3 = [r for r in r3_all if not r["liabilities"]]
    r3.sort(key=lambda r: -(0.5 * r["compl_TXLNA"] + 3.0 * r["spec_margin"]))
    r3_sel = diverse_top(r3, 25, 5)
    rounds.append({"round": "R3_counterselect", "library": len(r3_all),
                   "liability_free": len(r3),
                   "criterion": "0.5*compl + 3*margin, zero liabilities",
                   "best_margin": max((r["spec_margin"] for r in r3), default=None)})

    # R4 polish - add developability
    r4_all = [evaluate(v, spec) for v in
              pool_from([r["cdr3"] for r in r3_sel], 400, 1, indel=False)]
    r4 = [r for r in r4_all if not r["liabilities"] and r["developability"] >= 0.95]
    r4.sort(key=lambda r: -(0.5 * r["compl_TXLNA"] + 3.0 * r["spec_margin"]
                            + 2.0 * r["developability"]))
    final = diverse_top(r4, 8, 6)
    for i, r in enumerate(final, 1):
        r["candidate_id"] = f"{label}-{i:02d}"
    rounds.append({"round": "R4_polish", "library": len(r4_all),
                   "passed": len(r4),
                   "criterion": "R3 + developability>=0.95",
                   "best_margin": max((r["spec_margin"] for r in r4), default=None)})
    return rounds, final


def main():
    rng = random.Random(SEED)
    sdr_rows = json.loads((RESULTS / "step1b_sdr_map.json").read_text())["per_residue"]
    sdr_map = {r["pos"]: r for r in sdr_rows}
    parents = json.loads((RESULTS / "step2_parents.json").read_text())
    seeds = [r["cdr3"] for r in parents["parent_triage"]
             if r["id"] in parents["selected_parents"]]

    EPITOPES = {
        # the ordered coiled-coil surface picked on structural grounds (step 1b)
        "COILED-COIL aa367-382": {"range": (367, 382), "hotspots": (371, 372, 375),
                                  "len": (10, 18), "label": "CC"},
        # the C-terminal zone the US7622574 1C6/1F2 antibodies were raised against,
        # and the only region of TXLNA with reported functional blocking data
        "C-TERM aa493-523": {"range": (493, 523), "hotspots": (),
                             "len": (12, 20), "label": "CT"},
    }

    out = {"seed": SEED, "epitopes": {}}
    print("=" * 100)
    print("EPITOPE PARALOG-SELECTIVITY TEST  (4 directed rounds per epitope)")
    print("scorer = explicit biophysical heuristic; MARGIN is the interpretable")
    print("quantity, absolute complementarity is NOT an affinity estimate.")
    print("=" * 100)

    for name, cfg in EPITOPES.items():
        lo, hi = cfg["range"]
        spec = build_epitope_spec(lo, hi, sdr_map, cfg["hotspots"])
        base = [evaluate(s, spec) for s in seeds]
        rounds, final = campaign(spec, seeds, cfg["len"][0], cfg["len"][1],
                                 rng, cfg["label"])
        best_parent = max(r["spec_margin"] for r in base)
        best_final = max((r["spec_margin"] for r in final), default=None)
        out["epitopes"][name] = {
            "range": [lo, hi],
            "n_positions": len(spec),
            "n_discriminating_positions": sum(1 for d in spec if d["discriminating"]),
            "discriminating_positions": [
                f"{d['a']}{d['pos']} (B:{d['b']} G:{d['g']})"
                for d in spec if d["discriminating"]],
            "parent_baseline": base, "rounds": rounds, "final_panel": final,
            "best_parent_margin": best_parent, "best_matured_margin": best_final,
        }
        ndisc = sum(1 for d in spec if d["discriminating"])
        print(f"\n### {name}   ({ndisc} discriminating of {len(spec)} positions)")
        for d in spec:
            if d["discriminating"]:
                print(f"      {d['a']}{d['pos']}   beta:{d['b']}  gamma:{d['g']}  w={d['weight']}")
        print(f"  {'round':20} {'library':>8} {'best margin':>12}   criterion")
        print("  " + "-" * 92)
        for rd in rounds:
            bm = "n/a" if rd["best_margin"] is None else f"{rd['best_margin']:.2f}"
            print(f"  {rd['round']:20} {rd['library']:>8} {bm:>12}   {rd['criterion']}")
        print(f"  parent best margin {best_parent:+.2f}  ->  matured best "
              f"{'n/a' if best_final is None else format(best_final, '+.2f')}")
        if final:
            print(f"  {'candidate':10} {'CDR3':22} {'A':7} {'margin':8} {'netQ':6} {'dev':5}")
            for r in final[:5]:
                print(f"  {r['candidate_id']:10} {r['cdr3']:22} {r['compl_TXLNA']:<7.2f} "
                      f"{r['spec_margin']:<+8.2f} {r['net_charge']:<+6.1f} "
                      f"{r['developability']:<5.2f}")

    # --- verdict ---------------------------------------------------------
    print("\n" + "=" * 100)
    print("VERDICT")
    print("=" * 100)
    ranked = sorted(out["epitopes"].items(),
                    key=lambda kv: -(kv[1]["best_matured_margin"] or -99))
    for name, d in ranked:
        bm = d["best_matured_margin"]
        print(f"  {name:24} discriminating positions={d['n_discriminating_positions']:<3} "
              f"best achievable margin={'n/a' if bm is None else format(bm, '+.2f')}")
    out["recommended_epitope"] = ranked[0][0]
    print(f"\n  => higher achievable margin: {ranked[0][0]}")

    (RESULTS / "step3_epitope_selectivity.json").write_text(json.dumps(out, indent=2))
    print("\nwrote results/step3_epitope_selectivity.json")


if __name__ == "__main__":
    main()
