"""Step 18 - Does ProteinMPNN actually respond to WHICH peptide is bound?

Step 15 produced a result that has to be checked before any sequence from it
is called a TXLNA binder. The mutations it proposes for the TXLNA graft are
very nearly the mutations it proposes for the native MUC1 control on the same
scaffold - T28D, M34A, R52G, G100T, G102N, A105G, Y93F, H96T in both. If the
designs do not change when the epitope changes, then they are not designs
against TXLNA; they are generic improvements to SM3's paratope that would
look identical whatever peptide sat in the groove.

This is the same failure mode the Boltz work in this repo hit from the other
direction, where an identical paratope scored differently against two
epitopes and a scrambled peptide outscored the native one. Here the question
is the mirror image: different epitopes, identical output.

FOUR CONDITIONS, ONE BACKBONE
  A  native   the template's own MUC1 peptide          (step 15 control)
  B  graft    the TXLNA window on the same backbone    (step 15 lead)
  C  apo      the peptide chain deleted entirely
  D  scramble a composition-matched shuffle of B

Read it as follows. If B differs from C but not from A or D, ProteinMPNN is
responding to the presence of a peptide backbone and ignoring its sequence -
so the pipeline can place a paratope but cannot specialise it, and every
"TXLNA design" is an SM3 variant. If B differs from A and D as well, the
epitope sequence is doing work and the designs are epitope-directed.
"""
import json
import random
import shutil
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS
from step14_graft_and_design import ONE2THREE, THREE2ONE
from step15_mpnn_rounds import (CHIMERAS, WORK, build, expand, mpnn_round,
                                paratope, write_pdb)

NSEQ, TEMP, SEED = 80, 0.3, 4242


def scramble(seq, seed=20260927):
    rng = random.Random(seed)
    chars = list(seq)
    for _ in range(400):
        rng.shuffle(chars)
        out = "".join(chars)
        same = sum(a == b for a, b in zip(out, seq)) / len(seq)
        if same <= 0.25:
            return out
    return "".join(chars)


def consensus(recs, chains, sel):
    out = {}
    for ci, c in enumerate(chains):
        for p in sel[c]:
            cnt = Counter(r["chains"][ci][p - 1] for r in recs)
            out[f"{c}{p}"] = cnt
    return out


def agreement(a, b):
    """Fraction of designed positions whose consensus residue is the same,
    and the mean total-variation distance between the two distributions."""
    same, tv = 0, []
    for k in a:
        ca, cb = a[k], b[k]
        if ca.most_common(1)[0][0] == cb.most_common(1)[0][0]:
            same += 1
        na, nb = sum(ca.values()) or 1, sum(cb.values()) or 1
        keys = set(ca) | set(cb)
        tv.append(0.5 * sum(abs(ca[x] / na - cb[x] / nb) for x in keys))
    return same / len(a), sum(tv) / len(tv)


def main():
    ch = next(c for c in CHIMERAS if c["id"] == "SM3_514_522")
    ab, pep = build(ch)                       # pep already carries the graft
    chains = ch["ab"]
    sel = {c: expand(paratope(ab[c], pep), len(ab[c])) for c in chains}

    native_ch = next(c for c in CHIMERAS if c["window"] is None)
    ab_n, pep_n = build(native_ch)
    native_seq = "".join(THREE2ONE[r["comp"]] for r in pep_n)
    scr = scramble(ch["window"])
    print(f"epitope A native   {native_seq}")
    print(f"epitope B graft    {ch['window']}   (TXLNA {ch['txlna']})")
    print(f"epitope D scramble {scr}   "
          f"(composition-matched, "
          f"{sum(a==b for a,b in zip(scr, ch['window']))}/{len(scr)} positions "
          f"identical)")
    print(f"condition C: peptide chain removed entirely\n")

    conds = {}
    # A and B: same settings as each other, run fresh here so the comparison
    # is not across different seeds and temperatures from step 15.
    for tag, window, keep_pep in (("A_native", native_seq, True),
                                  ("B_graft", ch["window"], True),
                                  ("C_apo", None, False),
                                  ("D_scramble", scr, True)):
        pep_t = []
        if keep_pep:
            for r, aa in zip(pep, window):
                pep_t.append({"comp": ONE2THREE[aa], "resseq": r["resseq"],
                              "at": dict(r["at"])})
        pdb = WORK / f"resp_{tag}.pdb"
        blocks = [(c, ab[c]) for c in chains]
        if keep_pep:
            blocks.append((ch["pep"], pep_t))
        write_pdb(pdb, blocks)
        recs = mpnn_round(f"resp_{tag}", pdb, chains, sel, NSEQ, TEMP, SEED)
        recs = [r for r in recs if len(r["chains"]) == len(chains)]
        conds[tag] = consensus(recs, chains, sel)
        best = min(r["global_score"] for r in recs)
        print(f"  {tag:11} n={len(recs)}  best global_score={best:.4f}")

    print()
    print("=" * 96)
    print("AGREEMENT BETWEEN CONDITIONS (consensus identity / mean TV distance)")
    print("=" * 96)
    pairs = [("B_graft", "A_native", "graft vs native epitope"),
             ("B_graft", "D_scramble", "graft vs scrambled epitope"),
             ("B_graft", "C_apo", "graft vs NO peptide at all"),
             ("A_native", "C_apo", "native vs NO peptide at all")]
    res = {}
    for x, y, label in pairs:
        ident, tv = agreement(conds[x], conds[y])
        res[f"{x}|{y}"] = {"consensus_identity": round(ident, 4),
                           "mean_tv_distance": round(tv, 4)}
        print(f"  {label:32} identity {ident*100:5.1f}%   TV {tv:.3f}")

    bn = res["B_graft|A_native"]
    bs = res["B_graft|D_scramble"]
    bc = res["B_graft|C_apo"]
    seq_signal = 1 - min(bn["consensus_identity"], bs["consensus_identity"])
    pres_signal = 1 - bc["consensus_identity"]
    print()
    print(f"  sequence signal  (1 - identity to a DIFFERENT epitope) "
          f"{seq_signal:.3f}")
    print(f"  presence signal  (1 - identity to no peptide)          "
          f"{pres_signal:.3f}")
    if pres_signal > 0 and seq_signal / max(pres_signal, 1e-9) < 0.34:
        verdict = ("ProteinMPNN is responding to the PRESENCE of a peptide "
                   "backbone and largely IGNORING its sequence. The grafted "
                   "designs are SM3 variants, not TXLNA-directed paratopes.")
    elif seq_signal < 0.10:
        verdict = ("the epitope sequence changes almost nothing; the designs "
                   "are not epitope-directed")
    else:
        verdict = ("the epitope sequence is doing real work; the designs are "
                   "epitope-directed")
    print(f"\n  -> {verdict}")

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "step18_epitope_responsiveness.json").write_text(json.dumps({
        "epitopes": {"A_native": native_seq, "B_graft": ch["window"],
                     "D_scramble": scr, "C_apo": None},
        "n_designed_positions": sum(len(v) for v in sel.values()),
        "nseq": NSEQ, "temperature": TEMP, "seed": SEED,
        "pairwise": res,
        "sequence_signal": round(seq_signal, 4),
        "presence_signal": round(pres_signal, 4),
        "verdict": verdict,
        "consensus": {k: {p: dict(c) for p, c in v.items()}
                      for k, v in conds.items()},
    }, indent=2))
    print("\nwrote results/step18_epitope_responsiveness.json")


if __name__ == "__main__":
    main()
