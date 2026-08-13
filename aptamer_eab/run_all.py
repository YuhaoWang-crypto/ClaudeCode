"""Run the whole NPY / PP EAB aptamer-library design and write the outputs.

    python3 -m aptamer_eab.run_all

Writes:
    aptamer_eab/output/epitope_scan_{NPY,PP}.csv
    aptamer_eab/output/library_{A,B,C,D,S}_full.csv     (all scored candidates)
    aptamer_eab/output/ORDER_PANEL.csv                  (the curated shortlist)
    aptamer_eab/output/order_specs.csv                  (degenerate oligos)
    aptamer_eab/output/controls.csv
    figures/aptamer_eab_overview.png
"""

import json
import os

import pandas as pd

from . import known
from .targets import (NPY, PP, PYY, FAMILY, identity, diff_positions,
                      epitope_scan, best_epitopes, FOLD_NOTE)
from . import build_library as BL
from .eab_filters import eab_construct, anchor_check

OUT = os.path.join(os.path.dirname(__file__), "output")
FIG = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")


def section(t):
    print("\n" + "=" * 74)
    print(t)
    print("=" * 74)


def run_epitopes(window=12):
    section("1. TARGET / EPITOPE ANALYSIS  [DERIVED from UniProt sequences]")
    print(f"NPY  {NPY}")
    print(f"PP   {PP}")
    print(f"PYY  {PYY}")
    print()
    for a, b in [("NPY", "PP"), ("NPY", "PYY"), ("PP", "PYY")]:
        idp = identity(FAMILY[a], FAMILY[b])
        print(f"  identity {a} vs {b}: {idp:.1%} "
              f"({int(idp*36)}/36 identical)")
    print()

    # positional conservation map
    cons = "".join("*" if NPY[i] == PP[i] == PYY[i] else " " for i in range(36))
    print("  pos      1234567890123456789012345678901234567890"[:46])
    print(f"  NPY      {NPY}")
    print(f"  PP       {PP}")
    print(f"  PYY      {PYY}")
    print(f"  all-same {cons}")
    print(f"\n  fully conserved positions: "
          f"{[i+1 for i,c in enumerate(cons) if c=='*']}")
    print("  -> the conserved block is the C-TERMINUS. An aptamer raised "
          "against\n     the C-terminal segment will cross-react across the "
          "whole family.")

    results = {}
    for name, on, offs in [("NPY", NPY, [PP, PYY]), ("PP", PP, [NPY, PYY])]:
        scan = epitope_scan(on, offs, window=window)
        df = pd.DataFrame(scan)
        df.to_csv(f"{OUT}/epitope_scan_{name}.csv", index=False)
        top = best_epitopes(on, offs, window=window, top=3)
        print(f"\n  Top discriminating {window}-mer windows on {name} "
              f"(unique vs BOTH other family members):")
        for r in top:
            print(f"    {r['start']:>2}-{r['end']:<2}  {r['seq']}  "
                  f"score={r['score']:.2f}  unique={r['n_unique']}/{window} "
                  f"radical={r['n_radical']}")
        worst = sorted(scan, key=lambda r: r["score"])[0]
        print(f"    WORST window (do NOT target): {worst['start']}-{worst['end']}"
              f"  {worst['seq']}  score={worst['score']:.2f}")
        print(f"    fold context: {FOLD_NOTE[name]}")
        results[name] = top
    return results


def run_anchor():
    section("2. ANCHOR CHECK ON THE SCORING FUNCTION  [MEASURED outcomes]")
    print("  Two E-AB outcomes are published for this exact target family.")
    print("  Any scoring function that cannot separate them is not usable:\n")
    passed, rows = anchor_check(verbose=True)
    if not passed:
        raise SystemExit("ANCHOR CHECK FAILED - refusing to emit an order "
                         "panel from a scoring function that cannot reproduce "
                         "the two published outcomes.")
    return rows


def run_libraries():
    section("3. LIBRARY GENERATION + EAB SWITCHABILITY SCORING  [PREDICTED]")
    libs = {}
    print("  A  NPY knowledge-seeded (doped 4.31 core, 15% doping) ...")
    libs["A"] = BL.library_A_npy_seeded()
    print("  B  NPY de novo structured pool ...")
    libs["B"] = BL.library_B_npy_denovo()
    print("  C  PP  cross-seeded (doped 4.31 core, 30% doping) ...")
    libs["C"] = BL.library_C_pp_crossseeded()
    print("  D  PP  de novo structured pool ...")
    libs["D"] = BL.library_D_pp_denovo()
    print("  S  short unframed switches (28-40 nt) ...")
    libs["S"] = BL.library_short()
    print("  G  G-quadruplex track (ranked by G4Hunter, not by eab_score) ...")
    libs["G"] = BL.library_G_g4()

    for k, df in libs.items():
        df.to_csv(f"{OUT}/library_{k}_full.csv", index=False)
        if k == "G":
            print(f"    {k}: {len(df):>5} scored   "
                  f"g4_quality max={df.g4_quality.max():.1f} "
                  f"median={df.g4_quality.median():.1f}  "
                  f"G4Hunter median={df.g4hunter.median():.2f}  "
                  f"(TBA reference = "
                  f"{df[df.arch=='gquad_TBA'].g4_quality.iloc[0]:.1f} / "
                  f"{df[df.arch=='gquad_TBA'].g4hunter.iloc[0]:.2f}; "
                  f"eab_score NOT meaningful here)")
        else:
            print(f"    {k}: {len(df):>5} scored   "
                  f"eab_score max={df.eab_score.max():.1f} "
                  f"median={df.eab_score.median():.1f}")
    return libs


def run_controls(libs):
    section("4. CONTROLS + REFERENCE  [MEASURED sequence, PREDICTED metrics]")
    ctl = BL.controls()
    ctl.to_csv(f"{OUT}/controls.csv", index=False)
    for _, r in ctl.iterrows():
        print(f"  {r['id']:<42} len={r['len']:>3}  MFE={r['mfe']:>7.2f}  "
              f"p_mfe={r['p_mfe']:.3f}  termBP={r['terminal_bp']}  "
              f"eab_score={r['eab_score']:.1f}")
    ref = ctl[ctl.id == "REF_4.31_full80_published_binder"].iloc[0]
    print(f"\n  Reference bar: aptamer 4.31 scores {ref['eab_score']:.1f}. "
          f"It is a REAL binder with a REAL\n  230% E-AB signal change, so this "
          f"is the number every designed candidate\n  must be read against - "
          f"an in-silico score above it is a hypothesis, not an\n  improvement.")
    return ctl, float(ref["eab_score"])


def run_panel(libs, ref_score):
    section("5. CURATED ORDER PANEL  [DESIGN]")
    picks = []
    plan = [("A", 12, "NPY"), ("B", 10, "NPY"),
            ("C", 12, "PP"), ("D", 10, "PP"), ("S", 6, "either"),
            ("G", 8, "either")]
    seen_cores = []
    for key, k, target in plan:
        sub = BL.pick_panel(libs[key], k=k, min_edit=6, exclude=seen_cores)
        sub = sub.copy()
        sub["target"] = target
        picks.append(sub)
        seen_cores.extend([(r.get("core") or r["seq"]) for _, r in sub.iterrows()])
        rng = (f"g4_quality {sub.g4_quality.min():.1f}-{sub.g4_quality.max():.1f}, "
               f"G4Hunter {sub.g4hunter.min():.2f}-{sub.g4hunter.max():.2f}"
               if key == "G" else
               f"score {sub.eab_score.min():.1f}-{sub.eab_score.max():.1f}")
        print(f"  {key}: picked {len(sub):>2} for {target:<6} ({rng})")

    panel = pd.concat(picks, ignore_index=True)
    panel.insert(0, "id", [f"{r.sublib.split('_')[0]}{i+1:02d}"
                           for i, r in enumerate(panel.itertuples())])
    # eab_score is not a meaningful ranking for the G4 track (ViennaRNA cannot
    # see a G-quadruplex), so the comparison is blanked rather than misleading.
    panel["beats_4.31_in_silico"] = [
        "n/a (G4)" if str(r.sublib).startswith("G_quadruplex")
        else str(bool(r.eab_score > ref_score))
        for r in panel.itertuples()]
    cols = ["id", "target", "sublib", "arch", "len", "gc", "mfe", "p_mfe",
            "paired_frac", "terminal_bp", "tail3_free", "self_dimer_dg",
            "eab_score", "arm_arm_bp", "arm_core_bp", "keeps_arm_core_duplex",
            "g4_quality", "g4hunter", "g4_tracts",
            "longest_g_run", "beats_4.31_in_silico", "n_mut", "core",
            "seq", "struct"]
    panel = panel[[c for c in cols if c in panel.columns]]
    panel.to_csv(f"{OUT}/ORDER_PANEL.csv", index=False)
    print(f"\n  ORDER_PANEL.csv: {len(panel)} sequences "
          f"({(panel.target=='NPY').sum()} NPY, {(panel.target=='PP').sum()} PP, "
          f"{(panel.target=='either').sum()} target-agnostic scaffolds)")
    nb = int((panel["beats_4.31_in_silico"] == "True").sum())
    print(f"  {nb} of them out-score the 4.31 reference on the heuristic. That "
          f"is NOT evidence\n  of better binding - the heuristic was anchored "
          f"on 4.31's own structure.")

    framed = panel[panel.arm_core_bp.notna()]
    if len(framed):
        keep = framed.keeps_arm_core_duplex.fillna(False).astype(bool)
        print(f"\n  5'arm<->core duplex (4.31 has 5 bp at arm positions "
              f"14,15,18,19,20):\n    {int(keep.sum())}/{len(framed)} framed "
              f"panel members keep >=4 bp. The orthogonal contact-probability\n"
              f"    run puts ~1/3 of the predicted interface on the 5' arm "
              f"(strongest contacts\n    T14/C15/A16), so candidates that "
              f"break this duplex are deprioritised.")
        for sl in ["A_NPY_seeded", "C_PP_crossseeded"]:
            sub = framed[framed.sublib == sl]
            if len(sub):
                k = sub.keeps_arm_core_duplex.fillna(False).astype(bool)
                print(f"      {sl:<18} {int(k.sum())}/{len(sub)} keep it "
                      f"(median {sub.arm_core_bp.median():.0f} bp)")

    print("\n  Top 5 per target (id | eab_score | len | MFE | termBP):")
    for t in ["NPY", "PP"]:
        print(f"   {t}:")
        for _, r in panel[panel.target == t].nlargest(5, "eab_score").iterrows():
            print(f"     {r['id']:<6} {r['eab_score']:>5.1f}  {r['len']:>3} nt  "
                  f"{r['mfe']:>7.2f}  {r['terminal_bp']}")
    return panel


def run_specs(panel):
    section("6. ORDERABLE SPECS  [DESIGN]")
    specs = BL.order_specs()
    pd.DataFrame(specs).to_csv(f"{OUT}/order_specs.csv", index=False)
    for s in specs:
        print(f"  {s['name']}\n    {s['spec']}\n    use: {s['use']}\n"
              f"    scale: {s['scale']}")
    print("\n  Example EAB construct for the top NPY panel member:")
    top = panel[panel.target == "NPY"].nlargest(1, "eab_score").iloc[0]
    c = eab_construct(top["seq"], top["id"])
    for k, v in c.items():
        print(f"    {k:<11}: {v}")


def make_figure(libs, panel, ctl, ref_score):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import sys
    sys.path.insert(0, "/root/.claude/skills/synced/aptamer-design")
    import RNA

    os.makedirs(FIG, exist_ok=True)
    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.28)

    # (a) family conservation
    ax = fig.add_subplot(gs[0, 0])
    same_pp = [1 if NPY[i] == PP[i] else 0 for i in range(36)]
    same_pyy = [1 if NPY[i] == PYY[i] else 0 for i in range(36)]
    x = range(1, 37)
    ax.bar(x, same_pp, color="#c44e52", label="NPY = PP", width=0.85)
    ax.bar(x, [-v for v in same_pyy], color="#4c72b0", label="NPY = PYY",
           width=0.85)
    ax.axvspan(24.5, 36.5, color="grey", alpha=0.18)
    ax.text(30.5, 0.55, "conserved\nC-term\n(AVOID)", ha="center", fontsize=8)
    ax.set_xlabel("NPY residue"); ax.set_yticks([])
    ax.set_title("(a) family identity map", fontsize=10)
    ax.legend(fontsize=7, loc="lower left")

    # (b) epitope discrimination score
    ax = fig.add_subplot(gs[0, 1])
    for nm, on, offs, col in [("NPY", NPY, [PP, PYY], "#c44e52"),
                              ("PP", PP, [NPY, PYY], "#55a868")]:
        sc = epitope_scan(on, offs, window=12)
        ax.plot([r["start"] for r in sc], [r["score"] for r in sc],
                "-o", ms=3, color=col, label=nm)
    ax.set_xlabel("window start residue")
    ax.set_ylabel("discrimination score")
    ax.set_title("(b) 12-mer epitope scan", fontsize=10)
    ax.legend(fontsize=8)

    # (c) score distributions
    ax = fig.add_subplot(gs[0, 2])
    names = {"A": "A NPY seeded", "B": "B NPY de novo", "C": "C PP seeded",
             "D": "D PP de novo", "S": "S short"}   # G4 excluded: not comparable
    ax.boxplot([libs[k].eab_score for k in names], tick_labels=list(names.values()),
               showfliers=False)
    ax.axhline(ref_score, color="#c44e52", ls="--", lw=1.4,
               label=f"aptamer 4.31 = {ref_score:.0f}")
    ax.set_ylabel("EAB switchability score")
    ax.set_title("(c) library score distributions", fontsize=10)
    ax.tick_params(axis="x", rotation=30, labelsize=7)
    ax.legend(fontsize=7)

    # (d) switch landscape
    ax = fig.add_subplot(gs[1, 0])
    for k, col in zip("ABCDS", ["#c44e52", "#dd8452", "#55a868", "#4c72b0",
                                "#8172b3"]):
        d = libs[k]
        ax.scatter(d.mfe, d.p_mfe, s=5, alpha=0.25, color=col, label=k)
    ax.scatter([ctl[ctl.id.str.startswith("REF")].mfe.iloc[0]],
               [ctl[ctl.id.str.startswith("REF")].p_mfe.iloc[0]],
               s=140, marker="*", color="black", zorder=5, label="4.31")
    ax.set_xlabel("MFE (kcal/mol, DNA params, 37 C)")
    ax.set_ylabel("P(MFE structure)")
    ax.set_title("(d) stability vs plasticity", fontsize=10)
    ax.legend(fontsize=7, ncol=2)

    # (e) terminal stem = the switch element
    ax = fig.add_subplot(gs[1, 1])
    allc = pd.concat([libs[k] for k in "ABCDS"])
    ax.hexbin(allc.terminal_bp, allc.eab_score, gridsize=22, cmap="viridis",
              mincnt=1)
    ax.set_xlabel("terminal (5'-3') base pairs")
    ax.set_ylabel("EAB switchability score")
    ax.set_title("(e) the 5'anchor-3'MB stem drives the score", fontsize=10)

    # (f) MFE structure of the top NPY pick
    ax = fig.add_subplot(gs[1, 2])
    top = panel[panel.target == "NPY"].nlargest(1, "eab_score").iloc[0]
    coords = RNA.simple_xy_coordinates(top["struct"])
    # simple_xy_coordinates returns n+1 points (a trailing sentinel); keep n.
    xs = [c.X for c in coords][:len(top["struct"])]
    ys = [c.Y for c in coords][:len(top["struct"])]
    ax.plot(xs, ys, "-", color="0.6", lw=1)
    ax.scatter(xs, ys, c=["#c44e52" if ch != "." else "#4c72b0"
                          for ch in top["struct"]], s=9, zorder=3)
    ax.scatter([xs[0]], [ys[0]], s=110, marker="s", color="gold",
               edgecolor="k", zorder=4)
    ax.scatter([xs[-1]], [ys[-1]], s=110, marker="o", color="#1f77b4",
               edgecolor="k", zorder=4)
    ax.text(xs[0], ys[0], "  5'-SH", fontsize=7, va="center")
    ax.text(xs[-1], ys[-1], "  3'-MB", fontsize=7, va="center")
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"(f) MFE fold, top NPY pick {top['id']}\n"
                 f"MFE {top['mfe']:.1f} kcal/mol, {top['terminal_bp']} terminal bp",
                 fontsize=9)

    fig.suptitle("ssDNA aptamer library design for NPY / PP EAB sensors  "
                 "- in-silico prioritisation only, nothing experimentally tested",
                 fontsize=12)
    out = f"{FIG}/aptamer_eab_overview.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\n  figure -> {out}")


def main():
    os.makedirs(OUT, exist_ok=True)
    eps = run_epitopes()
    run_anchor()
    libs = run_libraries()
    ctl, ref_score = run_controls(libs)
    panel = run_panel(libs, ref_score)
    run_specs(panel)
    make_figure(libs, panel, ctl, ref_score)

    with open(f"{OUT}/summary.json", "w") as fh:
        json.dump({
            "reference_4.31_eab_score": ref_score,
            "library_sizes": {k: int(len(v)) for k, v in libs.items()},
            "g4_track_ranked_by": "g4_quality, anchored on TBA "
                                  "(ViennaRNA cannot model G4)",
            "panel_size": int(len(panel)),
            "top_epitopes": eps,
        }, fh, indent=2)

    section("DONE - MANDATORY CAVEATS")
    print("""  * Everything designed here is IN SILICO ONLY. No sequence in the
    panel has been shown to bind NPY or PP.
  * The EAB switchability score is an engineering heuristic over ViennaRNA
    secondary structure. It is NOT an affinity prediction and NOT a predicted
    signal change. Only aptamer 4.31 has measured numbers, and they come from
    the literature, not from this pipeline.
  * ViennaRNA DNA parameters (Mathews2004) do not model Mg2+/Na+ properly and
    do not model G-quadruplexes or any tertiary structure at all. G4 candidates
    in particular are under-modelled - fold them experimentally (CD, thermal
    melt) before trusting anything about them.
  * NPY, PP and PYY are 44-70% identical and share a near-identical C-terminus.
    Family specificity CANNOT be assumed; it has to be produced by explicit
    counter-selection and then measured.
  * A library is a starting pool for SELEX / screening, not a binder.""")


if __name__ == "__main__":
    main()
