"""Step 15 - Build the grafted complexes and redesign the paratope with
ProteinMPNN over four rounds.

THE PROTOCOL
------------
Round 1  free design of every paratope position, T = 0.3, wide sampling.
Round 2  positions that the top decile of round 1 agrees on are frozen; the
         rest are redesigned at T = 0.2. Freezing consensus is what makes the
         rounds do work - resampling the same unconstrained distribution four
         times is four times the compute for one round of information.
Round 3  same, T = 0.1, so the last round converges rather than wanders.
Round 4  no sampling. Every surviving sequence is scored, filtered for
         developability liabilities and scored for humanness against the human
         germline set, then ranked.

WHAT THE SCORE IS AND IS NOT
----------------------------
ProteinMPNN reports the negative log-likelihood of a sequence given a
backbone. That is a measure of how well the sequence fits the fold, NOT of
affinity for the antigen. It is the right objective for "make this paratope
buildable and well-packed against this epitope backbone" and the wrong one
for "rank these by Kd". Nothing downstream of here should be read as a
predicted affinity, and the ranking is therefore reported as a design score.
A fixed-backbone control is included so the numbers have a reference: the
same machinery run on the template's own epitope, where the answer - SM3's
real CDRs - is known.
"""
import json
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS
from step13_template_mining import SOLVENT, get_structure
from step14_graft_and_design import (THREE2ONE, ONE2THREE, chain_residues,
                                     read_atoms_named)

MPNN = Path("/tmp/claude-0/tools/ProteinMPNN")
WORK = Path("/tmp/claude-0/mpnn_work")
BB = ("N", "CA", "C", "O")

# Chimeras to build. Window choice comes from step 14: for the SM3 scaffold
# the top window by buried-position match is TXLNA 514-522, which aligns both
# of MUC1's prolines and puts the S515 phosphosite on a contacted position;
# 511-519 and 503-511 are the runners-up and are carried as alternatives so
# the result does not rest on one threading. For ACC4 only one window is
# backbone-admissible at all.
# Ordered by what the result depends on, because the run is long and a
# truncated run should lose the least important arm: the native control
# decides whether any of these numbers mean anything, the 514-522 graft is
# the lead, ACC4 is the independent scaffold, the last two are alternative
# windows that only refine the lead.
CHIMERAS = [
    {"id": "SM3_native_control", "pdb": "1SM3", "ab": ["H", "L"], "pep": "P",
     "window": None, "txlna": None, "phospho_pos": None},
    {"id": "SM3_514_522", "pdb": "1SM3", "ab": ["H", "L"], "pep": "P",
     "window": "SSPRVTEAP", "txlna": [514, 522], "phospho_pos": 2},
    {"id": "ACC4_502_510", "pdb": "2W65", "ab": ["A", "B"], "pep": "E",
     "window": "ERRPEGPGA", "txlna": [502, 510], "phospho_pos": None},
    {"id": "SM3_511_519", "pdb": "1SM3", "ab": ["H", "L"], "pep": "P",
     "window": "QAPSSPRVT", "txlna": [511, 519], "phospho_pos": 5},
    {"id": "SM3_503_511", "pdb": "1SM3", "ab": ["H", "L"], "pep": "P",
     "window": "RRPEGPGAQ", "txlna": [503, 511], "phospho_pos": None},
]

PARATOPE_CUTOFF = 5.0


def write_pdb(path, chains):
    """chains: [(chain_id, [ {comp, resseq, at:{name:xyz}} ])]"""
    n = 0
    with open(path, "w") as fh:
        for cid, res in chains:
            for i, r in enumerate(res, 1):
                for nm in BB:
                    if nm not in r["at"]:
                        continue
                    x, y, z = r["at"][nm]
                    n += 1
                    # PDB is a fixed-column format and being one column out
                    # is silent: the file still parses, the chain id lands in
                    # the resName field and ProteinMPNN reports a missing
                    # chain. Columns are 13-16 atom, 17 altLoc, 18-20 resName,
                    # 22 chain, 23-26 resSeq, 31-54 coordinates.
                    atom = nm if len(nm) == 4 else f" {nm:<3s}"
                    fh.write(
                        f"ATOM  {n:5d} {atom} {r['comp']:>3s} {cid}"
                        f"{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                        f"  1.00  0.00          {nm[0]:>2s}\n")
            fh.write("TER\n")
        fh.write("END\n")
    return n


def paratope(ab_res, pep_res, cutoff=PARATOPE_CUTOFF):
    """Antibody positions (1-based, renumbered) with an atom near the epitope.

    Contacts are measured against ALL template atoms, not just backbone, so
    the paratope is the real one; the PDB handed to ProteinMPNN is then
    backbone-only because that is all the model reads.
    """
    pep_xyz = [v for r in pep_res for v in r["at"].values()]
    c2 = cutoff ** 2
    sel = []
    for i, r in enumerate(ab_res, 1):
        if r["comp"] == "CYS":
            continue                      # never redesign a disulfide
        hit = any((x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2 <= c2
                  for (x, y, z) in r["at"].values()
                  for (px, py, pz) in pep_xyz)
        if hit:
            sel.append(i)
    return sel


def build(ch):
    atoms = read_atoms_named(get_structure(ch["pdb"]))
    pep = [r for r in chain_residues(atoms, ch["pep"])
           if r["comp"] not in SOLVENT and {"N", "CA", "C"} <= set(r["at"])]
    if ch["window"]:
        assert len(ch["window"]) == len(pep), (ch["id"], len(pep))
        for r, aa in zip(pep, ch["window"]):
            r["comp"] = ONE2THREE[aa]
    ab = {}
    for c in ch["ab"]:
        ab[c] = [r for r in chain_residues(atoms, c)
                 if r["comp"] not in SOLVENT and r["comp"] in THREE2ONE
                 and {"N", "CA", "C"} <= set(r["at"])]
    return ab, pep


def run(cmd, cwd=None):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:])
        print(r.stderr[-3000:])
        raise SystemExit(f"failed: {' '.join(str(c) for c in cmd)}")
    return r.stdout


def mpnn_round(tag, pdb, design_chains, designable, nseq, temp, seed):
    d = WORK / tag
    if d.exists():
        shutil.rmtree(d)
    (d / "pdbs").mkdir(parents=True)
    shutil.copy(pdb, d / "pdbs" / pdb.name)
    parsed, assigned, fixed = d / "p.jsonl", d / "a.jsonl", d / "f.jsonl"
    H = MPNN / "helper_scripts"
    run([sys.executable, str(H / "parse_multiple_chains.py"),
         f"--input_path={d/'pdbs'}", f"--output_path={parsed}"])
    run([sys.executable, str(H / "assign_fixed_chains.py"),
         f"--input_path={parsed}", f"--output_path={assigned}",
         f"--chain_list={' '.join(design_chains)}"])
    run([sys.executable, str(H / "make_fixed_positions_dict.py"),
         f"--input_path={parsed}", f"--output_path={fixed}",
         f"--chain_list={' '.join(design_chains)}",
         "--position_list=" + ", ".join(
             " ".join(str(p) for p in designable[c]) for c in design_chains),
         "--specify_non_fixed"])
    out = d / "out"
    run([sys.executable, str(MPNN / "protein_mpnn_run.py"),
         f"--jsonl_path={parsed}", f"--chain_id_jsonl={assigned}",
         f"--fixed_positions_jsonl={fixed}", f"--out_folder={out}",
         f"--num_seq_per_target={nseq}", f"--sampling_temp={temp}",
         f"--seed={seed}", "--batch_size=10"])
    fa = list((out / "seqs").glob("*.fa"))[0]
    return parse_fasta_scores(fa.read_text())


def parse_fasta_scores(text):
    recs, hdr = [], None
    for line in text.splitlines():
        if line.startswith(">"):
            hdr = line
        elif hdr:
            kv = dict(re.findall(r"(\w+)=([-\d.]+)", hdr))
            recs.append({"score": float(kv.get("score", "nan")),
                         "global_score": float(kv.get("global_score", "nan")),
                         "seq_recovery": float(kv.get("seq_recovery", "nan"))
                         if "seq_recovery" in kv else None,
                         "sample": kv.get("sample"),
                         "chains": line.strip().split("/")})
            recs[-1]["chains"] = line.strip().split("/")
            hdr = None
    return recs


def expand(sel, n):
    """Include each contact's sequence neighbours.

    A 9-residue epitope contacts only 14 SM3 positions at 5 A, which leaves
    ProteinMPNN almost nothing to move. The residues flanking a contact are in
    the same CDR loop and reshape it, so they are designable too; the +-1
    expansion is bounded by the chain and never reaches framework cores.
    """
    out = set()
    for p in sel:
        for q in (p - 1, p, p + 1):
            if 1 <= q <= n:
                out.add(q)
    return sorted(out)


def liabilities(seq):
    pats = {"deamidation_NG_NS": r"N[GS]", "isomerisation_DG_DP": r"D[GP]",
            "Nglyc_sequon": r"N[^P][ST]", "free_Met_oxidation_MM": r"MM",
            "poly_reactivity_RRR": r"[RK]{3,}"}
    return {k: len(re.findall(v, seq)) for k, v in pats.items()
            if re.search(v, seq)}


def consensus_freeze(recs, chains, designable, frac=0.10, thresh=0.8):
    """Positions the top decile agrees on get frozen for the next round."""
    top = sorted(recs, key=lambda r: r["global_score"])[:max(3, int(len(recs) * frac))]
    keep, frozen_aa = {}, {}
    for ci, c in enumerate(chains):
        counts = {p: Counter() for p in designable[c]}
        for r in top:
            s = r["chains"][ci]
            for p in designable[c]:
                counts[p][s[p - 1]] += 1
        still, froz = [], {}
        for p in designable[c]:
            aa, n = counts[p].most_common(1)[0]
            if n / len(top) >= thresh:
                froz[p] = aa
            else:
                still.append(p)
        keep[c], frozen_aa[c] = still, froz
    return keep, frozen_aa


def apply_frozen(pdb_seqs, chains, frozen):
    out = []
    for ci, c in enumerate(chains):
        s = list(pdb_seqs[ci])
        for p, aa in frozen[c].items():
            s[p - 1] = aa
        out.append("".join(s))
    return out


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    from humanness import germline_9mer_coverage, load_germline_kmers
    KMERS, kstats = load_germline_kmers()
    RESULTS.mkdir(exist_ok=True)
    OUT = RESULTS / "step15_mpnn_rounds.json"
    done = {}
    if OUT.exists():
        done = {e["chimera"]["id"]: e for e in json.loads(OUT.read_text())}
    print(f"humanness reference: {kstats}")
    allres = []
    for ch in CHIMERAS:
        if ch["id"] in done:
            print(f"{ch['id']}: already in {OUT.name}, skipping")
            allres.append(done[ch["id"]])
            continue
        ab, pep = build(ch)
        chains = ch["ab"]
        sel = {c: expand(paratope(ab[c], pep), len(ab[c])) for c in chains}
        pdb = WORK / f"{ch['id']}.pdb"
        write_pdb(pdb, [(c, ab[c]) for c in chains] + [(ch["pep"], pep)])
        native = {c: "".join(THREE2ONE[r["comp"]] for r in ab[c]) for c in chains}

        print("=" * 100)
        print(f"{ch['id']}   template {ch['pdb']}   epitope "
              f"{ch['window'] or 'NATIVE ' + ''.join(THREE2ONE[r['comp']] for r in pep)}"
              f"   designable {sum(len(v) for v in sel.values())} positions")
        print("=" * 100)

        designable, frozen = dict(sel), {c: {} for c in chains}
        rounds = []
        for rnd, (nseq, temp) in enumerate(
                [(80, 0.3), (80, 0.2), (80, 0.1)], start=1):
            if not any(designable.values()):
                print(f"  round {rnd}: converged, nothing left to design")
                break
            recs = mpnn_round(f"{ch['id']}_r{rnd}", pdb, chains, designable,
                              nseq, temp, 1000 + rnd)
            recs = [r for r in recs if len(r["chains"]) == len(chains)]
            for r in recs:
                r["chains"] = apply_frozen(r["chains"], chains, frozen)
            best = min(recs, key=lambda r: r["global_score"])
            print(f"  round {rnd}  n={len(recs):3}  T={temp}  "
                  f"designable={sum(len(v) for v in designable.values()):3}  "
                  f"best global_score={best['global_score']:.4f}  "
                  f"median={sorted(r['global_score'] for r in recs)[len(recs)//2]:.4f}")
            rounds.append({"round": rnd, "temp": temp, "n": len(recs),
                           "designable": {c: list(v) for c, v in designable.items()},
                           "best_global_score": best["global_score"],
                           "median_global_score": sorted(
                               r["global_score"] for r in recs)[len(recs) // 2]})
            keep, froz = consensus_freeze(recs, chains, designable)
            for c in chains:
                frozen[c].update(froz[c])
            nfroz = sum(len(v) for v in froz.values())
            print(f"           top-decile consensus froze {nfroz} position(s)"
                  f" -> {sum(len(v) for v in keep.values())} still free")
            designable = keep
            last = recs

        # round 4: no sampling, just judgement
        scored = []
        for r in last:
            seqs = r["chains"]
            lia = {}
            for c, s in zip(chains, seqs):
                for k, v in liabilities(s).items():
                    lia[f"{c}:{k}"] = v
            hcov, _ = germline_9mer_coverage(seqs[0], KMERS)
            mut = {c: [(p, native[c][p - 1], s[p - 1])
                       for p in sel[c] if s[p - 1] != native[c][p - 1]]
                   for c, s in zip(chains, seqs)}
            scored.append({"global_score": r["global_score"],
                           "sequences": dict(zip(chains, seqs)),
                           "n_mutations": sum(len(v) for v in mut.values()),
                           "mutations": {c: [f"{a}{p}{b}" for p, a, b in v]
                                         for c, v in mut.items()},
                           "liabilities": lia,
                           "heavy_germline_9mer_coverage": round(hcov, 4)})
        scored.sort(key=lambda r: (r["global_score"], len(r["liabilities"])))
        clean = [s for s in scored if not s["liabilities"]]
        print(f"  round 4  scored {len(scored)}; {len(clean)} carry no "
              f"sequence liability")
        for s in (clean or scored)[:3]:
            print(f"    score {s['global_score']:.4f}  {s['n_mutations']} muts"
                  f"  Hcov {s['heavy_germline_9mer_coverage']:.3f}"
                  f"  {s['mutations']}")
        allres.append({"chimera": ch, "paratope": {c: list(v) for c, v in sel.items()},
                       "native": native, "rounds": rounds,
                       "n_clean": len(clean),
                       "top": (clean or scored)[:10]})
        # written after EVERY chimera: the run is long enough that a timeout
        # partway through must not throw away the arms that finished.
        OUT.write_text(json.dumps(allres, indent=2))
        print(f"  [saved {len(allres)} chimera(s) to {OUT.name}]")
        print()

    OUT.write_text(json.dumps(allres, indent=2))
    print(f"wrote {OUT} ({len(allres)} chimeras)")


if __name__ == "__main__":
    main()
