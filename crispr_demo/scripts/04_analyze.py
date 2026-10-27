#!/usr/bin/env python3
"""
04_analyze.py - turn the pipeline output into the numbers that go in the report.

1. scan the reference for sgRNA near-cognate sites (a Cas-OFFinder style
   mismatch search) and rank them
2. classify every EDIT-only variant as on-target / off-target / background by
   its distance to a predicted cut site
3. re-estimate indel VAF from CIGAR strings (pileup AD under-counts long
   insertions because the reads get soft-clipped)
4. bin coverage and compute the EDIT/WT ratio, which exposes the large deletion
5. score everything against the ground truth

Writes results/*.tsv and results/analysis.json
"""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "ref"
WORK = ROOT / "work"
RES = ROOT / "results"
RES.mkdir(parents=True, exist_ok=True)

CHROM = "DEMO_chr1"
SAMPLES = ["WT_parent", "EDIT_pool"]
COMP = str.maketrans("ACGTN", "TGCAN")
SEED_LEN = 12          # PAM-proximal seed; mismatches here hurt the most
CUT_WINDOW = 12        # a variant this close to a predicted cut is attributed to it
DEPLETION_CUTOFF = 0.85  # below this, part of the population has lost the locus


def rc(s):
    return s.translate(COMP)[::-1]


def sh(cmd):
    return subprocess.run(cmd, shell=True, cwd=ROOT, capture_output=True,
                          text=True, check=True).stdout


def read_fasta(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


# ---------------------------------------------------------------------------
# 1. off-target site prediction
# ---------------------------------------------------------------------------
def scan_offtargets(genome, guide, max_mm=5):
    """Brute-force scan both strands for protospacer+PAM matches.

    Returns every site with <= max_mm mismatches to the 20 nt protospacer and a
    canonical NGG or near-canonical NAG PAM.
    """
    hits = []
    n = len(guide)
    for strand in ("+", "-"):
        seq = genome if strand == "+" else rc(genome)
        for i in range(len(seq) - n - 3):
            proto = seq[i:i + n]
            pam = seq[i + n:i + n + 3]
            if pam[1:] not in ("GG", "AG"):
                continue
            mm, seed_mm = 0, 0
            for k in range(n):
                if proto[k] != guide[k]:
                    mm += 1
                    if k >= n - SEED_LEN:
                        seed_mm += 1
                    if mm > max_mm:
                        break
            if mm > max_mm:
                continue
            # convert back to +strand 1-based coordinates
            if strand == "+":
                start = i + 1
                cut = start + 16
            else:
                start = len(genome) - (i + n + 3) + 1 + 3
                cut = start + (n - 16) - 1
            hits.append(dict(strand=strand, start=start, cut_site=cut,
                             protospacer=proto, pam=pam, mismatches=mm,
                             seed_mismatches=seed_mm,
                             canonical_pam=pam[1:] == "GG"))
    return hits


def risk_tier(hit):
    """A transparent heuristic ranking. NOT a validated CFD/MIT score."""
    if hit["mismatches"] == 0 and hit["canonical_pam"]:
        return "on-target"
    score = hit["mismatches"] + 1.5 * hit["seed_mismatches"]
    if not hit["canonical_pam"]:
        score += 2.0
    if score <= 3.0:
        return "high"
    if score <= 6.0:
        return "medium"
    return "low"


# ---------------------------------------------------------------------------
# 2. CIGAR-aware allele counting
# ---------------------------------------------------------------------------
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def wilson_ci(k, n, z=1.96):
    """95% Wilson score interval for a binomial proportion.

    Reported alongside every VAF: at realistic depth a sub-clonal allele
    frequency carries real uncertainty and a bare point estimate hides it.
    """
    if n == 0:
        return (None, None)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z / d * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5)
    return (round(max(0.0, centre - half), 4), round(min(1.0, centre + half), 4))


def count_allele_support(bam, pos, event_kind, event_len, anchor=20):
    """Count reads that genuinely span the event and whether they carry it.

    pileup AD treats a soft-clipped read as reference support, which badly
    under-estimates the VAF of a long insertion.  Requiring a read to span the
    whole event interval plus `anchor` bases on each side removes that bias.
    The window is centred on the event itself, not on the VCF anchor base, so
    that reference- and alt-carrying reads get the same landing zone.
    """
    ev_start = pos + 1                                   # first affected base
    ev_end = pos + event_len if event_kind == "del" else pos
    lo, hi = max(1, ev_start - anchor), ev_end + anchor
    region = "%s:%d-%d" % (CHROM, lo, hi)
    out = sh("samtools view -q 20 -F 0x900 %s %s" % (bam, region))
    alt = ref = 0
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        rpos, cigar = int(f[3]), f[5]
        ops = CIGAR_RE.findall(cigar)
        refp = rpos
        carries = False
        ref_start, ref_end = rpos, rpos - 1
        for ln, op in ops:
            ln = int(ln)
            if op in "MD=XN":
                if op == "D" and event_kind == "del" and ln == event_len \
                        and abs(refp - (pos + 1)) <= 3:
                    carries = True
                refp += ln
                ref_end = refp - 1
            elif op == "I":
                if event_kind == "ins" and ln == event_len and abs(refp - (pos + 1)) <= 3:
                    carries = True
        # must truly span the window, not just touch it
        if ref_start > lo or ref_end < hi:
            continue
        if carries:
            alt += 1
        else:
            ref += 1
    total = alt + ref
    lo_ci, hi_ci = wilson_ci(alt, total)
    return dict(alt=alt, ref=ref, depth=total,
                vaf=(alt / total) if total else None,
                ci95=[lo_ci, hi_ci])


# ---------------------------------------------------------------------------
# 3. coverage
# ---------------------------------------------------------------------------
def load_depth(depth_file, genome_len=60000):
    d = [0] * (genome_len + 1)
    for line in open(depth_file):
        _, p, v = line.split("\t")
        d[int(p)] = int(v)
    return d


def population_frequency(vaf_local, pos, depth_edit, depth_wt, mean_edit, mean_wt,
                         event_len=0, flank=300):
    """Convert a locally-measured VAF into a population frequency.

    A VAF is only ever measured among the reads that are actually there.  When
    part of the population has lost the whole region (here: the 1150 bp
    on-target deletion), those cells contribute no reads at all, so a variant
    carried by *every remaining* cell reads as VAF 1.0 even though it is only
    in ~55% of the population.  Scaling by the local coverage - normalised
    against the same window in the un-edited control, which cancels GC and
    mappability effects - recovers the population figure.
    """
    # two flanks that exclude the event itself, so the deleted/inserted bases
    # do not drag the local coverage estimate down
    gmax = len(depth_edit) - 1
    spans = [(max(1, pos - flank), max(1, pos - 1)),
             (min(gmax, pos + event_len + 1), min(gmax, pos + event_len + flank))]
    de_sum = dw_sum = n = 0
    for a, b in spans:
        if b < a:
            continue
        de_sum += sum(depth_edit[a:b + 1])
        dw_sum += sum(depth_wt[a:b + 1])
        n += b - a + 1
    if n == 0 or dw_sum <= 0 or mean_wt <= 0 or mean_edit <= 0:
        return None, None
    de, dw = de_sum / n, dw_sum / n
    ratio = (de / mean_edit) / (dw / mean_wt)
    if vaf_local is None:
        return None, round(ratio, 4)
    # Only rescale when coverage is genuinely depleted, i.e. a competing allele
    # has removed the locus from part of the population.  Applying the
    # correction everywhere would just feed ordinary coverage noise into the
    # frequency estimate.
    if ratio >= DEPLETION_CUTOFF:
        return round(vaf_local, 4), round(ratio, 4)
    return round(min(1.0, vaf_local * ratio), 4), round(ratio, 4)


def binned_coverage(depth_file, bin_size=200, genome_len=60000):
    sums = [0] * (genome_len // bin_size + 1)
    for line in open(depth_file):
        _, p, d = line.split("\t")
        sums[(int(p) - 1) // bin_size] += int(d)
    return [round(s / bin_size, 2) for s in sums[:genome_len // bin_size]]


# ---------------------------------------------------------------------------
def main():
    genome = read_fasta(REF / "DEMO_chr1.fasta")
    design = json.load(open(REF / "design_notes.json"))
    guides = design["guides"]

    # ---- off-target prediction -------------------------------------------
    ot_rows = []
    for g in guides:
        for hit in scan_offtargets(genome, g["protospacer"]):
            hit["guide"] = g["id"]
            hit["tier"] = risk_tier(hit)
            ot_rows.append(hit)
    ot_rows.sort(key=lambda h: (h["mismatches"], h["seed_mismatches"]))

    with open(RES / "offtarget_sites.tsv", "w") as fh:
        cols = ["guide", "tier", "strand", "start", "cut_site", "protospacer",
                "pam", "mismatches", "seed_mismatches", "canonical_pam"]
        fh.write("\t".join(cols) + "\n")
        for h in ot_rows:
            fh.write("\t".join(str(h[c]) for c in cols) + "\n")

    # per-base depth, needed to convert local VAF into a population frequency
    dep = {s: load_depth(WORK / "cov" / ("%s.depth.txt" % s)) for s in SAMPLES}
    mean_d = {s: sum(dep[s][1:]) / (len(dep[s]) - 1) for s in SAMPLES}

    # ---- EDIT-only small variants ----------------------------------------
    ann = sh("bcftools query -f '%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/DP\\t"
             "%INFO/MQ\\t%INFO/ANN[\\t%GT:%AD]\\n' work/vcf/joint.ann.vcf.gz")
    variants = []
    for line in ann.strip().splitlines():
        f = line.split("\t")
        pos, ref, alt, qual, dp, mq, annf = f[0], f[1], f[2], f[3], f[4], f[5], f[6]
        gt_wt, gt_ed = f[7], f[8]
        a = annf.split(",")[0].split("|")
        effect, impact, gene, hgvs_c, hgvs_p = a[1], a[2], a[3], a[9], a[10]

        def ad(field):
            g, adv = field.split(":")
            nums = [int(x) for x in adv.split(",")]
            return g, nums

        gwt, adwt = ad(gt_wt)
        ged, aded = ad(gt_ed)

        if len(ref) > len(alt):
            kind, elen = "del", len(ref) - len(alt)
        elif len(alt) > len(ref):
            kind, elen = "ins", len(alt) - len(ref)
        else:
            kind, elen = "snv", 0

        pos = int(pos)
        rec = dict(pos=pos, ref=ref, alt=alt, qual=float(qual), dp=int(dp),
                   mq=float(mq), kind=kind, event_len=elen, effect=effect,
                   impact=impact, gene=gene, hgvs_c=hgvs_c, hgvs_p=hgvs_p,
                   ad_wt=adwt, ad_edit=aded,
                   vaf_wt_pileup=round(adwt[1] / sum(adwt), 4) if sum(adwt) else None,
                   vaf_edit_pileup=round(aded[1] / sum(aded), 4) if sum(aded) else None)

        # CIGAR-corrected VAF for indels
        if kind in ("ins", "del"):
            for s, key in ((("WT_parent"), "wt"), (("EDIT_pool"), "edit")):
                c = count_allele_support("work/align/%s.bam" % s, pos, kind, elen)
                rec["span_%s" % key] = c
            rec["vaf_edit_span"] = (round(rec["span_edit"]["vaf"], 4)
                                    if rec["span_edit"]["vaf"] is not None else None)
            rec["vaf_wt_span"] = (round(rec["span_wt"]["vaf"], 4)
                                  if rec["span_wt"]["vaf"] is not None else None)
            rec["ci95_edit"] = rec["span_edit"]["ci95"]
            rec["informative_depth_edit"] = rec["span_edit"]["depth"]
        else:
            rec["vaf_edit_span"] = rec["vaf_edit_pileup"]
            rec["vaf_wt_span"] = rec["vaf_wt_pileup"]
            rec["ci95_edit"] = wilson_ci(aded[1], sum(aded))
            rec["informative_depth_edit"] = sum(aded)

        popf, covratio = population_frequency(
            rec["vaf_edit_span"], pos, dep["EDIT_pool"], dep["WT_parent"],
            mean_d["EDIT_pool"], mean_d["WT_parent"],
            event_len=(elen if kind == "del" else 0))
        rec["local_cov_ratio"] = covratio
        rec["population_freq"] = popf
        rec["competing_allele"] = (covratio is not None
                                   and covratio < DEPLETION_CUTOFF)
        # the interval has to ride the same rescaling as the point estimate,
        # otherwise the two are plotted on different scales
        ci = rec.get("ci95_edit") or [None, None]
        if ci[0] is None or covratio is None or covratio >= DEPLETION_CUTOFF:
            rec["ci95_population"] = ci
        else:
            rec["ci95_population"] = [round(min(1.0, ci[0] * covratio), 4),
                                      round(min(1.0, ci[1] * covratio), 4)]

        rec["present_wt"] = adwt[1] >= 3 and (adwt[1] / max(1, sum(adwt))) >= 0.05
        rec["present_edit"] = aded[1] >= 3 and (aded[1] / max(1, sum(aded))) >= 0.05
        rec["category"] = "shared (parental)" if rec["present_wt"] else "EDIT-specific"

        # attribute to a cut site
        best = None
        for h in ot_rows:
            d = abs(h["cut_site"] - pos)
            if d <= CUT_WINDOW and (best is None or d < best[0]):
                best = (d, h)
        if best:
            d, h = best
            rec["attributed_to"] = h["guide"]
            rec["cut_distance"] = d
            rec["site_mismatches"] = h["mismatches"]
            rec["origin"] = "on-target" if h["mismatches"] == 0 else "off-target"
        else:
            rec["attributed_to"] = None
            rec["cut_distance"] = None
            rec["site_mismatches"] = None
            rec["origin"] = "background" if rec["category"] == "EDIT-specific" else "-"
        variants.append(rec)

    # ---- structural variants ---------------------------------------------
    sv_raw = sh("bcftools query -f '%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t"
                "%FILTER\\t%INFO/PE\\t%INFO/SR[\\t%GT:%RV:%DV:%RR:%DR]\\n' "
                "work/sv/joint.bcf")
    svs = []
    for line in sv_raw.strip().splitlines():
        f = line.split("\t")
        chrom, start, end, svtype, filt, pe, sr = f[0], int(f[1]), int(f[2]), f[3], f[4], f[5], f[6]
        gts = {}
        for s, field in zip(SAMPLES, f[7:]):
            gt, rv, dv, rr, dr = field.split(":")
            rv, dv, rr, dr = int(rv), int(dv), int(rr), int(dr)
            gts[s] = dict(gt=gt, split_alt=rv, pair_alt=dv, split_ref=rr, pair_ref=dr,
                          vaf=round(rv / (rv + rr), 4) if (rv + rr) else 0.0,
                          ci95=list(wilson_ci(rv, rv + rr)),
                          informative_depth=rv + rr)
        rec = dict(chrom=chrom, start=start, end=end, svtype=svtype,
                   length=end - start, filter=filt, pe=int(pe), sr=int(sr),
                   genotypes=gts)
        cut = None
        for h in ot_rows:
            if h["mismatches"] == 0 and start <= h["cut_site"] <= end:
                cut = h
        rec["spans_cut_site"] = cut["cut_site"] if cut else None
        rec["attributed_to"] = cut["guide"] if cut else None
        rec["origin"] = "on-target" if cut else "unattributed"
        # which genes does it hit?
        hit_genes = [g["name"] for g in design["genes"]
                     if not (g["end"] < start or g["start"] > end)]
        rec["genes_affected"] = hit_genes
        svs.append(rec)

    # ---- coverage ---------------------------------------------------------
    cov = {}
    for s in SAMPLES:
        cov[s] = binned_coverage(WORK / "cov" / ("%s.depth.txt" % s))
    mean_wt = sum(cov["WT_parent"]) / len(cov["WT_parent"])
    mean_ed = sum(cov["EDIT_pool"]) / len(cov["EDIT_pool"])
    ratio = []
    for a, b in zip(cov["EDIT_pool"], cov["WT_parent"]):
        na, nb = a / mean_ed, b / mean_wt
        ratio.append(round(na / nb, 4) if nb > 0.05 else None)

    # ---- QC summary -------------------------------------------------------
    qc = {}
    for s in SAMPLES:
        fp = json.load(open(WORK / "qc" / ("%s.fastp.json" % s)))
        flag = (WORK / "align" / ("%s.flagstat.txt" % s)).read_text()
        depths = [int(l.split("\t")[2]) for l in
                  open(WORK / "cov" / ("%s.depth.txt" % s))]
        dup = int(re.search(r"(\d+) \+ \d+ duplicates", flag).group(1))
        total = int(re.search(r"(\d+) \+ \d+ in total", flag).group(1))
        mapped = int(re.search(r"(\d+) \+ \d+ mapped \(", flag).group(1))
        qc[s] = dict(
            raw_reads=fp["summary"]["before_filtering"]["total_reads"],
            raw_bases=fp["summary"]["before_filtering"]["total_bases"],
            q30_before=round(100 * fp["summary"]["before_filtering"]["q30_rate"], 2),
            clean_reads=fp["summary"]["after_filtering"]["total_reads"],
            q30_after=round(100 * fp["summary"]["after_filtering"]["q30_rate"], 2),
            adapter_trimmed=fp.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0),
            duplication_rate=round(100 * dup / total, 2),
            mapping_rate=round(100 * mapped / total, 2),
            mean_depth=round(sum(depths) / len(depths), 2),
            pct_ge_10x=round(100 * sum(1 for d in depths if d >= 10) / len(depths), 2),
            pct_ge_30x=round(100 * sum(1 for d in depths if d >= 30) / len(depths), 2),
        )

    # ---- score against ground truth --------------------------------------
    truth = []
    with open(REF / "truth_variants.tsv") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            truth.append(dict(zip(header, line.rstrip("\n").split("\t"))))

    def matched(t):
        """Match truth to a call, allowing for VCF left-alignment shifts.

        A deletion is matched by its interval, so that a large SV is paired
        with the truth record that sits inside it rather than with whichever
        breakpoint happens to be within 30 bp.
        """
        tp = int(t["pos"])
        for v in variants:
            if abs(v["pos"] - tp) <= 30:
                return ("small-variant", v["pos"], v.get("population_freq"),
                        v.get("ci95_population"), v.get("origin"))
        for s in svs:
            if abs(s["start"] - tp) <= 30 or (s["start"] <= tp <= s["end"]):
                return ("SV", "%d-%d" % (s["start"], s["end"]),
                        s["genotypes"]["EDIT_pool"]["vaf"],
                        s["genotypes"]["EDIT_pool"].get("ci95"), s.get("origin"))
        return (None, None, None, None, None)

    scoring = []
    for t in truth:
        how, where, obs, ci, origin = matched(t)
        exp = float(t["expected_vaf_edit"])
        scoring.append(dict(truth_id=t["id"], klass=t["klass"], pos=int(t["pos"]),
                            kind=t["kind"], detail=t["detail"],
                            expected_vaf_edit=exp,
                            observed_freq=obs, ci95=ci, called_origin=origin,
                            freq_error=(round(obs - exp, 4) if obs is not None else None),
                            detected=how is not None, detected_by=how,
                            detected_at=where, note=t["note"]))
    n_truth = len(scoring)
    n_found = sum(1 for s in scoring if s["detected"])

    # false positives: EDIT-specific calls that match no truth record
    truth_pos = [int(t["pos"]) for t in truth]
    fps = [v["pos"] for v in variants
           if v["category"] == "EDIT-specific"
           and not any(abs(v["pos"] - tp) <= 30 for tp in truth_pos)]

    # ---- write outputs ----------------------------------------------------
    with open(RES / "variants.tsv", "w") as fh:
        cols = ["pos", "ref", "alt", "kind", "qual", "dp", "mq", "category",
                "origin", "attributed_to", "cut_distance", "site_mismatches",
                "vaf_wt_pileup", "vaf_edit_pileup", "vaf_edit_span", "ci95_edit",
                "informative_depth_edit", "local_cov_ratio", "population_freq",
                "ci95_population", "competing_allele",
                "effect", "impact", "gene", "hgvs_c", "hgvs_p"]
        fh.write("\t".join(cols) + "\n")
        for v in variants:
            fh.write("\t".join(str(v.get(c, "")) for c in cols) + "\n")

    with open(RES / "structural_variants.tsv", "w") as fh:
        fh.write("chrom\tstart\tend\tsvtype\tlength\tfilter\tPE\tSR\torigin\t"
                 "genes_affected\tWT_vaf\tEDIT_vaf\n")
        for s in svs:
            fh.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%s\t%.3f\t%.3f\n" % (
                s["chrom"], s["start"], s["end"], s["svtype"], s["length"],
                s["filter"], s["pe"], s["sr"], s["origin"],
                ",".join(s["genes_affected"]),
                s["genotypes"]["WT_parent"]["vaf"], s["genotypes"]["EDIT_pool"]["vaf"]))

    with open(RES / "validation.tsv", "w") as fh:
        fh.write("truth_id\tklass\tpos\tkind\tdetail\texpected_freq\tobserved_freq\t"
                 "error\tdetected\tdetected_by\tdetected_at\tcalled_origin\n")
        for s in scoring:
            fh.write("%s\t%s\t%d\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                s["truth_id"], s["klass"], s["pos"], s["kind"], s["detail"],
                s["expected_vaf_edit"],
                "-" if s["observed_freq"] is None else "%.3f" % s["observed_freq"],
                "-" if s["freq_error"] is None else "%+.3f" % s["freq_error"],
                s["detected"], s["detected_by"] or "-", s["detected_at"] or "-",
                s["called_origin"] or "-"))

    analysis = dict(
        chrom=CHROM, genome_len=len(genome), design=design["guides"],
        qc=qc, variants=variants, svs=svs,
        offtargets=ot_rows[:40],
        coverage=dict(bin_size=200, wt=cov["WT_parent"], edit=cov["EDIT_pool"],
                      ratio=ratio, mean_wt=round(mean_wt, 2), mean_edit=round(mean_ed, 2)),
        validation=dict(rows=scoring, n_truth=n_truth, n_detected=n_found,
                        false_positives=fps),
    )
    with open(RES / "analysis.json", "w") as fh:
        json.dump(analysis, fh, indent=2)

    # ---- console summary --------------------------------------------------
    print("off-target scan: %d near-cognate sites (<=5 mm, NGG/NAG PAM)" % len(ot_rows))
    for h in ot_rows[:8]:
        print("   %-8s %-7s %s:%d  %d mm (seed %d)  PAM=%s  tier=%s"
              % (h["guide"], h["strand"], CHROM, h["start"], h["mismatches"],
                 h["seed_mismatches"], h["pam"], h["tier"]))
    print()
    print("%-8s %-5s %-18s %-11s %-7s %-15s %-7s %-7s %s" %
          ("POS", "KIND", "CATEGORY", "ORIGIN", "VAF", "95% CI", "cov", "pop%", "EFFECT"))
    for v in variants:
        ci = v.get("ci95_edit") or [None, None]
        print("%-8d %-5s %-18s %-11s %-7s %-15s %-7s %-7s %s|%s" %
              (v["pos"], v["kind"], v["category"], v["origin"],
               v["vaf_edit_span"], "%s-%s" % (ci[0], ci[1]),
               v.get("local_cov_ratio"), v.get("population_freq"),
               v["effect"], v["gene"]))
    print()
    for s in svs:
        print("SV %s:%d-%d %s %d bp  origin=%s  WT VAF=%.2f  EDIT VAF=%.2f  genes=%s"
              % (s["chrom"], s["start"], s["end"], s["svtype"], s["length"],
                 s["origin"], s["genotypes"]["WT_parent"]["vaf"],
                 s["genotypes"]["EDIT_pool"]["vaf"], ",".join(s["genes_affected"])))
    print()
    print("validation: %d/%d truth events recovered, %d unexplained EDIT-specific calls"
          % (n_found, n_truth, len(fps)))
    print("%-5s %-20s %-9s %-9s %s" % ("ID", "CLASS", "DESIGN", "OBSERVED", "ERROR"))
    for s in scoring:
        print("%-5s %-20s %-9.2f %-9s %s" % (
            s["truth_id"], s["klass"], s["expected_vaf_edit"],
            "-" if s["observed_freq"] is None else "%.3f" % s["observed_freq"],
            "-" if s["freq_error"] is None else "%+.3f" % s["freq_error"]))
    for s in scoring:
        if not s["detected"]:
            print("   MISSED: %s @%d (%s)" % (s["truth_id"], s["pos"], s["note"]))


if __name__ == "__main__":
    main()
