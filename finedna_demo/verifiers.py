"""Independent re-computation of FineDNA `rl_mix` gold answers.

Each checker takes the prompt text of a `basic_dna` row and returns our own
answer, computed from the sequence with plain Python. `compare()` then checks it
against the dataset's `answer_text`. None means "could not parse this prompt",
which is counted separately from a disagreement.
"""
import re

COMP = str.maketrans("ACGTU", "TGCAA")
CODON = {}
_bases = "TCAG"
_aa3 = ("Phe Phe Leu Leu Ser Ser Ser Ser Tyr Tyr Stop Stop Cys Cys Stop Trp "
        "Leu Leu Leu Leu Pro Pro Pro Pro His His Gln Gln Arg Arg Arg Arg "
        "Ile Ile Ile Met Thr Thr Thr Thr Asn Asn Lys Lys Ser Ser Arg Arg "
        "Val Val Val Val Ala Ala Ala Ala Asp Asp Glu Glu Gly Gly Gly Gly").split()
for i, a in enumerate(_bases):
    for j, b in enumerate(_bases):
        for k, c in enumerate(_bases):
            CODON[a + b + c] = _aa3[16 * i + 4 * j + k]


def seqs(text, alphabet="ACGTU", min_len=4):
    return re.findall(rf"[{alphabet}]{{{min_len},}}", text)


def longest(text, alphabet="ACGTU", min_len=4):
    s = seqs(text, alphabet, min_len)
    return max(s, key=len) if s else None


def count_overlapping(s, sub):
    return sum(1 for i in range(len(s) - len(sub) + 1) if s[i:i + len(sub)] == sub)


def gc(s):
    return s.count("G") + s.count("C")


ENZYMES = {"EcoRI": "GAATTC", "HindIII": "AAGCTT", "BamHI": "GGATCC", "PstI": "CTGCAG",
           "XhoI": "CTCGAG", "NotI": "GCGGCCGC", "KpnI": "GGTACC", "SmaI": "CCCGGG",
           "EcoRV": "GATATC", "SalI": "GTCGAC", "XbaI": "TCTAGA", "NcoI": "CCATGG",
           "NdeI": "CATATG", "SacI": "GAGCTC", "BglII": "AGATCT", "SpeI": "ACTAGT",
           "NheI": "GCTAGC", "ClaI": "ATCGAT", "ApaI": "GGGCCC", "HpaI": "GTTAAC"}


def recognition_site(text):
    m = (re.search(r"\(([ACGT]{4,})\)", text)
         or re.search(r"(?:recognition sequence|recognizes|targets|site|sequence)[:\s]+5?'?-?([ACGT]{4,8})\b", text))
    if m:
        return m.group(1)
    names = [e for e in ENZYMES if re.search(rf"\b{e}\b", text)]
    return ENZYMES[names[0]] if len(names) == 1 else None


# --- one function per task_type -------------------------------------------

def gc_content(p):
    s = longest(p); return f"{100 * gc(s) / len(s):.1f}%"

def reverse_complement(p):
    return longest(p).translate(COMP)[::-1]

def aligned_complement(p):
    return longest(p).translate(COMP)

def transcription_coding(p):
    return longest(p).replace("T", "U")

def transcription_template(p):
    s = longest(p)
    # template written 3'->5': its complement already reads 5'->3'
    if re.search(r"3'-" + s, p):
        return s.translate(COMP).replace("T", "U")
    return s.translate(COMP)[::-1].replace("T", "U")

def pyrimidine_count(p):
    s = longest(p); return str(s.count("C") + s.count("T"))

def purine_count(p):
    s = longest(p); return str(s.count("A") + s.count("G"))

def cpg_count(p):
    return str(count_overlapping(longest(p), "CG"))

def pam_count(p):
    s = longest(p); return str(count_overlapping(s[1:], "GG"))

def hydrogen_bond_count(p):
    s = longest(p); return str(2 * (len(s) - gc(s)) + 3 * gc(s))

def wallace_tm(p):
    s = longest(p); return f"{2 * (len(s) - gc(s)) + 4 * gc(s)}°C"

def long_oligo_tm(p):
    s = longest(p); return f"{64.9 + 41 * (gc(s) - 16.4) / len(s):.1f} °C"

def complete_codons(p):
    return str(len(longest(p)) // 3)

def cpg_obs_exp(p):
    s = longest(p)
    ratio = count_overlapping(s, "CG") * len(s) / (s.count("C") * s.count("G"))
    m = re.search(r"(?:greater than|exceed|above)\s+([0-9.]+)", p)
    return "Yes" if ratio > float(m.group(1).rstrip(".")) else "No"

def restriction_site_count(p):
    site = recognition_site(p)
    return str(count_overlapping(longest(p), site)) if site else None

def restriction_fragments_linear(p):
    site = recognition_site(p)
    return str(count_overlapping(longest(p), site) + 1) if site else None

def restriction_fragments_circular(p):
    site = recognition_site(p)
    return str(max(1, count_overlapping(longest(p), site))) if site else None

def translation_coding(p):
    s = longest(p); out = []
    for i in range(0, len(s) - 2, 3):
        aa = CODON[s[i:i + 3]]
        if aa == "Stop":
            break
        out.append(aa)
    return "-".join(out)

def orf_amino_acid_count(p):
    s = longest(p); start = s.find("ATG")
    if start < 0:
        return None
    n = 0
    for i in range(start, len(s) - 2, 3):
        if CODON[s[i:i + 3]] == "Stop":
            return f"{n} amino acids"
        n += 1
    return None

def polya_signal_position(p):
    i = longest(p).find("AATAAA"); return str(i + 1) if i >= 0 else None

def mirna_seed_extraction(p):
    rna = [x for x in seqs(p, "ACGU", 15) if "U" in x]
    return rna[0][1:8] if rna else None

def shine_dalgarno_spacing(p):
    s = longest(p); i = s.find("AGGAGG")
    j = s.find("ATG", i + 6) if i >= 0 else -1
    return str(j - (i + 6)) if j >= 0 else None

def telomere_repeat_count(p):
    m = re.search(r"(?:Repeat:|motif|repeat unit|repeat)\s+5?'?-?([ACGT]{4,8})\b", p)
    if not m:
        return None
    unit, s = m.group(1), longest(p.replace(m.group(0), ""))
    best = k = 1
    while unit * k in s:
        best = k; k += 1
    return str(best if unit in s else 0)

def cas9_cut_position(p):
    s = longest(p)
    pams = [i for i in range(len(s) - 2) if s[i + 1:i + 3] == "GG"]
    if len(pams) != 1:
        return None
    n1 = pams[0] + 1  # 1-based position of the N in NGG
    return f"Between positions {n1 - 4} and {n1 - 3}"


CHECKERS = {f.__name__: f for f in [
    gc_content, reverse_complement, aligned_complement, transcription_coding,
    transcription_template, pyrimidine_count, purine_count, cpg_count, pam_count,
    hydrogen_bond_count, wallace_tm, long_oligo_tm, complete_codons, cpg_obs_exp,
    restriction_site_count, restriction_fragments_linear, restriction_fragments_circular,
    translation_coding, orf_amino_acid_count, polya_signal_position,
    mirna_seed_extraction, shine_dalgarno_spacing, telomere_repeat_count,
    cas9_cut_position,
]}


def _norm(a):
    a = a.strip().lower().replace("°c", "").replace("%", "")
    a = a.replace(" ", "").replace("aminoacids", "")
    try:
        return float(a)
    except ValueError:
        return a


def compare(ours, gold):
    a, b = _norm(ours), _norm(gold)
    if isinstance(a, float) and isinstance(b, float):
        return abs(a - b) < 0.051
    return a == b


def check(task_type, prompt, gold):
    """Return (our_answer, status) with status in {agree, disagree, unparsed}."""
    fn = CHECKERS.get(task_type)
    if fn is None:
        return None, "no_checker"
    try:
        ours = fn(prompt)
    except Exception:  # unexpected prompt shape
        ours = None
    if ours is None:
        return None, "unparsed"
    return ours, "agree" if compare(ours, gold) else "disagree"
