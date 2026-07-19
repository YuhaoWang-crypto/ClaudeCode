#!/usr/bin/env python3
"""Design an RNA pool biased to the cTnI cardiac N-terminal epitope (res 4-32,
Arg/Pro-rich: ...RPRPAPAPIRRRSS...). Strategy: a structured stem holds a
recognition loop; loops seeded with pyrimidine/anionic motifs that pair shape+charge
against the basic RRR/KKK patch (not pure poly-anion, to avoid promiscuity). Fold-check
with ViennaRNA, keep diverse well-folded 33-mers. Deterministic."""
import random, RNA
random.seed(424242)

STEM5 = "GGCAG"; STEM3 = "CUGCC"      # 5 bp clamp (GC-rich, stable)
# recognition-loop motif fragments aimed at Arg/Pro epitope (mix of shape + limited charge)
MOTIFS = ["UCUCC","ACUCA","UGAUC","CUCAU","AUCUC","UCAUC","CAUCU","UACUG",
          "GAUCA","UCACU","AUGUC","CUAUC","UCUGA","ACUCU","UGUCA","CAUUC"]

def make(nloops=2):
    loop = "".join(random.choice(MOTIFS) for _ in range(nloops))
    # small internal spacer to reach ~33 nt
    spacer = "".join(random.choice("ACGU") for _ in range(max(0, 23 - len(loop))))
    seq = STEM5 + loop + spacer + STEM3
    return seq[:33].ljust(33, "A")

cands, seen = [], set()
tries = 0
while len(cands) < 6 and tries < 4000:
    tries += 1
    s = make(random.choice([2,3]))
    if s in seen: continue
    st, mfe = RNA.fold(s)
    npair = st.count("(")
    # keep only genuinely structured folds (>=5 bp) with a reasonable MFE
    if npair >= 6 and mfe <= -6.0:
        seen.add(s); cands.append((f"tni{len(cands)+1}", s, st, mfe))

print("Designed cTnI-epitope RNA candidates (33 nt):")
import json
out = {}
for name, s, st, mfe in cands:
    print(f"{name}: {s}")
    print(f"      {st}  MFE={mfe:.1f}")
    out[name] = s
# a scramble decoy: shuffle tni1 (destroys structure/motif) for the gate baseline
base = cands[0][1]
sc = list(base); random.shuffle(sc); scr = "".join(sc)
print(f"scramble: {scr}  (decoy baseline)")
out["scramble"] = scr
json.dump(out, open("tni_seqs.json","w"), indent=1)
print("\n-> tni_seqs.json")
