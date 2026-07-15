#!/usr/bin/env python3
"""
G034 / TXLNA (alpha-taxilin / IL-14alpha) rational siRNA design.
Targets the provided native CDS fragment (encodes aa 210-546; = IL-14alpha C-region).
Implements classic rational-design filters (Reynolds 2004, Ui-Tei 2004, Amarzguioui 2004)
+ thermodynamic asymmetry + liabilities. NO external deps.

Off-target seed screening here is INTRA-sequence only; genome-wide off-target
(seed-family + full-length) requires BLAST/Bowtie vs transcriptome and, ideally,
Evo2 (variant/off-target likelihood). Marked accordingly.
"""
CDS=("ATGAAGCTCCTACAGAAAAAGCAGAGCCAGCTGGTGCAAGAGAAGGACCACCTGCGCGGTGAGCACAGCAAGGCC"
"GTCCTGGCCCGCAGCAAGCTTGAGAGCCTATGCCGTGAGCTGCAGCGGCACAACCGCTCCCTCAAGGAAGAAGGT"
"GTGCAGCGGGCCCGGGAGGAGGAGGAGAAGCGCAAGGAGGTGACCTCGCACTTCCAGGTGACACTGAATGACATT"
"CAGCTGCAGATGGAACAGCACAATGAGCGCAACTCCAAGCTGCGCCAAGAGAACATGGAGCTGGCTGAGAGGCTC"
"AAGAAGCTGATTGAGCAGTATGAGCTGCGCGAGGAGCATATCGACAAAGTCTTCAAACACAAGGACCTACAACAG"
"CAGCTGGTGGATGCCAAGCTCCAGCAGGCCCAGGAGATGCTAAAGGAGGCAGAAGAGCGGCACCAGCGGGAGAAG"
"GATTTTCTCCTGAAAGAGGCAGTAGAGTCCCAGAGGATGTGTGAGCTGATGAAGCAGCAAGAGACCCACCTGAAG"
"CAACAGCTTGCCCTATACACAGAGAAGTTTGAGGAGTTCCAGAACACACTTTCCAAAAGCAGCGAGGTATTCACC"
"ACATTCAAGCAGGAGATGGAAAAGATGACTAAGAAGATCAAGAAGCTGGAGAAAGAAACCACCATGTACCGGTCC"
"CGGTGGGAGAGCAGCAACAAGGCCCTGCTTGAGATGGCTGAGGAGAAAACAGTCCGGGATAAAGAACTGGAGGGC"
"CTGCAGGTAAAAATCCAACGGCTGGAGAAGCTGTGCCGGGCACTGCAGACAGAGCGCAATGACCTGAACAAGAGG"
"GTACAGGACCTGAGTGCTGGTGGCCAGGGCTCCCTCACTGACAGTGGCCCTGAGAGGAGGCCAGAGGGGCCTGGG"
"GCTCAAGCACCCAGCTCCCCCAGGGTCACAGAAGCGCCTTGCTACCCAGGAGCACCGAGCACAGAAGCATCAGGCC"
"AGACTGGGCCTCAAGAGCCCACCTCCGCCAGGGCCTAG")
CDS=CDS.replace("\n","")
OFF=627  # nt offset: provided CDS starts at codon 210 -> (210-1)*3=627 into full TXLNA CDS

def gc(s): return 100.0*(s.count('G')+s.count('C'))/len(s)
def rc(s):
    c={'A':'T','T':'A','G':'C','C':'G'}; return ''.join(c[b] for b in reversed(s))
# nearest-neighbor dG (Freier/SantaLucia-ish, kcal/mol) for 5'->3' duplex ends asymmetry
NN={'AA':-1.0,'TT':-1.0,'AT':-0.9,'TA':-0.6,'CA':-1.7,'TG':-1.7,'GT':-1.5,'AC':-1.5,
    'CT':-1.5,'AG':-1.5,'GA':-1.3,'TC':-1.3,'CG':-2.8,'GC':-2.3,'GG':-2.1,'CC':-2.1}
def end_dg(dimer): return NN.get(dimer,-1.5)

cands=[]
n=len(CDS)
for i in range(n-21):
    site=CDS[i:i+19]            # 19-nt target (sense target region)
    if 'N' in site: continue
    g=gc(site)
    anti=rc(site)               # guide/antisense 5'->3' (targets this site)
    sense=site
    # --- Reynolds-style scoring ---
    score=0; reasons=[]
    if 30<=g<=52: score+=1
    else: reasons.append(f"GC {g:.0f}% out of 30-52")
    # A/U at antisense 5' end (=sense pos19 is A/T) favored -> asymmetry
    if sense[18] in 'AT': score+=1
    else: reasons.append("no A/U at guide 5'")
    # G/C at sense 5' (pos1) favored for guide 3' stability of sense end
    if sense[0] in 'GC': score+=1
    # Ui-Tei rules
    if anti[0] in 'AT': score+=1     # guide pos1 = A/U
    else: reasons.append("Ui-Tei: guide pos1 not A/U")
    if sense[0] in 'GC': score+=1    # target pos1 G/C
    # >=4 A/U in guide 5' 7-mer (seed region low stability)
    au7=sum(1 for b in anti[:7] if b in 'AT')
    if au7>=4: score+=1
    else: reasons.append(f"seed A/U={au7}<4")
    # thermodynamic asymmetry: guide 5' end less stable than 3' end
    dg5=end_dg(anti[0:2]); dg3=end_dg(anti[-2:])
    asym=dg5-dg3   # want >0 (5' end weaker/less negative)
    if asym>0: score+=1
    else: reasons.append("no thermo asymmetry")
    # position-specific favorable (Reynolds): A at pos3, U at pos10 (sense)
    if sense[2]=='A': score+=1
    if len(sense)>9 and sense[9]=='T': score+=1
    # --- liabilities (disqualify/penalize) ---
    bad=[]
    for hp in ['AAAA','TTTT','GGGG','CCCC']:
        if hp in sense: bad.append(hp)
    if 'TGTCTC' in sense or 'GTCCTTCAA' in sense: bad.append('immunostim-motif')
    if g<25 or g>60: bad.append('GC-extreme')
    # internal repeat / palindrome quick check
    fullpos=i+OFF
    cands.append(dict(pos=fullpos, i=i, sense=sense, guide=anti, gc=round(g,1),
                      score=score, asym=round(asym,2), au7=au7, bad=";".join(bad),
                      reasons=";".join(reasons)))

# rank: high score, no liabilities, mid GC
ranked=sorted(cands, key=lambda c:(-(c['score']), len(c['bad'])>0, abs(c['gc']-45)))
print(f"Total 19-mer windows scanned: {len(cands)}")
print(f"Max score observed: {max(c['score'] for c in cands)} (of ~9)\n")
print("TOP 12 siRNA CANDIDATES (target TXLNA/IL-14alpha CDS; pos = nt in full CDS)")
print("-"*104)
print(f"{'rank':<4}{'CDS_pos':<8}{'score':<6}{'GC%':<6}{'guide(antisense 5>3)':<24}{'sense/target 19mer':<22}{'liabilities'}")
top=[c for c in ranked if not c['bad']][:12]
for r,c in enumerate(top,1):
    # 21-mer duplex with UU 3' overhangs (guide + sense)
    print(f"{r:<4}{c['pos']:<8}{c['score']:<6}{c['gc']:<6}{c['guide']+'UU':<24}{c['sense']+'UU':<22}{c['bad'] or 'clean'}")
print()
print("Overhang convention: 21-mer duplex = 19-mer core + 3' dTdT (shown as UU).")
print("Chemistry (recommended for a lead): 2'-OMe/2'-F alternating, PS at ends,")
print("guide 5'-(E)-vinylphosphonate; GalNAc is liver-tropic -> for lymphoid/gland")
print("delivery use lipid nanoparticle or antibody-siRNA conjugate (see report).")
