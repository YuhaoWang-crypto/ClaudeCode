#!/usr/bin/env python3
"""
G034 / TXLNA(IL-14alpha) — finalize top-5 siRNA duplexes + U6 shRNA constructs.
Emits: siRNA duplex (guide/passenger, dTdT overhangs, 2'chemistry map, LNA option),
and U6-driven shRNA cloning oligos (sense-loop-antisense-TTTTTT) with strand annotation.
Off-target ranking is finalized separately from BLAST (sirna_offtarget_blast.txt).
"""
CANDS=[  # (id, CDS_pos, sense19)
 ("siG034_1",956,"CGAGGAGCATATCGACAAA"),
 ("siG034_2",929,"GAAGCTGATTGAGCAGTAT"),
 ("siG034_3",857,"GCAGATGGAACAGCACAAT"),
 ("siG034_4",1111,"GGATGTGTGAGCTGATGAA"),
 ("siG034_5",1158,"CTTGCCCTATACACAGAGA"),
]
def rc(s):
    c={'A':'T','T':'A','G':'C','C':'G'}; return ''.join(c[b] for b in reversed(s))
def to_rna(s): return s.replace('T','U')
LOOP="TTCAAGAGA"   # canonical 9-nt shRNA loop (Brummelkamp)
TERM="TTTTTT"      # Pol III (U6) terminator

print("="*78)
print("G034 siRNA DUPLEXES (target TXLNA/IL-14alpha; 19-mer core + 3' dTdT)")
print("="*78)
for cid,pos,sense in CANDS:
    guide=rc(sense)  # antisense = guide (loaded into RISC), targets mRNA
    gc=100*(sense.count('G')+sense.count('C'))/19
    print(f"\n[{cid}]  TXLNA CDS pos {pos}  GC={gc:.0f}%")
    print(f"  passenger(sense) 5'-{to_rna(sense)}dTdT-3'")
    print(f"  guide(antisense) 3'-dTdT{to_rna(rc(guide))}-5'  == 5'-{to_rna(guide)}dTdT-3'")
    # chemistry map: alternating 2'-OMe(o)/2'-F(f); PS(*) at 2 terminal linkages each end; guide 5'-VP
    chem=[]
    for i,b in enumerate(guide):
        chem.append(('f' if i%2 else 'o'))
    print(f"  guide 2'-chem (o=2'OMe f=2'F): {''.join(chem)}   5'-(E)-vinylphosphonate; PS x2 both ends")
    # LNA option: LNA at guide positions 2 and 12-14 lowers seed-off-target & boosts potency (Bramsen 2009)
    print(f"  LNA option: place LNA at guide positions 2, 12, 14 (1-indexed) to suppress seed off-target")

print("\n"+"="*78)
print("U6 shRNA CONSTRUCTS  (5'-sense-LOOP-antisense-TTTTTT-3'; loop=%s)"%LOOP)
print("clone into pLKO.1 / U6 cassette; annotate: sense=passenger, antisense=guide"%())
print("="*78)
for cid,pos,sense in CANDS:
    anti=rc(sense)
    hairpin=sense+LOOP+anti+TERM
    # sense-strand cloning oligo with AgeI(5')/EcoRI(3') style overhangs (pLKO.1)
    top="CCGG"+sense+"CTCGAG"+anti+"TTTTTG"   # AgeI compatible 5', with XhoI internal(alt loop) - show canonical pLKO
    print(f"\n[{cid}_shRNA] TXLNA CDS {pos}")
    print(f"  hairpin transcript (RNA sense): 5'-{to_rna(sense)}-[loop {to_rna(LOOP)}]-{to_rna(anti)}-UUUUUU-3'")
    print(f"  pLKO.1 sense oligo   5'-CCGG {sense} CTCGAG {anti} TTTTTG-3'")
    print(f"  pLKO.1 antisense oligo 5'-AATTCAAAAA {sense} CTCGAG {anti} -3'  (rev-comp, EcoRI)")
    print(f"  strand roles: guide(=antisense, loaded to RISC) = 5'-{to_rna(anti)}-3'  -> targets TXLNA mRNA @CDS{pos}")
