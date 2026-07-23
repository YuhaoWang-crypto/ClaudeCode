#!/usr/bin/env python3
"""ORM1 vs ORM2 MSA -> ORM1-divergent epitope, avoiding N-glycosylation sequons
(N-X-[S/T], X!=P) since AGP is ~45% carbohydrate (glycan-shielded surface)."""
import re
from Bio.Align import PairwiseAligner, substitution_matrices

# mature chains (signal peptide 1-18 removed; mature residue 1 = Q)
ORM1 = "QIPLCANLVPVPITNATLDRITGKWFYIASAFRNEEYNKSVQEIQATFFYFTPNKTEDTIFLREYQTRQDQCIYNTTYLNVQRENGTISRYVGGQEHFAHLLILRDTKTYMLAFDVNDEKNWGLSVYADKPETTKEQLGEFYEALDCLRIPKSDVVYTDWKKDKCEPLEKQHEKERKQEEGES"
ORM2 = "QIPLCANLVPVPITNATLDRITGKWFYIASAFRNEEYNKSVQEIQATFFYFTPNKTEDTIFLREYQTRQNQCFYNSSYLNVQRENGTVSRYEGGREHVAHLLFLRDTKTLMFGSYLDDEKNWGLSFYADKPETTKEQLGEFYEALDCLCIPRSDVMYTDWKKDKCEPLEKQHEKERKQEEGES"
print(f"mature ORM1 len={len(ORM1)}  ORM2 len={len(ORM2)}")

aln = PairwiseAligner(); aln.substitution_matrix = substitution_matrices.load("BLOSUM62")
aln.open_gap_score=-11; aln.extend_gap_score=-1; aln.mode="global"
a = aln.align(ORM1, ORM2)[0]
s1, s2 = str(a[0]), str(a[1])
idn = sum(1 for x,y in zip(s1,s2) if x==y and x!='-')
col = sum(1 for x,y in zip(s1,s2) if x!='-' and y!='-')
print(f"ORM1 vs ORM2 identity: {idn}/{col} = {100*idn/col:.0f}%\n")

# N-glyc sequons on ORM1 (mature numbering)
glyc = set()
for m in re.finditer(r"N[^P][ST]", ORM1):
    glyc.add(m.start())  # the Asn position (0-based)
glyc_sites = sorted(g+1 for g in glyc)
print(f"ORM1 N-glyc sequons (Asn, mature #): {glyc_sites}  -> glycan-shielded, AVOID")

# per-ORM1-residue divergence
m2 = {}; ci=0
for x,y in zip(s1,s2):
    if x!='-':
        m2[ci]=y; ci+=1
div = [i for i,r in enumerate(ORM1) if r!=m2.get(i,'-')]
# contiguous divergent blocks >=3
blocks=[]; st=None
for i in range(len(ORM1)):
    if i in set(div):
        if st is None: st=i
    else:
        if st is not None: blocks.append((st,i-1)); st=None
if st is not None: blocks.append((st,len(ORM1)-1))
print("\nORM1-divergent blocks (mature #):")
for aa,bb in blocks:
    if bb-aa>=2:
        seg=ORM1[aa:bb+1]
        # glyc Asn within +-2 of block?
        shielded = any(aa-2<=g<=bb+2 for g in glyc)
        tag = "  <-- near glycan (caution)" if shielded else "  <-- clean (good epitope)"
        print(f"  {aa+1:>3}-{bb+1:<3} {seg}{tag}")
