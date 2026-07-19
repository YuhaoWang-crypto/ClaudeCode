#!/usr/bin/env python3
"""Locate cardiac-specific (divergent, surface) epitope of cTnI vs skeletal TnI isoforms."""
from Bio.Align import PairwiseAligner, substitution_matrices

cTnI  = "MADGSSDAAREPRPAPAPIRRRSSNYRAYATEPHAKKKSKISASRKLQLKTLLLQIAKQELEREAEERRGEKGRALSTRCQPLELAGLGFAELQDLCRQLHARVDKVDEERYDIEAKVTKNITEIADLTQKIFDLRGKFKRPTLRRVRISADAMMQALLGARAKESLDLRAHLKQVKKEDTEKENREVGDWRKNIDALSGMEGRKKKFES"
TNNI2 = "MGDEEKRNRAITARRQHLKSVMLQIAATELEKEESRREAEKQNYLAEHCPPLHIPGSMSEVQELCKQLHAKIDAAEEEKYDMEVRVQKTSKELEDMNQKLFDLRGKFKRPPLRRVRMSADAMLKALLGSKHKVCMDLRANLKQVKKEDTEKERDLRDVGDWRKNIEEKSGMEGRKKMFESES"
TNNI1 = "MPEVERKPKITASRKLLLKSLMLAKAKECWEQEHEEREAEKVRYLAERIPTLQTRGLSLSALQDLCRELHAKVEVVDEERYDIEAKCLHNTREIKDLKLKVMDLRGKFKRPPLRRVRVSADAMLRALLGSKHKVSMDLRANLKSVKKEDTEKERPVEVGDWRKNVEAMSGMEGRKKMFDAAKSPTSQ"

aln = PairwiseAligner()
aln.substitution_matrix = substitution_matrices.load("BLOSUM62")
aln.open_gap_score = -11; aln.extend_gap_score = -1; aln.mode = "global"

def report(name, other):
    a = aln.align(cTnI, other)[0]
    print(f"\n### cTnI vs {name}   score={a.score:.0f}")
    # identity over aligned (non-gap) columns
    s1, s2 = str(a[0]), str(a[1])
    idn = sum(1 for x,y in zip(s1,s2) if x==y and x!='-')
    col = sum(1 for x,y in zip(s1,s2) if x!='-' and y!='-')
    print(f"identity: {idn}/{col} aligned = {100*idn/col:.0f}%")
    return s1, s2

s1_2, s2_2 = report("TNNI2 (fast skel)", TNNI2)
s1_1, s2_1 = report("TNNI1 (slow skel)", TNNI1)

# Map: how many N-terminal cTnI residues have NO skeletal counterpart (gap in BOTH)?
# Use the TNNI1 alignment (slow skeletal, closest to cardiac) to find the conserved anchor.
print("\n### cTnI vs TNNI1 head (first 60 aln columns): cardiac-unique N-term extension")
print("cTnI :", s1_1[:70])
print("TNNI1:", s2_1[:70])

# Per-position cardiac-specificity on the ORIGINAL cTnI numbering:
# a cTnI residue is 'cardiac-divergent' if it differs from BOTH skeletal isoforms at its aligned column.
def col_map(s1, s2):
    """map cTnI residue index (0-based) -> aligned partner char ('-' if gap)."""
    m = {}; ci = 0
    for x, y in zip(s1, s2):
        if x != '-':
            m[ci] = y; ci += 1
    return m
m2 = col_map(s1_2, s2_2)   # cTnI->TNNI2
m1 = col_map(s1_1, s2_1)   # cTnI->TNNI1
div = []
for i, r in enumerate(cTnI):
    p2 = m2.get(i, '-'); p1 = m1.get(i, '-')
    if r != p2 and r != p1:          # differs from both skeletal (or no counterpart)
        div.append(i)
# contiguous blocks
blocks = []; start = None
for i in range(len(cTnI)):
    if i in set(div):
        if start is None: start = i
    else:
        if start is not None: blocks.append((start, i-1)); start = None
if start is not None: blocks.append((start, len(cTnI)-1))
print("\n### Cardiac-divergent blocks (1-based cTnI residue numbering), differs from BOTH skeletal:")
for a,b in blocks:
    if b-a >= 3:  # blocks >=4 residues
        print(f"  {a+1:>3}-{b+1:<3}  {cTnI[a:b+1]}")
