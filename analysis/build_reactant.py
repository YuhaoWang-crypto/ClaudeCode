import numpy as np

atoms = []
for line in open("data/8WT9.pdb"):
    if line.startswith(("ATOM", "HETATM")):
        atoms.append(dict(ch=line[21], resn=line[17:20].strip(), resi=int(line[22:26]),
                          an=line[12:16].strip(),
                          xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
def g(ch, resi, an):
    for a in atoms:
        if a['ch'] == ch and a['resi'] == resi and a['an'] == an: return a['xyz'].copy()
    return None
def u(v): return v/np.linalg.norm(v)

P = g('H',18,'P'); O5 = g('H',18,"O5'"); OP1 = g('H',18,'OP1'); OP2 = g('H',18,'OP2')
O3 = g('H',17,"O3'"); C3 = g('H',17,"C3'"); C5 = g('H',18,"C5'")
Mg1 = g('B',401,'MG'); w1 = g('B',501,'O'); w2 = g('B',502,'O')
D105OD2 = g('B',105,'OD2')

# 1) re-ligate: move O3' to 1.62 A from P (intact pre-cleavage diester)
u3 = u(O3 - P)
O3_lig = P + 1.62*u3
# 2) in-line nucleophile opposite the (future) leaving O3'
Wnuc = P - 3.0*u3
# 3) second metal (canonical two-metal B site): coordinate OP2 + nucleophile + D105,
#    ~3.3 A from P and ~3.8 A from Mg1. Solve numerically for a clash-free position.
def obj(m):
    return ((np.linalg.norm(m-OP2)-2.1)**2 + (np.linalg.norm(m-Wnuc)-2.2)**2
            + 0.5*(np.linalg.norm(m-D105OD2)-2.2)**2
            + 8*max(0, 3.1-np.linalg.norm(m-P))**2
            + 8*max(0, 3.4-np.linalg.norm(m-Mg1))**2)
best, bo = None, 1e9
rng = np.random.default_rng(0)
seed = OP2 + 2.1*u(OP2 - (P+Mg1)/2)
for _ in range(4000):
    cand = seed + rng.normal(0, 1.2, 3)
    o = obj(cand)
    if o < bo: bo, best = o, cand
# local polish
for _ in range(2000):
    cand = best + rng.normal(0, 0.15, 3)
    o = obj(cand)
    if o < bo: bo, best = o, cand
Mg2 = best

def rec(el, xyz, name, resn, resi):
    return (f"HETATM{9000+resi:5d} {name:<4s}{resn:>3s} X{resi:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {el:>2s}\n")

# assemble cluster atom list (element, coord, atomname, resname, resi)
members = [
    ('Mg', Mg1, 'MG', 'MG', 401), ('Mg', Mg2, 'MG', 'MG', 402),
    ('O', w1, 'O', 'HOH', 501), ('O', w2, 'O', 'HOH', 502), ('O', Wnuc, 'O', 'HOH', 900),
    ('C', g('B',11,'CB'),'CB','ASP',11), ('C', g('B',11,'CG'),'CG','ASP',11),
    ('O', g('B',11,'OD1'),'OD1','ASP',11), ('O', g('B',11,'OD2'),'OD2','ASP',11),
    ('C', g('B',60,'CG'),'CG','GLU',60), ('C', g('B',60,'CD'),'CD','GLU',60),
    ('O', g('B',60,'OE1'),'OE1','GLU',60), ('O', g('B',60,'OE2'),'OE2','GLU',60),
    ('C', g('B',105,'CB'),'CB','ASP',105), ('C', g('B',105,'CG'),'CG','ASP',105),
    ('O', g('B',105,'OD1'),'OD1','ASP',105), ('O', g('B',105,'OD2'),'OD2','ASP',105),
    ('C', C3,'C3p','DC',17), ('O', O3_lig,'O3p','DC',17),
    ('P', P,'P','DC',18), ('O', OP1,'OP1','DC',18), ('O', OP2,'OP2','DC',18),
    ('O', O5,'O5p','DC',18), ('C', C5,'C5p','DC',18),
    ('C', g('B',13,'CA'),'CA','ALA',13), ('C', g('B',13,'CB'),'CB','ALA',13),
]
def write(fn, ser=False):
    lines = [rec(el, xyz, nm, rn, ri) for (el, xyz, nm, rn, ri) in members]
    if ser:
        CA = g('B',13,'CA'); CB = g('B',13,'CB')
        OG = CB + 1.42*u(CB - CA)
        lines.append(rec('O', OG, 'OG', 'SER', 13))
    open(fn, "w").writelines(lines + ["END\n"])
write("mdqm/reactant_wt.pdb", ser=False)
write("mdqm/reactant_ser.pdb", ser=True)

# sanity geometry
def d(a,b): return round(float(np.linalg.norm(a-b)),2)
print("=== pre-cleavage Michaelis reactant (built) ===")
print("P-O3'(re-ligated):", d(P,O3_lig), " P-O5'(ester):", d(P,O5))
print("P-Wnuc (nucleophile, in-line):", d(P,Wnuc))
print("angle Wnuc-P-O3' :", round(float(np.degrees(np.arccos(np.dot(u(Wnuc-P),u(O3_lig-P))))),1), "deg (want ~180)")
print("Mg1-P:", d(Mg1,P), " Mg2-P:", d(Mg2,P), " Mg1-Mg2:", d(Mg1,Mg2), "(want ~3.5-4.5)")
print("Mg2-OP2:", d(Mg2,OP2), " Mg2-D105OD2:", d(Mg2,D105OD2), " Mg2-Wnuc:", d(Mg2,Wnuc))
print("Mg1-w1:", d(Mg1,w1), " Mg1-w2:", d(Mg1,w2), " Mg1-OD1(D11):", d(Mg1,g('B',11,'OD1')))
print("clash check Wnuc to any heavy <1.5?:",
      [ (m[2], d(m[1],Wnuc)) for m in members if d(m[1],Wnuc)<2.0 and not np.allclose(m[1],Wnuc)])
