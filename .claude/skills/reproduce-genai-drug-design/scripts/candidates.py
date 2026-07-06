"""Candidate molecules to dock / analyse — the SMILES from the paper appendix.
REPLACE these with your target paper's candidates (keep the reference/hit first).
The example below is the DprE1 set (TCA1 + GTD_9.1..9.10)."""

CANDIDATES = {
    "TCA1":    "O=C(OCC)NC(C1=C(NC(C2=NC3=CC=CC=C3S2)=O)SC=C1)=O",
    "GTD_9.1": "O=C(NCC(F)F)CCC1=C(C2=CC(CCNC(C)=O)=CC(N)=C2C)C3=C(O[C@@]4([H])[C@]3([H])CCCC4)C=N1",
    "GTD_9.2": "CN1[C@@H](C[C@@H]2OCCC[C@@H]12)CNC(C(C)(c3cc(c(cn4)cc5c4occ5)c(NCC(F)F)c(N)c3)C)=O",
    "GTD_9.3": "Cc1c(C2=C(N3CCC(OCCC[C@@H]4C[C@H]5CCC[C@@H]5N4)CC3)COC2=O)nc(C(NCC(F)F)=O)cc1N",
    "GTD_9.4": "CCC(Cc1cc(c2c3c(O[C@@H]4CC[C@@H](C[C@@H]43)C)cnc2CCC(NCC(F)F)=O)c(C)c(CN)c1)=O",
    "GTD_9.5": "Cc(c(C(NCC(F)F)=O)nc1)c2c1[nH]cc2CCO[C@@H]3CCN(C4=Cc5c(CC4)c6c(CCO6)cc5N)C3",
    "GTD_9.6": "CC1=CC(N2CC[C@H](NC[C@@H]3[C@@H](C(C)(c4cc(CC(F)F)c(CO)c(N)c4)C)Cn5c3ncc5)C2)=NC1",
    "GTD_9.7": "CC(CC1CCOCC1)(c2cc(C3=NC(C(N4CC[C@H](C4)CNCc5ccccc5)=C3)=O)c(F)c(N)c2)C",
    "GTD_9.8": "Cc1nc2c(CNC2)c(C[C@@](C(NCC(F)F)=O)(c3cc(CC[C@H]4Cc5c(O4)cccc5)c(C)c(N)c3)C)n1",
    "GTD_9.9": "CN1[C@@H](C[C@H]2CCCC[C@H]12)CNC(C(C)(c3cc(c4ncc5c(CCO5)c4)c(CCC(F)F)c(N)c3)C)=O",
    "GTD_9.10": "Cc1c(C(NCC(F)F)=O)cc(C2=C(N3CCC(OCCC[C@H]4Cc5c(O4)cccc5)CC3)COC2=O)nc1N",
}
