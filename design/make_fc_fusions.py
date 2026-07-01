#!/usr/bin/env python3
"""
Build humanized Fc-fusion ("peptibody") constructs for the lead binders.

Format: [binder peptide] - (G4S)x3 flexible linker - human IgG1 hinge-CH2-CH3.
The hinge cysteines drive homodimerization -> a BIVALENT reagent (avidity gain,
long serum half-life via FcRn). An IgG4-S228P Fc option is provided when
effector-silent neutralization is preferred.

These are ready for gene synthesis + mammalian (HEK/CHO) expression as secreted
Fc-fusion homodimers, or the bare peptides can be made by SPPS for primary assays.
"""
import json

# Human IgG1 Fc (hinge-CH2-CH3, EU ~216-447). Standard Fc-fusion cassette.
IGG1_FC = ("EPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVE"
           "VHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTL"
           "PPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQ"
           "GNVFSCSVMHEALHNHYTQKSLSLSPGK")
# IgG4 Fc with S228P (stabilized hinge), effector-attenuated
IGG4_FC = ("ESKYGPPCPPCPAPEFLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPEVQFNWYVDGVEVHN"
           "AKTKPREEQFNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKGLPSSIEKTISKAKGQPREPQVYTLPPS"
           "QEEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSRLTVDKSRWQEGNV"
           "FSCSVMHEALHNHYTQKSLSLSLGK")
LINKER = "GGGGSGGGGSGGGGS"

LEADS = {
  # id : (sequence, note)
  "TZP-P1_R2_008": ("SLSTLENELSTLENEISTIENEFSTG", "lead de novo peptide (F-anchor); best-balanced"),
  "TZP-P2_R2_002": ("SLSTLENELSTLENEISTIENEWSTG", "lead de novo peptide (W-anchor); cleanest"),
  "TZP-P3_R2_035": ("SLSTLENELSTLENEISTIENEWSTGLENEISTG", "extended, highest binding_confidence"),
  "TZP-P4_R2_010": ("SLSTLENELSTLENEISTIENEYSTG", "Tyr-anchor variant"),
}

def build(pep):
    return {
      "peptide_only": pep,
      "IgG1_Fc_fusion": pep + LINKER + IGG1_FC,
      "IgG4S228P_Fc_fusion": pep + LINKER + IGG4_FC,
    }

out={}
for name,(seq,note) in LEADS.items():
    out[name]={"note":note, **build(seq)}

json.dump(out, open("results/fc_fusion_constructs.json","w"), indent=1)

# also write a FASTA
with open("results/wetlab_constructs.fasta","w") as fh:
    for name,(seq,note) in LEADS.items():
        c=build(seq)
        fh.write(f">{name}__peptide  {note}\n{c['peptide_only']}\n")
        fh.write(f">{name}__IgG1Fc_fusion\n{c['IgG1_Fc_fusion']}\n")
        fh.write(f">{name}__IgG4S228P_Fc_fusion\n{c['IgG4S228P_Fc_fusion']}\n")
print("Wrote results/fc_fusion_constructs.json and results/wetlab_constructs.fasta")
for name,(seq,note) in LEADS.items():
    print(f"  {name}: peptide {len(seq)} aa -> IgG1 fusion {len(build(seq)['IgG1_Fc_fusion'])} aa")
