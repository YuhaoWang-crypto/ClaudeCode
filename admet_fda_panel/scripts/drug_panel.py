"""FDA-approved drug panel for the ADMET / safety-panel demonstration.

30 marketed (or formerly marketed) small-molecule drugs chosen to span
therapeutic classes, physicochemical space, and — importantly for the
secondary-pharmacology mapping exercise — to include known positive controls
for several safety-panel targets:

  hERG / QT            : terfenadine, cisapride, amiodarone, haloperidol,
                         quinidine, sertindole-like antipsychotics (risperidone)
  CYP3A4 inhibition    : ketoconazole, ritonavir-like azoles
  P-glycoprotein       : verapamil, quinidine, digoxin (substrate)
  Estrogen receptor    : tamoxifen (SERM)
  Androgen receptor    : bicalutamide (antagonist)
  Aromatase (CYP19A1)  : anastrozole (inhibitor)
  PPAR-gamma           : rosiglitazone (agonist)

SMILES are canonical, drawn from public sources (DrugBank / PubChem / ChEMBL
conventions) and validated against reference molecular formulae in
`REFERENCE_FORMULA` by `load_example_drugs()`.
"""

from __future__ import annotations

# name, SMILES, ATC-ish class, expected molecular formula (for validation)
_DRUGS = [
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O", "NSAID / antiplatelet", "C9H8O4"),
    ("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "NSAID", "C13H18O2"),
    ("Acetaminophen", "CC(=O)Nc1ccc(O)cc1", "Analgesic / antipyretic", "C8H9NO2"),
    ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "CNS stimulant", "C8H10N4O2"),
    ("Warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O", "Anticoagulant", "C19H16O4"),
    ("Atorvastatin",
     "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
     "Statin / lipid-lowering", "C33H35FN2O5"),
    ("Metformin", "CN(C)C(=N)NC(=N)N", "Antidiabetic (biguanide)", "C4H11N5"),
    ("Omeprazole", "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1",
     "Proton-pump inhibitor", "C17H19N3O3S"),
    ("Simvastatin",
     "CCC(C)(C)C(=O)O[C@H]1C[C@H](C)C=C2C=C[C@H](C)[C@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@@H]12",
     "Statin / lipid-lowering", "C25H38O5"),
    ("Diazepam", "CN1c2ccc(Cl)cc2C(=NCC1=O)c1ccccc1", "Benzodiazepine", "C16H13ClN2O"),
    ("Fluoxetine", "CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1", "SSRI antidepressant",
     "C17H18F3NO"),
    ("Sertraline", "CN[C@H]1CC[C@@H](c2ccc(Cl)c(Cl)c2)c2ccccc21",
     "SSRI antidepressant", "C17H17Cl2N"),
    ("Haloperidol", "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1",
     "Typical antipsychotic", "C21H23ClFNO2"),
    ("Risperidone",
     "Cc1nc2CCCCn2c(=O)c1CCN1CCC(CC1)c1noc2cc(F)ccc12",
     "Atypical antipsychotic", "C23H27FN4O2"),
    ("Amiodarone", "CCCCc1oc2ccccc2c1C(=O)c1cc(I)c(OCCN(CC)CC)c(I)c1",
     "Class III antiarrhythmic", "C25H29I2NO3"),
    ("Terfenadine",
     "CC(C)(C)c1ccc(cc1)C(O)CCCN1CCC(CC1)C(O)(c1ccccc1)c1ccccc1",
     "Antihistamine (withdrawn)", "C32H41NO2"),
    ("Cisapride",
     "COc1cc(C(=O)NC2CCN(CCCOc3ccc(F)cc3)CC2OC)c(N)cc1Cl",
     "Prokinetic (withdrawn)", "C23H29ClFN3O4"),
    ("Verapamil",
     "COc1ccc(CCN(C)CCCC(C#N)(C(C)C)c2ccc(OC)c(OC)c2)cc1OC",
     "Calcium-channel blocker", "C27H38N2O4"),
    ("Propranolol", "CC(C)NCC(O)COc1cccc2ccccc12", "Beta-blocker", "C16H21NO2"),
    ("Ketoconazole",
     "CC(=O)N1CCN(CC1)c1ccc(OC[C@@H]2CO[C@](Cn3ccnc3)(O2)c2ccc(Cl)cc2Cl)cc1",
     "Azole antifungal", "C26H28Cl2N4O4"),
    ("Ciprofloxacin",
     "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
     "Fluoroquinolone antibiotic", "C17H18FN3O3"),
    ("Imatinib",
     "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
     "BCR-ABL kinase inhibitor", "C29H31N7O"),
    ("Gefitinib",
     "COc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OCCCN1CCOCC1",
     "EGFR kinase inhibitor", "C22H24ClFN4O3"),
    ("Sildenafil",
     "CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1",
     "PDE5 inhibitor", "C22H30N6O4S"),
    ("Tamoxifen",
     "CC/C(=C(\\c1ccccc1)/c1ccc(OCCN(C)C)cc1)/c1ccccc1",
     "SERM / antineoplastic", "C26H29NO"),
    ("Rosiglitazone",
     "CN(CCOc1ccc(CC2SC(=O)NC2=O)cc1)c1ccccn1",
     "PPAR-gamma agonist (TZD)", "C18H19N3O3S"),
    ("Anastrozole",
     "CC(C)(C#N)c1cc(Cn2cncn2)cc(C(C)(C)C#N)c1",
     "Aromatase inhibitor", "C17H19N5"),
    ("Bicalutamide",
     "CC(O)(CS(=O)(=O)c1ccc(F)cc1)C(=O)Nc1ccc(C#N)c(c1)C(F)(F)F",
     "Androgen-receptor antagonist", "C18H14F4N2O4S"),
    ("Quinidine",
     "COc1ccc2nccc([C@@H](O)[C@H]3C[C@@H]4CC[N@]3C[C@@H]4C=C)c2c1",
     "Class Ia antiarrhythmic", "C20H24N2O2"),
    ("Loratadine",
     "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc32)CC1",
     "Second-generation antihistamine", "C22H23ClN2O2"),
]


def load_example_drugs(validate: bool = True):
    """Return a DataFrame with columns name / smiles / drug_class.

    When `validate` is True the parsed molecular formula of every SMILES is
    checked against the curated reference formula and a ValueError is raised
    on any mismatch, so a typo in a SMILES string cannot silently propagate
    into the ADMET results.
    """
    import pandas as pd
    from rdkit import Chem, RDLogger
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula

    RDLogger.DisableLog("rdApp.*")

    rows, problems = [], []
    for name, smi, cls, formula in _DRUGS:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            problems.append(f"{name}: SMILES failed to parse")
            continue
        got = CalcMolFormula(mol)
        if validate and got != formula:
            problems.append(f"{name}: formula {got} != expected {formula}")
        rows.append({"name": name, "smiles": Chem.MolToSmiles(mol),
                     "drug_class": cls, "formula": got})

    if problems:
        raise ValueError("Drug panel validation failed:\n  " + "\n  ".join(problems))

    df = pd.DataFrame(rows)
    print(f"✓ Data loaded successfully! {len(df)} molecules")
    return df


if __name__ == "__main__":
    print(load_example_drugs().to_string())
