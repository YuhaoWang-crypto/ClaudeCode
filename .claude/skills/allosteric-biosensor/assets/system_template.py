"""
system_template.py -- copy-paste skeleton for a new biosensor system.

Fill in the receptor sequence + analyte SMILES + structure-derived loop sites,
drop the Receptor/System into biosensor_pipeline/systems.py, and the rest of the
pipeline (library build, Boltz payloads, analysis) works unchanged.

Steps:
  1. seq   -> RCSB FASTA https://www.rcsb.org/fasta/entry/<PDB> (strip His tags)
  2. SMILES-> PubChem  .../compound/name/<name>/property/CanonicalSMILES/JSON
  3. loops -> from biosensor_pipeline.structure import annotate
              _, offset, _, loop_centers = annotate("<PDB>", "A", full_seq)
"""

from biosensor_pipeline.systems import Receptor, System, TEM1

MY_RECEPTOR = Receptor(
    name="MYBINDER",
    pdb="XXXX",
    seq="...FILL_IN...",                 # tag-stripped receptor sequence
    ligand_name="my_analyte",
    ligand_smiles="C...",               # analyte SMILES (with stereo if known)
    ligand_formula="C..H..",
    loop_sites=[],                       # <10 structure-derived loop centers
    modeled_offset=0,                    # offset of modeled struct within seq
    note="what this binder is / where it came from.",
)

MY_SYSTEM = System(
    key="mytag",
    role="validation",                   # or "reproduction"
    receptor=MY_RECEPTOR,
    reporter=TEM1,                       # or another Reporter you defined
    primary_insertion="253",
    description="one-line description of the detection target.",
)

# Then in systems.py:  RECEPTORS[MY_RECEPTOR.name] = MY_RECEPTOR
#                      SYSTEMS[MY_SYSTEM.key]      = MY_SYSTEM
