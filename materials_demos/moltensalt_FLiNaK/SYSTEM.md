# FLiNaK molten-salt initial structure

- composition (atom counts): Li:23  Na:6  K:21  F:50
- total atoms: 100
- target density: 2.02 g/cm^3 (ref T = 773 K)
- cubic cell edge: 11.936 Angstrom
- files: POSCAR (VASP), init.xyz (CP2K), cell edge above

Next: run NVT AIMD at the reference T (see hpc_templates/), sample uncorrelated frames, label E/forces, package with dpdata into deepmd/npy (see moltensalt_mlp_pipeline.py pack_deepmd()).
