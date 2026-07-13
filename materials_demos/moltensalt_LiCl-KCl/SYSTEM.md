# LiCl-KCl molten-salt initial structure

- composition (atom counts): Li:28  K:20  Cl:48
- total atoms: 96
- target density: 1.6 g/cm^3 (ref T = 723 K)
- cubic cell edge: 14.060 Angstrom
- files: POSCAR (VASP), init.xyz (CP2K), cell edge above

Next: run NVT AIMD at the reference T (see hpc_templates/), sample uncorrelated frames, label E/forces, package with dpdata into deepmd/npy (see moltensalt_mlp_pipeline.py pack_deepmd()).
