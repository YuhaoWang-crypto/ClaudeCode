# Building the Au(111) slab + adsorbate (SERS)

## Slab
- Au fcc, a = 4.08 Å. Cleave (111).
- 4×4 surface supercell, 4-5 layers, fix bottom 2-3 layers to bulk.
- Vacuum ≥ 15 Å along z (plus dipole correction, LDIPOL/IDIPOL=3).
- k-points: e.g. 3×3×1 Monkhorst-Pack for a 4×4 cell (converge it).

```python
from ase.build import fcc111, add_adsorbate
slab = fcc111('Au', size=(4,4,5), a=4.08, vacuum=15.0)
# add_adsorbate(slab, molecule, height=2.3, position=...) for each binding mode
```

## Adsorption modes to compare
- anchor via the amino N, the carboxylate O, or (if present) thiol S;
- try top / bridge / fcc-hollow registries;
- keep the LOWEST-energy pose per enantiomer.

## The chirality result
For L and D separately: relax (INCAR.relax) -> E_ads, then frequencies
(INCAR.raman_freq). Because the free-molecule spectra are identical (demo #6),
any L/D difference in the ADSORBED spectra comes purely from the different
adsorption geometry / footprint on the surface — that is the SERS observable.
