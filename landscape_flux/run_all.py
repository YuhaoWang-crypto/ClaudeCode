"""Run every landscape–flux demo in order and write all figures."""
import warnings

from . import (lf1_maxent, lf2_landscape_flux, lf3_flux_realdata, lf4_hamiltonian, lf8_jacobian_slowing,
               lf5_nuclear_texture, lf6_dielectric, lf7_calorimetry)


def main():
    warnings.filterwarnings("ignore")
    for m in (lf1_maxent, lf2_landscape_flux, lf3_flux_realdata, lf4_hamiltonian,
              lf5_nuclear_texture, lf6_dielectric, lf7_calorimetry, lf8_jacobian_slowing):
        m.run()
        print()


if __name__ == "__main__":
    main()
