"""Run the whole conformational-flux pipeline."""
from . import (figures, m1_dataset, m2_cryptic_pocket, m3_msm_diagnostics,
               m4_tpt_flux, m5_variant_allostery, m6_traditional_contrast,
               m7_druggability, m8_positive_control,
               m9_adaptive_benchmark, m10_tem_projection)

MODULES = [m1_dataset, m2_cryptic_pocket, m3_msm_diagnostics, m4_tpt_flux,
           m5_variant_allostery, m6_traditional_contrast, m7_druggability,
           m8_positive_control, m9_adaptive_benchmark, m10_tem_projection]


def main():
    for m in MODULES:
        m.main()
        print()
    print("=" * 78)
    print("figures")
    print("=" * 78)
    figures.main()


if __name__ == "__main__":
    main()
