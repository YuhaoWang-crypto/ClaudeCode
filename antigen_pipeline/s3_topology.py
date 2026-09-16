"""Step 3 -- extracellular topology gate.

An antibody, an ADC or a CAR can only engage a domain that faces the outside
of an intact cell.  This step removes proteins whose annotation places them in
the ER, Golgi, an endomembrane compartment or the cytoplasm, proteins that are
purely secreted or extracellular-matrix, and proteins whose parsed topology
leaves no extracellular domain long enough to bind.

The gate is applied to annotation, not to expression, so it is independent of
the tumour data and cannot be tuned by the ranking.
"""

from __future__ import annotations

import pandas as pd

from . import config as C


def _location_verdict(row) -> tuple[bool, str]:
    """Does any UniProt location entry put this protein on the plasma membrane?

    UniProt writes locations as a semicolon-separated list where the organelle
    comes first: 'Endoplasmic reticulum membrane (Single-pass type I membrane
    protein)' is intracellular, while a bare 'Membrane (Single-pass type I
    membrane protein)' means the plasma membrane. Matching substrings across
    the whole string confuses the two, so each entry is judged on its own.
    """
    raw = str(row.uniprot_subcellular)
    if raw.lower() in ("nan", ""):
        # no annotation: SURFY's own surface call and the parsed topology decide
        return True, ""
    entries = [e.strip().lower() for e in raw.split(";") if e.strip()]
    surface = [e for e in entries
               if any(e.startswith(p) for p in C.SURFACE_LOCATION_PREFIXES)]
    if surface:
        return True, ""
    shown = "; ".join(e.split(" (")[0] for e in entries[:3])
    return False, f"localisation: no plasma-membrane entry ({shown})"


def run(surfaceome: pd.DataFrame | None = None) -> pd.DataFrame:
    if surfaceome is None:
        surfaceome = pd.read_csv(C.RESULTS / "s1_surfaceome.csv")
    df = surfaceome.copy()
    print(f"[S3] extracellular topology gate over {len(df)} surface genes")

    verdicts = df.apply(_location_verdict, axis=1)
    df["loc_pass"] = [v[0] for v in verdicts]
    df["loc_reason"] = [v[1] for v in verdicts]

    df["ecd_pass"] = df.max_ecd_segment >= C.MIN_ECD_LENGTH
    df.loc[~df.ecd_pass, "ecd_reason"] = (
        "no extracellular segment >= " + str(C.MIN_ECD_LENGTH) + " aa"
    )
    df["ecd_reason"] = df.get("ecd_reason", pd.Series("", index=df.index)).fillna("")

    df["antibody_accessible"] = df.loc_pass & df.ecd_pass
    df["exclusion_reason"] = df.apply(
        lambda r: r.loc_reason if not r.loc_pass else (r.ecd_reason if not r.ecd_pass else ""),
        axis=1,
    )

    # Accessibility score: how much target surface an antibody actually sees.
    # Long single-pass ectodomains with many solvent-exposed glycosylation and
    # aromatic anchor sites score highest; short multi-pass loops score lowest.
    ecd = df.max_ecd_segment.clip(0, 600) / 600.0
    sites = (df.accessible_sites.fillna(0).clip(0, 20)) / 20.0
    topo_bonus = df.topology_class.map(
        lambda t: 1.0 if str(t).startswith("type I") else
                  0.9 if str(t).startswith("type II") else
                  0.85 if "GPI" in str(t) else 0.55
    ).astype(float)
    df["accessibility_score"] = (0.5 * ecd + 0.2 * sites + 0.3 * topo_bonus).clip(0, 1)
    df.loc[~df.antibody_accessible, "accessibility_score"] = 0.0

    # Provisional surface-confirmation tier; step 4 can upgrade ML-only genes
    # that Open Targets independently localises to the plasma membrane.
    df["surface_confirmation"] = df.apply(
        lambda r: "experimental (CSPA/GPI)" if r.surface_evidence == "experimental"
        else "unconfirmed (ML prediction only)", axis=1)

    passed = df[df.antibody_accessible].copy()
    excluded = df[~df.antibody_accessible].copy()
    passed.to_csv(C.RESULTS / "s3_accessible_candidates.csv", index=False)
    excluded[["gene", "protein_name", "topology_class", "max_ecd_segment",
              "uniprot_subcellular", "exclusion_reason"]].to_csv(
        C.RESULTS / "s3_excluded.csv", index=False)

    by_reason = (excluded.exclusion_reason.str.split(":").str[0]
                 .value_counts().to_dict())
    detail = excluded.exclusion_reason.value_counts().head(12).to_dict()
    print(f"  excluded {len(excluded)} genes; {len(passed)} antibody-accessible remain")
    for k, v in by_reason.items():
        print(f"    {k}: {v}")

    C.write_provenance("s3", {
        "input_genes": int(len(df)),
        "excluded": int(len(excluded)),
        "accessible": int(len(passed)),
        "exclusions_by_class": {str(k): int(v) for k, v in by_reason.items()},
        "exclusions_detail": {str(k): int(v) for k, v in detail.items()},
        "min_ecd_length": C.MIN_ECD_LENGTH,
        "negative_control_status": {
            g: ("not in surfaceome search space" if not (df.gene == g).any()
                else ("excluded: " + str(df.loc[df.gene == g, "exclusion_reason"].iloc[0])
                      if not bool(df.loc[df.gene == g, "antibody_accessible"].iloc[0])
                      else "passed the topology gate"))
            for g in C.NEGATIVE_CONTROLS
        },
    })
    return passed


if __name__ == "__main__":
    run()
