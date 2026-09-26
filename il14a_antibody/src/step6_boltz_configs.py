"""Step 6 - Emit ready-to-run, epitope-CONSTRAINED Boltz-2 configs.

WHY THIS EXISTS
---------------
The earlier campaign ran 24 Boltz-2 complexes plus an OpenFold3 cross-check and
could not obtain a reproducible antigen interface: for the TMEL1011 parent the
target-contact Jaccard between two replicate samples was 0.012, and the
residue-pair Jaccard was 0. That is a failure to converge, and the supplementary
report correctly refused to rank candidates on it.

The single most likely cause is visible in that report's own methods: "抗体重/轻
链采用各自单序列输入，没有施加表位或接触约束" - no epitope or contact constraints
were applied. An unconstrained co-fold against a 546-residue, largely
low-confidence, coiled-coil antigen has an enormous search space and will scatter
the binder over the surface.

The Boltz protein-design and protein-screen APIs both accept `epitope_residues`,
`non_binding_residues` and explicit `pocket` constraints. These configs use them:
the binder is steered onto the chosen epitope and steered AWAY from the
paralog-conserved surface. That is a materially different experiment from the one
that failed, and it is the one worth paying for.

COST (estimated from the live API at the time of writing)
--------------------------------------------------------
    protein design  $0.050 per design
    protein screen  $0.025 per candidate
So the staged plan below is roughly:
    stage 1  screen 12 existing candidates ............  ~$0.30
    stage 2  design 200 nanobodies, constrained .......  ~$10.00
    stage 3  screen top 50 re-ranked ..................  ~$1.25
    stage 4  structure+binding on 10 finalists ........  (per-prediction pricing)

Nothing here is submitted. Each file is a payload for the corresponding
Boltz MCP tool.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS

AF_CIF = "https://alphafold.ebi.ac.uk/files/AF-P40222-F1-model_v6.cif"

# Boltz uses 0-INDEXED residue indices; TXLNA numbering is 1-based.
def idx(first, last):
    return list(range(first - 1, last))


# Two epitopes, both carried forward deliberately (see REPORT.md):
#   CC  aa367-382  ordered coiled-coil, highest computed paralog-selectivity
#                  headroom, but no functional blocking evidence
#   CT  aa493-523  the 1C6/1F2 zone, the only region with reported functional
#                  blocking data, but disordered and inside US7622574
EPITOPES = {
    "CC_367_382": {
        "epitope": idx(367, 382),
        "crop": idx(300, 460),
        # paralog-conserved exposed positions on the same helix: steer away so
        # the paratope is pushed onto the TXLNA-unique face
        "non_binding": [p - 1 for p in (355, 357, 358, 359, 361, 362, 363, 364,
                                        377, 379, 380, 383, 384, 387, 389, 390)],
    },
    "CT_493_523": {
        "epitope": idx(493, 523),
        "crop": idx(430, 546),
        "non_binding": [p - 1 for p in (433, 435, 440, 445, 450, 455, 460, 465)],
    },
}


def target_block(cfg):
    return {
        "type": "structure_template",
        "structure": {"type": "url", "url": AF_CIF},
        "chain_selection": {
            "A": {
                "chain_type": "polymer",
                "crop_residues": cfg["crop"],
                "epitope_residues": cfg["epitope"],
                "non_binding_residues": [p for p in cfg["non_binding"]
                                         if p in set(cfg["crop"])
                                         and p not in set(cfg["epitope"])],
            }
        },
    }


DESIGN_RULES = {
    # no new cysteine in a CDR; no PTM/glycosylation motifs; cap hydrophobicity
    "excluded_amino_acids": ["C"],
    "excluded_sequence_motifs": ["NG", "NS", "DG", "DP", "NXS", "NXT"],
    "max_hydrophobic_fraction": 0.45,
}


def main():
    out_dir = RESULTS / "boltz_configs"
    out_dir.mkdir(parents=True, exist_ok=True)

    constructs = json.loads((RESULTS / "step5_constructs.json").read_text())
    hum_vhh = constructs["variable_domains"]["humanised_VHH"]
    hum_vh = constructs["variable_domains"]["humanised_TMEL1011_VH"]
    hum_vl = constructs["variable_domains"]["humanised_TMEL1011_VL_lambda"]

    parents = json.loads((RESULTS / "step2_parents.json").read_text())
    vhh_pre = ("QVQLVESGGGLVQAGGSLRLSCAASGFTFSSYAMGWFRQAPGKEREFVSAISGSGGSTYYA"
               "DSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYC")
    vhh_post = "WGQGTQVTVSS"

    written = []
    for tag, cfg in EPITOPES.items():
        tgt = target_block(cfg)

        # --- stage 1: screen what already exists, constrained to the epitope --
        proteins = [
            {"id": "humanised_VHH_C05",
             "entities": [{"type": "protein", "chain_ids": ["H"], "value": hum_vhh}]},
            {"id": "TMEL1011_humanised_Fv",
             "entities": [{"type": "protein", "chain_ids": ["H"], "value": hum_vh},
                          {"type": "protein", "chain_ids": ["L"], "value": hum_vl}]},
        ]
        for r in parents["parent_triage"]:
            if not r["usable_as_parent"]:
                continue
            proteins.append({
                "id": f"parent_{r['id']}",
                "entities": [{"type": "protein", "chain_ids": ["H"],
                              "value": vhh_pre + r["cdr3"] + vhh_post}]})
        p = out_dir / f"stage1_screen_{tag}.json"
        p.write_text(json.dumps({
            "_tool": "mcp__Boltz_API__boltz_start_protein_screen",
            "_note": "Run boltz_estimate_protein_screen with the same body first.",
            "target": tgt, "proteins": proteins,
            "idempotency_key": f"il14a-{tag}-stage1-screen-v1",
        }, indent=2))
        written.append(p)

        # --- stage 2: constrained de novo nanobody design --------------------
        p = out_dir / f"stage2_design_nanobody_{tag}.json"
        p.write_text(json.dumps({
            "_tool": "mcp__Boltz_API__boltz_start_protein_design",
            "_note": "Estimated ~$10.00 at $0.05/design for num_proteins=200.",
            "num_proteins": 200,
            "target": tgt,
            "binder_specification": {"type": "boltz_curated",
                                     "binder": "boltz_nanobody",
                                     "rules": DESIGN_RULES},
            "idempotency_key": f"il14a-{tag}-stage2-nb200-v1",
        }, indent=2))
        written.append(p)

        # --- stage 3: CDR3-only redesign of the humanised parent -------------
        # Keeps the humanised framework and both other CDRs fixed, redesigning
        # only CDR3 within the length window measured in step 2. This is the
        # "optimise the existing antibody" path rather than de novo.
        cdr3_start = len(vhh_pre)                       # 0-indexed
        cdr3_end = len(hum_vhh) - len(vhh_post) - 1
        p = out_dir / f"stage3_cdr3_redesign_{tag}.json"
        p.write_text(json.dumps({
            "_tool": "mcp__Boltz_API__boltz_start_protein_design",
            "_note": ("Parent-anchored: only CDR3 is redesigned; the humanised "
                      "framework and the VHH hallmark tetrad stay fixed."),
            "num_proteins": 150,
            "target": tgt,
            "binder_specification": {
                "type": "no_template",
                "modality": "nanobody",
                "entities": [{
                    "type": "designed_protein", "chain_ids": ["H"],
                    # fixed framework + CDR1/CDR2, designed CDR3 of length 10-18
                    "value": f"{hum_vhh[:cdr3_start]}10..18{vhh_post}",
                }],
                "rules": DESIGN_RULES,
            },
            "idempotency_key": f"il14a-{tag}-stage3-cdr3-v1",
        }, indent=2))
        written.append(p)

    readme = """# Boltz-2 configs (not submitted)

Each file is a request body for the named Boltz MCP tool. Always run the
matching `boltz_estimate_*` tool with the same body first and check the cost.

Order of execution:
  1. stage1_screen_*        - does anything we already have touch the epitope?
  2. stage2_design_*        - constrained de novo nanobody generation
  3. stage3_cdr3_redesign_* - parent-anchored CDR3 redesign (preferred path)

The critical difference from the earlier campaign: every target block sets
`epitope_residues` AND `non_binding_residues`. The earlier runs applied no
epitope or contact constraints and their replicate contact sets did not
reproduce (Jaccard 0.012). Treat any run whose replicate contact sets disagree
as a failed run, not as a low-affinity result.

Acceptance gate before believing any ranking:
  * two independent samples of the SAME input must agree on the contact set
    (target-residue Jaccard >= 0.6), and
  * the contacted residues must lie within the intended epitope.
Only then are the relative scores worth reading.
"""
    (out_dir / "README.md").write_text(readme)

    print("=" * 84)
    print("EPITOPE-CONSTRAINED BOLTZ CONFIGS WRITTEN (nothing submitted)")
    print("=" * 84)
    for p in written:
        d = json.loads(p.read_text())
        print(f"  {p.name:40} -> {d['_tool'].split('__')[-1]}")
    print(f"  {'README.md':40} -> acceptance gates and run order")
    print(f"\nall under: {out_dir}")
    print("\nkey change vs the campaign that failed: epitope_residues and")
    print("non_binding_residues are now set on every target block.")


if __name__ == "__main__":
    main()
