#!/usr/bin/env python3
"""
Rank the ab2 scFv CDR affinity-maturation screen against the fully-modified
Tirzepatide (peptide + Aib2/Aib13 + co-folded K20 acyl chain).

Screen job: prot_scr_SAlb4E3KrEdF3qUOg5bz  (52 scFv variants, epitope-directed)

Primary signal = binding_confidence (Boltz-2 binder head; the honest
discriminator under epitope forcing). ipTM / protein_ipTM report interface
geometry; structure_confidence reports fold plausibility of the scFv+complex.

We rank by binding_confidence and report the DELTA vs the ab2 wild-type
(M00_ab2_WT) so improvements from each CDR point mutation are explicit.
"""
import json, sys

WT_ID = "M00_ab2_WT"

def load_screen(path):
    d = json.load(open(path))
    out = {}
    for r in d["results"]:
        m = r.get("metrics", {}) or {}
        art = r.get("artifacts") or {}
        struct = (art.get("structure") or {}) if isinstance(art, dict) else {}
        out[r["external_id"]] = {
            "iptm": m.get("iptm"),
            "protein_iptm": m.get("protein_iptm"),
            "bind_conf": m.get("binding_confidence"),
            "struct_conf": m.get("structure_confidence"),
            "ptm": m.get("ptm"),
            "min_ipae": m.get("min_interaction_pae"),
            "cif": struct.get("url"),
        }
    return out

def mut_label(ext_id):
    # ext_id like "M23_H3_L6W" -> region/mutation
    parts = ext_id.split("_", 1)
    return parts[1] if len(parts) > 1 else ext_id

if __name__ == "__main__":
    raw = sys.argv[1] if len(sys.argv) > 1 else "results/maturation_raw.json"
    lib = {e["lib_id"]: e for e in json.load(open("design/library_ab2_maturation.json"))}
    sc = load_screen(raw)

    wt = sc.get(WT_ID, {})
    wt_bind = wt.get("bind_conf") or 0.0

    rows = []
    for ext_id, m in sc.items():
        lib_id = ext_id.split("_", 1)[0]
        e = lib.get(lib_id, {})
        rows.append({
            "ext_id": ext_id,
            "lib_id": lib_id,
            "mut": e.get("mut", mut_label(ext_id)),
            **m,
            "d_bind": (m["bind_conf"] - wt_bind) if isinstance(m.get("bind_conf"), (int, float)) else None,
            "seq": e.get("seq"),
        })

    rows.sort(key=lambda r: (-(r["bind_conf"] or -1), -(r["iptm"] or -1)))
    json.dump(rows, open("results/maturation_ranked.json", "w"), indent=1, default=str)

    print(f"WT ab2 binding_confidence = {wt_bind:.3f}\n")
    print(f"{'rank':>4} {'id':13} {'mut':10} {'bind':>5} {'Δbind':>6} {'iptm':>5} {'pIptm':>5} {'struct':>6} {'ptm':>5}")
    for i, r in enumerate(rows, 1):
        db = f"{r['d_bind']:+.3f}" if r["d_bind"] is not None else "   -  "
        pit = r["protein_iptm"]
        print(f"{i:>4} {r['lib_id']:13} {str(r['mut']):10} "
              f"{(r['bind_conf'] or 0):5.3f} {db:>6} {(r['iptm'] or 0):5.2f} "
              f"{(pit if pit is not None else 0):5.2f} {(r['struct_conf'] or 0):6.2f} {(r['ptm'] or 0):5.2f}")

    improved = [r for r in rows if r["d_bind"] and r["d_bind"] > 0 and r["lib_id"] != "M00"]
    print(f"\n{len(improved)} variants improved over ab2 WT (binding_confidence).")
