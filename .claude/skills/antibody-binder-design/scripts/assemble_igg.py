#!/usr/bin/env python3
"""Assemble a full humanized IgG from VH/VL, and emit the Fab chains used for co-fold.

Usage:
  python assemble_igg.py --vh VHSEQ --vl VLSEQ --light lambda|kappa [--name AB] [--fasta out.fasta] [--json out.json]

Constant regions are human IgG1 CH1-hinge-CH2-CH3, human C-lambda and C-kappa.
Emits: full heavy chain, full light chain, and the Fab chains (Fd=VH+CH1, light=VL+CL)
which are what you co-fold with the antigen (the Fc never contacts antigen).
"""
import argparse, json, sys

CH1 = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV"
HINGE_CH2_CH3 = ("EPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEY"
                 "KCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK")
CL = {
 "lambda": "QPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS",
 "kappa":  "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC",
}

def assemble(vh, vl, light="lambda"):
    light = light.lower()
    if light not in CL:
        raise SystemExit(f"--light must be lambda or kappa, got {light}")
    fd = vh + CH1
    lc = vl + CL[light]
    return {
        "heavy_chain_IgG1": vh + CH1 + HINGE_CH2_CH3,
        "light_chain": lc,
        "fab_Fd": fd,          # VH-CH1  -> co-fold chain H
        "fab_light": lc,       # VL-CL   -> co-fold chain L
        "light_type": light,
        "lengths": {"VH": len(vh), "VL": len(vl), "Fd": len(fd), "LC": len(lc),
                    "HC_full": len(vh)+len(CH1)+len(HINGE_CH2_CH3)},
    }

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--vh", required=True); ap.add_argument("--vl", required=True)
    ap.add_argument("--light", default="lambda"); ap.add_argument("--name", default="ab")
    ap.add_argument("--fasta"); ap.add_argument("--json")
    a = ap.parse_args()
    r = assemble(a.vh.strip().upper(), a.vl.strip().upper(), a.light)
    r["name"] = a.name
    if a.fasta:
        with open(a.fasta, "w") as f:
            f.write(f">{a.name}_HeavyChain_IgG1 {r['lengths']['HC_full']}aa\n{r['heavy_chain_IgG1']}\n")
            f.write(f">{a.name}_LightChain_{r['light_type']} {r['lengths']['LC']}aa\n{r['light_chain']}\n")
    if a.json: json.dump(r, open(a.json, "w"), indent=1)
    print(json.dumps({k: (v if not isinstance(v, str) or len(v) < 40 else v[:20]+f"...({len(v)}aa)") for k, v in r.items()}, indent=1))
