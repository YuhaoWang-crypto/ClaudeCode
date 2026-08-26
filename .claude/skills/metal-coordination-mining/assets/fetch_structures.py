#!/usr/bin/env python3
"""
fetch_structures — build the structure set that mcmine mines.

Two steps, each resumable and safe to re-run:

    accessions   InterPro / Pfam / SSF family  ->  UniProt accession table
    download     accession table               ->  <accession>.pdb from AFDB

Filtering between the two steps matters as much as the download: keep only
entries that HAVE an AlphaFold model and are long enough not to be fragments
(the published run used length > 150 aa), and deduplicate — one accession can be
annotated to several families and would otherwise be mined repeatedly.

Requires only the standard library plus pandas.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

INTERPRO_ROOT = "https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction"
HEADERS = {"Accept": "application/json", "User-Agent": "mcmine/1.0"}


def _family_url(family_id: str, page_size: int = 200) -> str:
    fid = family_id.strip()
    if fid.upper().startswith("IPR"):
        db = "InterPro"
    elif fid.upper().startswith("PF"):
        db = "PFAM"
    elif fid.upper().startswith("SSF"):
        db = "ssf"
    else:
        raise ValueError(f"Unrecognized family id {fid!r} (expected IPR…, PF…, SSF…)")
    return f"{INTERPRO_ROOT}/{db}/{fid}/?page_size={page_size}"


def _get_json(url: str, timeout: int = 60, retries: int = 3):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=timeout) as res:
                if getattr(res, "status", 200) == 204:
                    return None
                return json.loads(res.read().decode())
        except urllib.error.HTTPError as exc:
            if exc.code in (408, 429, 500, 502, 503):  # transient / throttled
                time.sleep(20 * (attempt + 1))
                continue
            raise
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            time.sleep(5 * (attempt + 1))
    return None


def cmd_accessions(args: argparse.Namespace) -> None:
    families = [f.strip() for f in args.families.split(",") if f.strip()]
    if args.families_file:
        families += [
            ln.strip() for ln in Path(args.families_file).read_text().splitlines() if ln.strip()
        ]
    if not families:
        raise SystemExit("Give --families and/or --families-file")

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    frames = []

    for family in families:
        out_csv = outdir / f"{family}_accessions.csv"
        if out_csv.exists() and out_csv.stat().st_size > 0 and not args.force:
            print(f"SKIP {family} (already retrieved)")
            frames.append(pd.read_csv(out_csv))
            continue

        rows, url, page = [], _family_url(family, args.page_size), 0
        while url:
            payload = _get_json(url)
            if not payload:
                break
            for item in payload.get("results", []):
                md = item.get("metadata", {})
                rows.append(
                    {
                        "AccessionID": md.get("accession"),
                        "Name": md.get("name"),
                        "Length": md.get("length"),
                        "Gene": md.get("gene"),
                        "in_alphafold": md.get("in_alphafold"),
                        "family": family,
                    }
                )
            url = payload.get("next")
            page += 1
            if page % 10 == 0:
                print(f"  {family}: {len(rows):,} accessions…", flush=True)
            time.sleep(args.delay)
            if args.max_pages and page >= args.max_pages:
                break

        df = pd.DataFrame(rows)
        df.to_csv(out_csv, index=False)
        print(f"{family}: {len(df):,} accessions -> {out_csv.name}")
        frames.append(df)

    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    combined = combined.dropna(subset=["AccessionID"])
    combined["Length"] = pd.to_numeric(combined["Length"], errors="coerce")
    combined["in_alphafold"] = (
        combined["in_alphafold"].astype(str).str.strip().str.lower().eq("true")
    )

    # one row per accession, keeping every family it came from
    grouped = (
        combined.groupby("AccessionID", as_index=False)
        .agg(
            Name=("Name", "first"),
            Length=("Length", "max"),
            in_alphafold=("in_alphafold", "any"),
            families=("family", lambda s: ";".join(sorted(set(s)))),
        )
    )
    kept = grouped[grouped["in_alphafold"] & (grouped["Length"] > args.min_length)]

    all_csv = outdir / "accessions_all.csv"
    kept_csv = outdir / f"accessions_af2_len_gt_{args.min_length}.csv"
    grouped.to_csv(all_csv, index=False)
    kept.to_csv(kept_csv, index=False)
    print(
        f"\nunique accessions: {len(grouped):,}\n"
        f"with AF2 model and length > {args.min_length}: {len(kept):,}\n"
        f"-> {kept_csv}"
    )


def _download_one(accession: str, outdir: Path, timeout: int) -> tuple[str, bool, str]:
    dest = outdir / f"{accession}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return accession, True, "cached"
    payload = _get_json(f"{AFDB_API}/{accession}", timeout=timeout)
    if not payload:
        return accession, False, "no AFDB record"
    record = payload[0] if isinstance(payload, list) and payload else payload
    url = record.get("pdbUrl") if isinstance(record, dict) else None
    if not url:
        return accession, False, "no pdbUrl"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": HEADERS["User-Agent"]})
        with urllib.request.urlopen(req, timeout=timeout) as res:
            data = res.read()
        if not data:
            return accession, False, "empty body"
        dest.write_bytes(data)
        return accession, True, "downloaded"
    except Exception as exc:  # network hiccup on one file must not kill the batch
        return accession, False, f"{type(exc).__name__}: {exc}"


def cmd_download(args: argparse.Namespace) -> None:
    src = Path(args.accessions).expanduser().resolve()
    if src.suffix.lower() == ".csv":
        df = pd.read_csv(src)
        column = args.column if args.column in df.columns else df.columns[0]
        accessions = df[column].astype(str).str.strip().tolist()
    else:
        accessions = [ln.strip() for ln in src.read_text().splitlines() if ln.strip()]
    if args.limit:
        accessions = accessions[: args.limit]

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {len(accessions):,} AFDB models -> {outdir}")

    ok = 0
    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(_download_one, acc, outdir, args.timeout): acc for acc in accessions
        }
        for n, future in enumerate(as_completed(futures), 1):
            acc, success, note = future.result()
            if success:
                ok += 1
            else:
                failures.append(f"{acc}\t{note}")
            if n % args.report_every == 0 or n == len(accessions):
                print(f"[{n}/{len(accessions)}] ok={ok} failed={len(failures)}", flush=True)

    log = outdir.parent / "download_failures.txt"
    log.write_text("\n".join(failures))
    print(f"\ndownloaded/cached: {ok:,}\nfailed: {len(failures):,} (see {log})")
    if failures and ok == 0:
        sys.exit(1)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="fetch_structures", description=__doc__.split("\n")[1])
    sub = p.add_subparsers(dest="command", required=True)

    sp = sub.add_parser("accessions", help="InterPro/Pfam/SSF family -> UniProt accessions")
    sp.add_argument("--families", default="", help="Comma list, e.g. PF02311,IPR013096")
    sp.add_argument("--families-file", help="File with one family id per line")
    sp.add_argument("--outdir", required=True)
    sp.add_argument("--min-length", type=int, default=150, help="Fragment cutoff (default 150 aa)")
    sp.add_argument("--page-size", type=int, default=200)
    sp.add_argument("--delay", type=float, default=0.5, help="Politeness delay between pages")
    sp.add_argument("--max-pages", type=int, default=None, help="Stop early (for a pilot run)")
    sp.add_argument("--force", action="store_true", help="Re-retrieve families already on disk")
    sp.set_defaults(func=cmd_accessions)

    sp = sub.add_parser("download", help="Accession list/CSV -> AFDB PDB files")
    sp.add_argument("--accessions", required=True, help="CSV or one-per-line text file")
    sp.add_argument("--column", default="AccessionID", help="Accession column if CSV")
    sp.add_argument("--outdir", required=True)
    sp.add_argument("--workers", type=int, default=8, help="Keep modest; AFDB throttles")
    sp.add_argument("--timeout", type=int, default=60)
    sp.add_argument("--limit", type=int, default=None)
    sp.add_argument("--report-every", type=int, default=500)
    sp.set_defaults(func=cmd_download)

    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    args.func(args)
