#!/usr/bin/env python3
"""
Build paths_setup_a.txt and paths_setup_b.txt for gwas_calibration_qc.wdl.

Each output file must contain one gs:// path per line, in the same protein pairing
order implied by filenames (the QC script pairs by protein key from the stem).

This script lists *.tsv.gz, *.tsv, *.csv.gz, *.csv under two GCS prefixes using
gsutil. Requires google-cloud SDK (gsutil) on PATH.

Example:
  python3 scripts/generate_path_lists.py \\
    --prefix-a gs://bucket/run_a/regenie_summary_statistics/ \\
    --prefix-b gs://bucket/run_b/regenie_summary_statistics/ \\
    --out-a setup_a_paths.txt \\
    --out-b setup_b_paths.txt

Then upload the two .txt files to GCS and reference them in the Cromwell inputs JSON.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def gsutil_ls(prefix: str) -> list[str]:
    prefix = prefix.rstrip("/") + "/"
    if not prefix.startswith("gs://"):
        print("Expected gs:// prefix", file=sys.stderr)
        sys.exit(1)
    r = subprocess.run(
        ["gsutil", "-m", "ls", "-r", prefix + "**"],
        capture_output=True,
        text=True,
        check=False,
    )
    if r.returncode != 0:
        r = subprocess.run(
            ["gsutil", "ls", prefix + "*.gz"],
            capture_output=True,
            text=True,
            check=False,
        )
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        sys.exit(r.returncode)
    lines = [ln.strip() for ln in r.stdout.splitlines() if ln.strip()]
    allowed = (".tsv.gz", ".tsv", ".csv.gz", ".csv")
    out = [ln for ln in lines if ln.endswith(allowed) and not ln.endswith("/")]
    out.sort()
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--prefix-a", required=True, help="GCS prefix for setup A sumstats")
    p.add_argument("--prefix-b", required=True, help="GCS prefix for setup B sumstats")
    p.add_argument("--out-a", type=Path, default=Path("setup_a_paths.txt"))
    p.add_argument("--out-b", type=Path, default=Path("setup_b_paths.txt"))
    args = p.parse_args()

    a = gsutil_ls(args.prefix_a)
    b = gsutil_ls(args.prefix_b)
    args.out_a.write_text("\n".join(a) + ("\n" if a else ""))
    args.out_b.write_text("\n".join(b) + ("\n" if b else ""))
    print(f"Wrote {len(a)} paths to {args.out_a}")
    print(f"Wrote {len(b)} paths to {args.out_b}")


if __name__ == "__main__":
    main()
