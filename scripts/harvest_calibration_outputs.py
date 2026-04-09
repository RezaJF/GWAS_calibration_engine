#!/usr/bin/env python3
# Reza Jabal; PhD; rjabal@broadinstitute.org
"""
Copy selected Cromwell workflow outputs to a destination GCS prefix (optional helper).

After a successful run, Cromwell delocalises outputs to the workflow execution
directory in the configured bucket. Use this script to copy the tarball and key
CSVs to a stable gs:// prefix for downstream consumers.

Requires gsutil.

Example:
  python3 scripts/harvest_calibration_outputs.py \\
    --source gs://cromwell-executions/.../gwas_calibration_qc/.../call-run_calibration_qc/ \\
    --dest gs://your-bucket/calibration/results/run_2026-04-08/
"""

from __future__ import annotations

import argparse
import subprocess
import sys


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True, help="GCS prefix to the run_calibration_qc shard output")
    p.add_argument("--dest", required=True, help="Destination gs:// prefix (trailing slash)")
    args = p.parse_args()
    src = args.source.rstrip("/") + "/"
    dest = args.dest.rstrip("/") + "/"
    cmd = ["gsutil", "-m", "cp", "-r", src + "*", dest]
    print("Running:", " ".join(cmd), flush=True)
    r = subprocess.run(cmd)
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
