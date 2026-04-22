#!/usr/bin/env python3
"""Emit 11 x workflow inputs and 11 x Cromwell options JSONs for fg3 pQTL calibration + gwaslab mqq."""

from __future__ import annotations

import json
from pathlib import Path

# Chunk id → input_sumstats path list stem (as used under calibration_qc/input/)
CHUNKS = [
    ("01.chunk_0001_0500", "0001_0500", "0001-0500"),
    ("02.chunk_0501_1000", "0501_1000", "0501-1000"),
    ("03.chunk_1001_1500", "1001_1500", "1001-1500"),
    ("04.chunk_1501_2000", "1501_2000", "1501-2000"),
    ("05.chunk_2001_2500", "2001_2500", "2001-2500"),
    ("06.chunk_2501_3000", "2501_3000", "2501-3000"),
    ("07.chunk_3001_3500", "3001_3500", "3001-3500"),
    ("08.chunk_3501_4000", "3501_4000", "3501-4000"),
    ("09.chunk_4001_4500", "4001_4500", "4001-4500"),
    ("10.chunk_4501_5000", "4501_5000", "4501-5000"),
    ("11.chunk_5001_5416", "5001_5416", "5001-5416"),
]

BASE = "gs://r13-data/fg3_aggregate_pQTLs"
# Build and push this image (Dockerfile under ../docker) before submit.
DOCKER = (
    "europe-west1-docker.pkg.dev/finngen-refinery-dev/fg-refinery-registry/"
    "gwas-calibration-qc:v3"
)


def main() -> None:
    out = Path(__file__).resolve().parent
    for chunk_dir, file_tag, product_chunk in CHUNKS:
        label = f"input_sumstats_chunk_{file_tag}"
        inputs: dict = {
            "gwas_calibration_qc.paths_setup_a": (
                f"{BASE}/{chunk_dir}/calibration_qc/input/input_sumstats_paths_chunk_{file_tag}.txt"
            ),
            "gwas_calibration_qc.two_setups": False,
            "gwas_calibration_qc.setup_label_a": label,
            "gwas_calibration_qc.setup_label_b": "_unused",
            "gwas_calibration_qc.protein_id_mode": "first_segment",
            "gwas_calibration_qc.n_jobs": 2,
            "gwas_calibration_qc.memory_gb": 64,
            "gwas_calibration_qc.lead_window_bp": 1_500_000,
            "gwas_calibration_qc.diagnostic_plots": False,
            "gwas_calibration_qc.top_n_trans": 100,
            "gwas_calibration_qc.probability_rho": 1.0,
            "gwas_calibration_qc.gwaslab_mqq_plots": True,
            "gwas_calibration_qc.gwaslab_mqq_n_jobs": 2,
            "gwas_calibration_qc.gwaslab_mqq_build": "38",
            "gwas_calibration_qc.gwaslab_mqq_extra_disk_gb": 120,
            "gwas_calibration_qc.docker": DOCKER,
        }
        opts: dict = {
            "final_workflow_outputs_dir": f"{BASE}/{chunk_dir}/calibration_qc/cromwell_final_outputs",
            "final_workflow_outputs_mode": "copy",
            "google_labels": {
                "product": "core-analysis-r13",
                "fg3_pqtl_calib_mqq": "true",
                "chunk": product_chunk,
            },
        }
        stem = f"{chunk_dir}.calibration_qc"
        (out / f"{stem}.inputs.json").write_text(
            json.dumps(inputs, indent=2) + "\n", encoding="utf-8"
        )
        (out / f"{stem}.cromwell_workflow_options.json").write_text(
            json.dumps(opts, indent=2) + "\n", encoding="utf-8"
        )
    print(f"Wrote {2 * len(CHUNKS)} JSON files to {out}")


if __name__ == "__main__":
    main()
