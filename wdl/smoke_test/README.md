# WDL smoke test — 2 files

Runs the full [`gwas_calibration_qc.wdl`](../gwas_calibration_qc.wdl) pipeline
(calibration metrics **plus** the GWASlab mqq step) on **two** sumstats files
so that regressions — especially Cromwell/WDL localization regressions like
the `write_lines` declaration-time serialization bug — surface in < 10 minutes
instead of failing an 11-chunk production batch.

## What it catches

- **`write_lines` declaration-time evaluation** that serialises `gs://` URIs
  instead of local worker paths (the exact bug that killed round 1 of the
  full submission). Since 2026-04-22, `gwas_calibration_qc.py:resolve_input`
  hard-fails with exit code 2 when it sees any cloud URI prefix in a file
  list (unless `--allow-remote-paths` is passed). This test will terminate
  quickly and unambiguously when that bug returns.
- Missing Python modules in the image (e.g. `gwaslab`, `matplotlib`).
- Dockerfile COPY regressions that leave scripts out of `/opt/gwas_calibration_utils/`.
- `sep` / `build` / column-detection regressions in `gwaslab_mqq_plots.py`.
- Basic metadata-harvesting regressions (the `results/` tarball must contain
  all 5 declared outputs).

## What it does NOT catch

- Batch-global statistical properties (with only 2 traits, λ, tail excess,
  GOF, rare-vs-common metrics are noisy). Use the full-chunk batch for that.
- Scatter–gather WDL (separate smoke test needed for
  `gwas_calibration_qc_scattered.wdl`).

## Inputs (templates)

Two template JSONs are tracked in this directory:

- `smoke_2files.calibration_qc.inputs.example.json` — workflow inputs with
  `gs://your-bucket/...` and image-tag placeholders.
- `smoke_2files.cromwell_workflow_options.example.json` — Cromwell options
  (`final_workflow_outputs_dir`, `google_labels`) with placeholders.

Make your own real-bucket copies (these file names are ignored by git):

```bash
cp wdl/smoke_test/smoke_2files.calibration_qc.inputs.example.json \
   wdl/smoke_test/smoke_2files.calibration_qc.inputs.json
cp wdl/smoke_test/smoke_2files.cromwell_workflow_options.example.json \
   wdl/smoke_test/smoke_2files.cromwell_workflow_options.json
# then edit the two new files to replace placeholders with your bucket / image tag
```

Upload (once per test fixture) a 2-line path list that references two small
sumstats your org already has:

```bash
# Example: take the first two URIs from a production chunk-01 path list
gsutil cat gs://YOUR-BUCKET/<...>/input_sumstats_paths_chunk_0001_0500.txt \
  | head -n 2 > /tmp/smoke_2files.txt
gsutil cp /tmp/smoke_2files.txt \
  gs://YOUR-BUCKET/tmp/gwas_calibration_smoke/input_sumstats_paths_smoke_2files.txt
```

## Run

```bash
# 1. Make sure the tunnel is up (port 5000 → Cromwell).
# 2. Submit:
python3 Github_clones/CromwellInteract/cromwell_interact.py \
  --port 5000 --http_port 80 \
  submit \
    --wdl Github_clones/GWAS_calibration_engine/wdl/gwas_calibration_qc.wdl \
    --inputs Github_clones/GWAS_calibration_engine/wdl/smoke_test/smoke_2files.calibration_qc.inputs.json \
    --options Github_clones/GWAS_calibration_engine/wdl/smoke_test/smoke_2files.cromwell_workflow_options.json \
    --label gwas_calib_mqq_smoke_2files
```

Expected wall time: **5–8 min** on a 2-core, 16 GB preemptible VM (includes
image pull + 2-file localization + 2 calibration workers + 2 mqq workers +
tarball upload). Expected artefacts in
`$final_workflow_outputs_dir`:

- `calibration_outputs.tgz`
  - `results/calibration_compare.metrics.long.csv`  (2 proteins × metrics)
  - `results/calibration_compare.summary.json`
  - `results/calibration_compare.tiered_report.txt`
  - `results/gwas_calibration_qc.log`
  - `results/gwaslab_mqq/A1BG_cross_batch_normalised_4144_gwaslab_mqq.png`
  - `results/gwaslab_mqq/A1CF_cross_batch_normalised_4144_gwaslab_mqq.png`

## Pass criteria (for automation)

```bash
# After workflow status == Succeeded:
WF=<workflow_id>
gsutil cat "gs://cromwell-fg-3/gwas_calibration_qc/${WF}/call-run_calibration_one_setup/stderr" \
  | grep -qE '^ERROR: No input files resolved' \
  && { echo "REGRESSION: write_lines localization bug is back"; exit 1; } || true

gsutil ls "$final_workflow_outputs_dir" | \
  grep -q calibration_outputs.tgz || { echo "Missing tarball"; exit 1; }
```
