# fg3 pQTL calibration + GWASlab mqq — chunk 01..11 (`input_sumstats`)

Submission pack for the **single-task** workflow
[`wdl/gwas_calibration_qc.wdl`](../gwas_calibration_qc.wdl), running over the **11 fg3 pQTL
chunks** whose per-trait summary statistics live under:

```
gs://r13-data/fg3_aggregate_pQTLs/<NN>.chunk_<start>_<end>/input_sumstats/*.gz
```

Each `gwas_calibration_qc.paths_setup_a` points at the pre-built path list

```
gs://r13-data/fg3_aggregate_pQTLs/<NN>.chunk_<start>_<end>/calibration_qc/input/input_sumstats_paths_chunk_<start>_<end>.txt
```

one `gs://` URI per line, which Cromwell localises for the calibration + mqq task.

## Layout

**Tracked (reusable assets):**

- `emit_chunk_input_jsons.py` — regenerates 22 JSON files from one constant table; edit `DOCKER` / `BASE` / `CHUNKS` to adapt to a different batch.
- `submit_all_chunks.sh` — loops over every `*.calibration_qc.inputs.json` and submits via `CromwellInteract/cromwell_interact.py submit`.
- `README.md` — this file.

**Ignored (run-specific, regenerable):**

- `<NN>.chunk_<start>_<end>.calibration_qc.inputs.json` ×11 — produced by `emit_chunk_input_jsons.py`; embeds internal bucket paths.
- `<NN>.chunk_<start>_<end>.calibration_qc.cromwell_workflow_options.json` ×11 — per-chunk `final_workflow_outputs_dir`, Google labels.
- `submit_all_chunks.log*` — stdout captured from a submission run.
- `SUBMISSION_REPORT*.md` — per-submission report with workflow IDs and image digests.

## Pre-flight (on your workstation / tunnel host)

1. **Build and push** the docker image (tag `v3` is referenced by the JSONs). Requirements now include **`gwaslab>=4.0.8`** and **matplotlib** in both
   `Github_clones/GWAS_calibration_engine/requirements.txt` and
   `Python_scripts/gwas_calibration_utils/pyproject.toml`.

```
cd REPO_ROOT  # parent of Github_clones/
docker build \
  -f Github_clones/GWAS_calibration_engine/docker/Dockerfile \
  -t europe-west1-docker.pkg.dev/finngen-refinery-dev/fg-refinery-registry/gwas-calibration-qc:v3 \
  .
docker push  europe-west1-docker.pkg.dev/finngen-refinery-dev/fg-refinery-registry/gwas-calibration-qc:v3
```

2. **Tunnel** to the Cromwell server (SOCKS on port 5000):

```
python3 Github_clones/CromwellInteract/cromwell_interact.py connect <server>
```

3. **Submit all 11 chunks**:

```
bash Github_clones/GWAS_calibration_engine/wdl/fg3_input_sumstats_mqq/submit_all_chunks.sh
```

The submit script uses `--options <...cromwell_workflow_options.json>` (carries
`google_labels.product=core-analysis-r13`) and writes workflow IDs into
`CromwellInteract/workflows.log`.

## Outputs per chunk

`final_workflow_outputs_dir` is
`gs://r13-data/fg3_aggregate_pQTLs/<chunk>/calibration_qc/cromwell_final_outputs/`.
Artefacts include the standard calibration tarball (metrics CSV, summary JSON,
tiered report, log) and, with `gwaslab_mqq_plots = true` (default), a set of
`results/gwaslab_mqq/<trait>_gwaslab_mqq.png` panels — each a GWASlab
`mode="mqq"` figure with **Manhattan + MAF-stratified QQ** using
`stratified=True`, `cut=14`, `skip=3`, `marker_size=(5, 10)`.

## Changing docker tag / resources

Edit the constants at the top of `emit_chunk_input_jsons.py` (`DOCKER`, etc.)
and rerun it — the 22 JSONs will be overwritten in place.
