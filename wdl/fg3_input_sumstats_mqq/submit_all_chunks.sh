#!/usr/bin/env bash
# Resubmit gwas_calibration_qc for all 11 fg3 pQTL chunks (input_sumstats path lists on gs://r13-data/fg3_aggregate_pQTLs/.../calibration_qc/input/)
#
# Prerequisite: build and push the Docker image referenced in the inputs (default tag v3), then run emit_chunk_input_jsons.py to refresh the JSONs if you change the tag.
# Prerequisite: SSH tunnel to the Cromwell server, e.g. cromwell_interact.py connect ... with SOCKS on port 5000.
#
# Usage: from this directory, after activating the same Python that has cromwell_interact dependencies,
#   bash submit_all_chunks.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$here/../.." && pwd)"
wdl="${repo_root}/wdl/gwas_calibration_qc.wdl"
interact="${repo_root}/../CromwellInteract/cromwell_interact.py"
port="${CROMWELL_SOCKS_PORT:-5000}"
http_port="${CROMWELL_HTTP_PORT:-80}"

if [[ ! -f "$interact" ]]; then
  echo "CromwellInteract not found at $interact" >&2
  exit 1
fi
if [[ ! -f "$wdl" ]]; then
  echo "WDL not found at $wdl" >&2
  exit 1
fi

shopt -s nullglob
inputs=( "$here"/*.calibration_qc.inputs.json )
if [[ ${#inputs[@]} -eq 0 ]]; then
  echo "No *calibration_qc.inputs.json in $here — run: python3 emit_chunk_input_jsons.py" >&2
  exit 1
fi

for inp in "${inputs[@]}"; do
  base="${inp%%.inputs.json}"
  opt="${base}.cromwell_workflow_options.json"
  label="$(basename "$base")"
  echo "Submitting $label"
  python3 "$interact" \
    --port "$port" \
    --http_port "$http_port" \
    submit \
    --wdl "$wdl" \
    --inputs "$inp" \
    --options "$opt" \
    --label "gwas_calib_mqq_${label}"
  sleep 2
done
echo "Done. Check workflows.log under CromwellInteract and Cromwell for job IDs."
