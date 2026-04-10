version 1.0

# =============================================================================
# GWAS calibration QC — scatter / gather (maximum Cromwell parallelism)
# =============================================================================
#
# Splits master path list(s) into N shards, runs --phase scatter per shard in
# parallel, then --phase gather to compute batch-global calibration_score and
# quality_prob_good (median/MAD and GOF min–max require the full batch).
# =============================================================================

workflow gwas_calibration_qc_scattered {

  input {
    File    master_paths_setup_a
    File?   master_paths_setup_b
    Boolean two_setups = true

    Int     n_shards = 10

    File?   lead_variants_json
    File?   cis_json

    String  setup_label_a = "setup_A"
    String  setup_label_b = "setup_B"

    String  protein_id_mode = "first_segment"
    Int     n_jobs = 4
    Int     memory_gb = 16
    Boolean diagnostic_plots = true
    Int     top_n_trans = 100
    Float   probability_rho = 1.0
    Int     lead_window_bp = 1500000

    String docker

    Int split_cpu = 1
    Int split_memory_mb = 1024
    Int gather_cpu = 1
    Int gather_memory_gb = 4
  }

  Int two_setups_int = if (two_setups) then 1 else 0

  call split_path_lists {
    input:
      master_paths_a   = master_paths_setup_a,
      master_paths_b   = master_paths_setup_b,
      n_shards         = n_shards,
      two_setups_int   = two_setups_int,
      cpu              = split_cpu,
      memory_mb        = split_memory_mb,
      docker_image     = docker
  }

  scatter (idx in range(length(split_path_lists.shard_paths_a))) {
    call compute_shard_metrics {
      input:
        paths_setup_a      = split_path_lists.shard_paths_a[idx],
        paths_setup_b      = split_path_lists.shard_paths_b[idx],
        two_setups         = two_setups,
        lead_variants_json = lead_variants_json,
        cis_json           = cis_json,
        setup_label_a      = setup_label_a,
        setup_label_b      = setup_label_b,
        protein_id_mode    = protein_id_mode,
        n_jobs             = n_jobs,
        memory_gb          = memory_gb,
        diagnostic_plots   = diagnostic_plots,
        top_n_trans        = top_n_trans,
        probability_rho    = probability_rho,
        lead_window_bp     = lead_window_bp,
        docker             = docker
    }
  }

  call aggregate_and_score {
    input:
      partial_csvs    = compute_shard_metrics.partial_metrics_csv,
      setup_label_a   = setup_label_a,
      setup_label_b   = setup_label_b,
      probability_rho = probability_rho,
      gather_cpu      = gather_cpu,
      gather_memory_gb = gather_memory_gb,
      docker          = docker
  }

  output {
    Array[File] shard_partial_metrics_csv = compute_shard_metrics.partial_metrics_csv
    Array[File] shard_logs                = compute_shard_metrics.shard_log
    Array[File] shard_plots_tgz           = compute_shard_metrics.shard_plots_tgz
    Array[File] shard_lead_variants_json  = compute_shard_metrics.shard_lead_variants_json

    File metrics_long_csv  = aggregate_and_score.metrics_long_csv
    File summary_json      = aggregate_and_score.summary_json
    File tiered_report_txt = aggregate_and_score.tiered_report_txt
    File pivot_csv         = aggregate_and_score.pivot_csv
    File gather_log        = aggregate_and_score.gather_log
  }
}


task split_path_lists {

  input {
    File  master_paths_a
    File? master_paths_b
    Int   n_shards
    Int   two_setups_int

    Int    cpu
    Int    memory_mb
    String docker_image
  }

  command <<<
    set -euo pipefail
    python3 << 'PY'
    import math

    n_shards_req = max(1, int("~{n_shards}"))
    path_a = r"~{master_paths_a}"
    path_b_str = r"~{if defined(master_paths_b) then master_paths_b else ""}"
    two = int("~{two_setups_int}") == 1

    with open(path_a, encoding="utf-8") as f:
        lines_a = [ln.rstrip("\n") for ln in f if ln.strip() and not ln.strip().startswith("#")]

    lines_b = None
    if two:
        if not path_b_str:
            raise SystemExit("two_setups is true but master_paths_setup_b was not provided.")
        with open(path_b_str, encoding="utf-8") as f:
            lines_b = [ln.rstrip("\n") for ln in f if ln.strip() and not ln.strip().startswith("#")]
        if len(lines_b) != len(lines_a):
            raise SystemExit(
                "Path lists must have the same length for paired two-setup mode: "
                f"{len(lines_a)} vs {len(lines_b)}."
            )

    L = len(lines_a)
    if L == 0:
        raise SystemExit("Master path list A is empty.")

    n_eff = min(n_shards_req, L)
    shard_size = math.ceil(L / n_eff)

    for k in range(n_eff):
        start = k * shard_size
        end = min(L, start + shard_size)
        chunk_a = lines_a[start:end]
        with open(f"shard_a_{k:04d}.txt", "w", encoding="utf-8") as out:
            out.write("\n".join(chunk_a) + ("\n" if chunk_a else ""))
        if lines_b is not None:
            chunk_b = lines_b[start:end]
            with open(f"shard_b_{k:04d}.txt", "w", encoding="utf-8") as out:
                out.write("\n".join(chunk_b) + ("\n" if chunk_b else ""))

    # Single-setup scatter still passes a dummy second list to satisfy File typing.
    if lines_b is None:
        with open("shard_b_0000.txt", "w", encoding="utf-8") as out:
            out.write("")
        for k in range(1, n_eff):
            with open(f"shard_b_{k:04d}.txt", "w", encoding="utf-8") as out:
                out.write("")
    PY
  >>>

  output {
    Array[File] shard_paths_a = glob("shard_a_*.txt")
    Array[File] shard_paths_b = glob("shard_b_*.txt")
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: "~{memory_mb} MB"
    disks: "local-disk 10 SSD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 1
    noAddress: true
  }
}


task compute_shard_metrics {

  input {
    File    paths_setup_a
    File    paths_setup_b
    Boolean two_setups

    File?   lead_variants_json
    File?   cis_json

    String  setup_label_a
    String  setup_label_b
    String  protein_id_mode
    Int     n_jobs
    Int     memory_gb
    Boolean diagnostic_plots
    Int     top_n_trans
    Float   probability_rho
    Int     lead_window_bp
    String  docker
  }

  Array[File] files_a = read_lines(paths_setup_a)
  Array[File] files_b = if (two_setups) then read_lines(paths_setup_b) else []

  Int disk_pad_gb = 20
  Int disk_gb = if (two_setups) then ceil(size(files_a, "GB") + size(files_b, "GB") + disk_pad_gb) else ceil(size(files_a, "GB") + disk_pad_gb)

  command <<<
    set -euo pipefail

    mkdir -p results

    python3 /opt/gwas_calibration_utils/gwas_calibration_qc.py \
      --phase scatter \
      --file-list ~{write_lines(files_a)} \
      ~{if (two_setups) then "--file-list " + write_lines(files_b) else ""} \
      ~{if defined(lead_variants_json) then "--lead-variants-json " + lead_variants_json else "--lead-variants-json results/lead_variants.json"} \
      ~{"--cis-json " + cis_json} \
      --setup-labels "~{setup_label_a}" "~{setup_label_b}" \
      --protein-id-mode ~{protein_id_mode} \
      --outdir results \
      --n-jobs ~{n_jobs} \
      ~{if diagnostic_plots then "--diagnostic-plots" else ""} \
      --top-n-trans ~{top_n_trans} \
      --probability-rho ~{probability_rho} \
      --lead-window-bp ~{lead_window_bp}

    tar -czf shard_plots.tgz -C results .

    if [ -f results/lead_variants.json ]; then
      cp results/lead_variants.json shard_lead_variants.json
    else
      echo '{}' > shard_lead_variants.json
    fi

    cp results/gwas_calibration_qc.log shard_gwas_calibration_qc.log
  >>>

  output {
    File partial_metrics_csv     = "results/calibration_compare.partial_metrics.csv"
    File shard_log               = "shard_gwas_calibration_qc.log"
    File shard_plots_tgz         = "shard_plots.tgz"
    File shard_lead_variants_json = "shard_lead_variants.json"
  }

  runtime {
    docker: docker
    cpu: n_jobs
    memory: "~{memory_gb} GB"
    disks: "local-disk " + disk_gb + " SSD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
    noAddress: true
  }
}


task aggregate_and_score {

  input {
    Array[File] partial_csvs
    String        setup_label_a
    String        setup_label_b
    Float         probability_rho

    Int    gather_cpu
    Int    gather_memory_gb
    String docker
  }

  command <<<
    set -euo pipefail

    mkdir -p results

    python3 /opt/gwas_calibration_utils/gwas_calibration_qc.py \
      --phase gather \
      --partial-csv-list ~{write_lines(partial_csvs)} \
      --outdir results \
      --probability-rho ~{probability_rho} \
      --setup-labels "~{setup_label_a}" "~{setup_label_b}"

    if [ ! -f results/calibration_compare.metrics.pivot.csv ]; then
      echo "# single_setup_or_no_paired_pivot" > results/calibration_compare.metrics.pivot.csv
    fi

    cp results/gwas_calibration_qc.log gather_gwas_calibration_qc.log
  >>>

  output {
    File metrics_long_csv  = "results/calibration_compare.metrics.long.csv"
    File summary_json      = "results/calibration_compare.summary.json"
    File tiered_report_txt = "results/calibration_compare.tiered_report.txt"
    File pivot_csv         = "results/calibration_compare.metrics.pivot.csv"
    File gather_log        = "gather_gwas_calibration_qc.log"
  }

  runtime {
    docker: docker
    cpu: gather_cpu
    memory: "~{gather_memory_gb} GB"
    disks: "local-disk 10 SSD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 1
    noAddress: true
  }
}
