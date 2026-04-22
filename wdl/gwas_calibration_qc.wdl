version 1.0

# =============================================================================
# GWAS calibration QC — one or two setups
# =============================================================================
#
# Runs gwas_calibration_qc.py on Cromwell (GCP). Sumstats are listed as gs://
# paths (one per line per setup). Use two_setups=false for a single batch (e.g.
# one intermediary/ folder). Set two_setups=true and provide paths_setup_b for
# A-vs-B comparison.
# =============================================================================

workflow gwas_calibration_qc {

  input {
    File    paths_setup_a
    File?   paths_setup_b
    Boolean two_setups = true

    # Optional. If omitted, the task passes --lead-variants-json results/lead_variants.json; Python
    # creates it from localised sumstats (genome-wide min p per protein). Supply a File to use a fixed JSON.
    File?   lead_variants_json
    File?   cis_json

    String setup_label_a = "setup_A"
    String setup_label_b = "setup_B"

    String protein_id_mode = "first_segment"
    Int    n_jobs = 4
    Int    memory_gb = 16
    Boolean diagnostic_plots = true
    Int    top_n_trans = 100
    Float  probability_rho = 1.0
    # Half-width around each lead (bp): exclude |pos - lead_pos| <= this on the lead chromosome (default ±1.5 Mb → 3 Mb cis-like span).
    Int    lead_window_bp = 1500000

    # gwaslab Manhattan+MAF-stratified QQ (mqq) figures for each localised sumstats file; written under results/gwaslab_mqq/
    Boolean gwaslab_mqq_plots = true
    Int     gwaslab_mqq_n_jobs = 2
    String  gwaslab_mqq_build = "38"
    # Extra local SSD when gwaslab mqq is enabled (plots are I/O- and memory-heavy).
    Int     gwaslab_mqq_extra_disk_gb = 120

    String docker
  }

  if (two_setups) {
    call run_calibration_two_setups {
      input:
        paths_setup_a      = paths_setup_a,
        paths_setup_b      = select_first([paths_setup_b]),
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
        gwaslab_mqq_plots   = gwaslab_mqq_plots,
        gwaslab_mqq_n_jobs  = gwaslab_mqq_n_jobs,
        gwaslab_mqq_build  = gwaslab_mqq_build,
        gwaslab_mqq_extra_disk_gb = gwaslab_mqq_extra_disk_gb,
        docker             = docker
    }
  }

  if (!two_setups) {
    call run_calibration_one_setup {
      input:
        paths_setup_a      = paths_setup_a,
        lead_variants_json = lead_variants_json,
        cis_json           = cis_json,
        setup_label_a      = setup_label_a,
        protein_id_mode    = protein_id_mode,
        n_jobs             = n_jobs,
        memory_gb          = memory_gb,
        diagnostic_plots   = diagnostic_plots,
        top_n_trans        = top_n_trans,
        probability_rho    = probability_rho,
        lead_window_bp     = lead_window_bp,
        gwaslab_mqq_plots   = gwaslab_mqq_plots,
        gwaslab_mqq_n_jobs  = gwaslab_mqq_n_jobs,
        gwaslab_mqq_build  = gwaslab_mqq_build,
        gwaslab_mqq_extra_disk_gb = gwaslab_mqq_extra_disk_gb,
        docker             = docker
    }
  }

  output {
    File calibration_outputs_tgz = select_first([run_calibration_one_setup.calibration_outputs_tgz, run_calibration_two_setups.calibration_outputs_tgz])
    File metrics_long_csv        = select_first([run_calibration_one_setup.metrics_long_csv, run_calibration_two_setups.metrics_long_csv])
    File summary_json            = select_first([run_calibration_one_setup.summary_json, run_calibration_two_setups.summary_json])
    File tiered_report_txt       = select_first([run_calibration_one_setup.tiered_report_txt, run_calibration_two_setups.tiered_report_txt])
    File qc_log                  = select_first([run_calibration_one_setup.qc_log, run_calibration_two_setups.qc_log])
  }
}


task run_calibration_one_setup {

  input {
    File    paths_setup_a
    File?   lead_variants_json
    File?   cis_json

    String  setup_label_a
    String  protein_id_mode
    Int     n_jobs
    Int     memory_gb
    Boolean diagnostic_plots
    Int     top_n_trans
    Float   probability_rho
    Int     lead_window_bp
    Boolean gwaslab_mqq_plots
    Int     gwaslab_mqq_n_jobs
    String  gwaslab_mqq_build
    Int     gwaslab_mqq_extra_disk_gb
    String  docker
  }

  Array[File] files_a = read_lines(paths_setup_a)

  Int disk_pad_gb = 20
  Int disk_gb = if gwaslab_mqq_plots then ceil(size(files_a, "GB") + disk_pad_gb + gwaslab_mqq_extra_disk_gb) else ceil(size(files_a, "GB") + disk_pad_gb)

  command <<<
    set -euo pipefail

    mkdir -p results
    LST="~{write_lines(files_a)}"

    gwas-calibration-qc \
      --file-list "$LST" \
      ~{if defined(lead_variants_json) then "--lead-variants-json " + lead_variants_json else "--lead-variants-json results/lead_variants.json"} \
      ~{if defined(cis_json) then "--cis-json " + cis_json else ""} \
      --setup-labels "~{setup_label_a}" "_unused" \
      --protein-id-mode ~{protein_id_mode} \
      --outdir results \
      --n-jobs ~{n_jobs} \
      ~{if diagnostic_plots then "--diagnostic-plots" else ""} \
      --top-n-trans ~{top_n_trans} \
      --probability-rho ~{probability_rho} \
      --lead-window-bp ~{lead_window_bp}

    if ~{gwaslab_mqq_plots}; then
      mkdir -p results/gwaslab_mqq
      gwaslab-mqq-plots \
        --file-list "$LST" \
        --outdir results/gwaslab_mqq \
        --n-jobs ~{gwaslab_mqq_n_jobs} \
        --build ~{gwaslab_mqq_build} \
        --log-level INFO
    fi

    tar -czf calibration_outputs.tgz -C results .
  >>>

  output {
    File calibration_outputs_tgz = "calibration_outputs.tgz"
    File metrics_long_csv        = "results/calibration_compare.metrics.long.csv"
    File summary_json            = "results/calibration_compare.summary.json"
    File tiered_report_txt       = "results/calibration_compare.tiered_report.txt"
    File qc_log                  = "results/gwas_calibration_qc.log"
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


task run_calibration_two_setups {

  input {
    File    paths_setup_a
    File    paths_setup_b
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
    Boolean gwaslab_mqq_plots
    Int     gwaslab_mqq_n_jobs
    String  gwaslab_mqq_build
    Int     gwaslab_mqq_extra_disk_gb
    String  docker
  }

  Array[File] files_a = read_lines(paths_setup_a)
  Array[File] files_b = read_lines(paths_setup_b)

  Int disk_pad_gb = 20
  Int disk_extra_mqq = if gwaslab_mqq_plots then gwaslab_mqq_extra_disk_gb else 0
  Int disk_gb = ceil(size(files_a, "GB") + size(files_b, "GB") + disk_pad_gb + disk_extra_mqq)

  command <<<
    set -euo pipefail

    mkdir -p results
    LSTA="~{write_lines(files_a)}"
    LSTB="~{write_lines(files_b)}"

    gwas-calibration-qc \
      --file-list "$LSTA" \
      --file-list "$LSTB" \
      ~{if defined(lead_variants_json) then "--lead-variants-json " + lead_variants_json else "--lead-variants-json results/lead_variants.json"} \
      ~{if defined(cis_json) then "--cis-json " + cis_json else ""} \
      --setup-labels "~{setup_label_a}" "~{setup_label_b}" \
      --protein-id-mode ~{protein_id_mode} \
      --outdir results \
      --n-jobs ~{n_jobs} \
      ~{if diagnostic_plots then "--diagnostic-plots" else ""} \
      --top-n-trans ~{top_n_trans} \
      --probability-rho ~{probability_rho} \
      --lead-window-bp ~{lead_window_bp}

    if ~{gwaslab_mqq_plots}; then
      mkdir -p results/gwaslab_mqq
      gwaslab-mqq-plots \
        --file-list "$LSTA" \
        --outdir results/gwaslab_mqq \
        --n-jobs ~{gwaslab_mqq_n_jobs} \
        --build ~{gwaslab_mqq_build} \
        --log-level INFO
    fi

    tar -czf calibration_outputs.tgz -C results .
  >>>

  output {
    File calibration_outputs_tgz = "calibration_outputs.tgz"
    File metrics_long_csv        = "results/calibration_compare.metrics.long.csv"
    File summary_json            = "results/calibration_compare.summary.json"
    File tiered_report_txt       = "results/calibration_compare.tiered_report.txt"
    File qc_log                  = "results/gwas_calibration_qc.log"
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
