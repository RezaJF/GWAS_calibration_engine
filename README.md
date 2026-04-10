# GWAS calibration QC — Cromwell workflow

A reproducible **[WDL](https://openwdl.org/)** workflow for large-scale **statistical calibration quality control** of genome-wide association study (GWAS) or protein quantitative trait locus (pQTL) **summary statistics** on **Google Cloud Platform**. The workflow localises per-trait sumstats from object storage, runs a containerised calibration analysis in parallel where configured, and packages metrics, reports, and optional diagnostic figures for downstream triage or method comparison.

**Licence:** The **WDL**, **Dockerfile**, **helper scripts**, and **documentation** in this repository are released under the **[MIT Licence](LICENSE)**. The bundled **gwas-calibration-utils** dependency remains under its own **MIT** terms (see that package’s metadata).

This repository provides **orchestration only** (WDL, Docker build metadata, and operational examples). The numerical methods are implemented in the **gwas-calibration-utils** Python package, which is installed inside the workflow container image.

---

## Purpose and scope

**Calibration** here means checking whether *p*-values on variants that are **not** expected to carry a true genetic signal behave like **null *p*-values** — i.e. approximately **uniform on (0, 1)** after appropriate masking. Strong cis associations and lead signals violate that null on nearby variants; the analysis therefore **masks** (excludes) configurable **cis regions** (from optional JSON) and a **symmetric window around each lead variant** (workflow default **±1.5 Mb** on the lead chromosome, **3 Mb** total — input **`lead_window_bp`**) before evaluating the *trans* remainder.

Typical uses:

- **Batch QC** across hundreds or thousands of traits after a production GWAS or pQTL pipeline.
- **Comparing two processing setups** (e.g. normalisation or batch-correction choices) on the same traits using paired path lists and shared masking configuration.

The workflow does **not** replace association testing, fine-mapping, or colocalisation; it **summarises distributional behaviour** of *p*-values on the masked SNP set and emits ranked **composite calibration scores** for relative ranking within a run.

---

## Scientific rationale

Under the **global null** for a given SNP (no association), the two-sided test *p*-value is (asymptotically) **Uniform(0, 1)**. Transforming *p* to a one-degree-of-freedom chi-square statistic,

$$
X = F_{\chi^2(1)}^{-1}(1 - p),
$$

gives **$X \sim \chi^2(1)$** under that null. **Genomic control** and **quantile inflation** summarise whether realised *X* (or equivalently *p*) is **stochastically larger** than the null (inflation) or **smaller** (deflation).

After **removing cis and lead-adjacent variants**, the remaining SNPs are treated as **null-enriched** for the purpose of calibration diagnostics: under a well-calibrated pipeline, their *p*-values should not show **systematic** bulk or **tail** inflation relative to Uniform(0, 1). Deviations often indicate **misspecification**, **residual confounding**, **batch artefacts**, **allele-frequency quirks**, or **overfitting** in secondary processing — hence the value of a standardised QC pass before release or cross-cohort meta-analysis.

---

## Mathematical foundation (summary)

All metrics below are computed on the **masked** SNP table for each trait (same chromosome/position conventions as the input files).


| Family                  | Construction                                                                                                                                                                   | Null reference                                                  |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------- |
| **Genomic control λGC** | Let $X_i$ be $\chi^2(1)$ transforms of *p*-values. $\lambda_{\mathrm{GC}} = \mathrm{median}_i(X_i) / \mathrm{median}(\chi^2(1))$. The median of $\chi^2(1)$ is ≈ **0.4549**.   | $\lambda_{\mathrm{GC}} \approx 1$                               |
| **Quantile λ**          | For quantile *q*, ratio of the empirical quantile of $X_i$ to the theoretical $\chi^2(1)$ quantile at *q* (e.g. *q* ∈ {0.5, 0.9, 0.99, 0.999}).                                | ≈ 1 at each *q*                                                 |
| **Tail excess (K/E)**   | At threshold α, **K** = count of SNPs with *p* ≤ α; **E** = *m*α for *m* SNPs. Report **K/E** for several α (e.g. 10⁻⁴ … 10⁻⁷).                                                | **K/E** ≈ 1                                                     |
| **Goodness-of-fit**     | Histogram test of *p* against Uniform(0, 1) on a fixed bin grid (including a fine bin near 0). Yields a $\chi^2$ statistic and *p*-value for global departure from uniformity. | Non-small GOF *p*                                               |
| **Rare vs common**      | Same K/E-style tail ratio restricted to **rare** vs **common** SNPs (using allele-frequency column).                                                                           | Similar ratios; large gaps suggest frequency-specific artefacts |


**Composite calibration score** (per trait, **lower is better**): a **weighted sum of penalties** measuring distance from the null for tail quantile λ, median λ, tail excess, rare–common discrepancy, and a **capped, down-weighted** GOF term so that tail behaviour drives the score. **Robust *z*-scores** of that composite are formed **within each workflow run** (median and MAD), then mapped to a **run-relative** “quality probability” via a logistic transform; this is a **ranking aid within the batch**, not a calibrated Bayesian posterior.

Full formula detail is documented alongside the **gwas-calibration-utils** package (see **Attribution**).

---

## End-to-end flow

The diagram links **Cromwell orchestration**, **data localisation**, **masking**, **metric and score computation**, and **artefacts**. Solid steps run on the worker VM inside the container; dashed notes are conceptual (null expectations).

```mermaid
flowchart TB
  subgraph wf["Cromwell workflow"]
    IN_JSON[("Inputs JSON<br/>paths, labels, docker URI")]
    BR{two_setups?}
    T1[["Task: single setup"]]
    T2[["Task: two setups"]]
    IN_JSON --> BR
    BR -->|false| T1
    BR -->|true| T2
  end

  subgraph task["Each task (VM + container)"]
    PL[["Read path list(s)<br/>one gs:// URI per line"]]
    LOC[["Cromwell localises objects<br/>to ephemeral disk"]]
    RES[["Resolve trait keys<br/>from filenames"]]
    CIS{{"cis JSON<br/>provided?"}}
    MASK_CIS[["Drop SNPs in cis intervals<br/>or listed positions"]]
    LD{{"lead JSON<br/>provided?"}}
    GEN_LD[["Else target path under results/<br/>engine writes leads from<br/>genome-wide min p per trait"]]
    LOAD_LD[["Load lead chrom / pos"]]
    MASK_LD[["Drop SNPs with<br/>|pos - lead| <= lead_window_bp<br/>(default ±1.5 Mb)"]]
    TAB[["Masked trans table<br/>per trait"]]
    MET[["chi-square transforms,<br/>lambda_GC, lambda_q,<br/>K/E, GOF, rare vs common"]]
    SC[["Composite score<br/>robust z, quality prob"]]
    OUT[["results/ CSV, JSON,<br/>tiered report, logs"]]
    PLOTS{{diagnostic_plots?}}
    FIG[["QQ, tail-excess, λ curves<br/>(optional)"]]
    TGZ[["tar → calibration_outputs.tgz"]]

    PL --> LOC
    LOC --> RES
    RES --> CIS
    CIS -->|yes| MASK_CIS
    CIS -->|no| LD
    MASK_CIS --> LD
    LD -->|no| GEN_LD
    LD -->|yes| LOAD_LD
    GEN_LD --> LOAD_LD
    LOAD_LD --> MASK_LD
    MASK_LD --> TAB
    TAB --> MET
    MET --> SC
    SC --> OUT
    OUT --> PLOTS
    PLOTS -->|yes| FIG
    PLOTS -->|no| TGZ
    FIG --> TGZ
  end

  T1 --> PL
  T2 --> PL

  subgraph null["Interpretation (conceptual)"]
    N1["Under null: p ~ Uniform(0,1)<br/>on masked SNPs"]
    N2["λ ≈ 1, K/E ≈ 1<br/>GOF not systematically significant"]
  end

  MET -.-> N1
  SC -.-> N2

  subgraph deloc["Outputs"]
    GCS[["Declared Files +<br/>optional final_workflow_outputs_dir"]]
  end

  TGZ --> GCS
```



**Lead variants:** If the workflow input **`lead_variants_json`** is **omitted**, the task passes an output path under **`results/`** that **does not exist yet**. The engine **creates** it by scanning the localised sumstats and taking the **genome-wide minimum *p*-value** row per trait, then applies the lead-window mask (**default `lead_window_bp` = 1 500 000** → **±1.5 Mb** on the lead chromosome, **3 Mb** total cis-like span). Override via **`gwas_calibration_qc.lead_window_bp`** in inputs JSON. Supply **`lead_variants_json`** only when you require a **fixed**, pre-computed lead file from object storage.

**Cis regions:** Optional JSON defines per-trait intervals and/or explicit positions to exclude before metrics.

---

## Repository layout

```
GWAS_calibration_engine/
├── LICENSE
├── docker/
│   └── Dockerfile
├── scripts/
│   ├── generate_path_lists.py
│   └── harvest_calibration_outputs.py
├── wdl/
│   ├── gwas_calibration_qc.wdl
│   └── gwas_calibration_qc.example.json
├── requirements.txt
└── README.md
```

Chunk-specific path lists and Cromwell inputs that point at internal buckets are **not** tracked in this repository (see **`.gitignore`**). Build path lists with **`scripts/generate_path_lists.py`** and start from **`wdl/gwas_calibration_qc.example.json`**.

---

## Workflow inputs (summary)


| Input                            | Type            | Description                                                                                               |
| -------------------------------- | --------------- | --------------------------------------------------------------------------------------------------------- |
| `paths_setup_a`                  | `File`          | Text file: one `gs://` URI per line (sumstats per trait).                                                 |
| `two_setups`                     | `Boolean`       | `true` = compare two setups; `false` = single batch.                                                      |
| `paths_setup_b`                  | `File?`         | Second path list; required when `two_setups` is `true`.                                                   |
| `lead_variants_json`             | `File?`         | Optional fixed lead file. Omit to auto-generate under `results/` on the worker.                           |
| `cis_json`                       | `File?`         | Optional cis masking specification.                                                                       |
| `setup_label_a`, `setup_label_b` | `String`        | Human-readable setup names (second label unused when `two_setups` is `false`).                            |
| `protein_id_mode`                | `String`        | How trait IDs are derived from filenames (`first_segment` or `stem`).                                     |
| `n_jobs`, `memory_gb`            | `Int`           | Parallel worker count and worker RAM (gigabytes; formatted as `"{memory_gb} GB"` in the backend runtime). |
| `diagnostic_plots`               | `Boolean`       | Emit per-trait diagnostic figures when dependencies are present in the image.                             |
| `top_n_trans`, `probability_rho` | `Int` / `Float` | Analysis tuning (tail reporting depth; sigmoid sharpness for run-relative quality).                       |
| `lead_window_bp`                 | `Int`           | Half-width in bp around each lead on the lead chromosome (default **1500000** → ±1.5 Mb, **3 Mb** total span). **0** disables lead masking while still loading leads. |
| `docker`                         | `String`        | Full container image URI (including tag).                                                                 |


---

## Parallelism, scattering, and jobs per task

This workflow **does not use Cromwell `scatter`** over traits. Each successful submission runs **exactly one** backend job: either **`run_calibration_one_setup`** or **`run_calibration_two_setups`**. There is **no** fan-out of one virtual machine per protein; Cromwell schedules **one task instance** per workflow run (plus the usual retries).

**What the path list represents:** Every line in **`paths_setup_a`** (and, when comparing, **`paths_setup_b`**) is a **`gs://`** object that Cromwell **localises to the same worker** before the command starts. All listed traits are therefore processed **inside that single task**, sharing CPU, RAM, and ephemeral disk.

**`n_jobs` (shards inside one job):** The input **`n_jobs`** is passed through to the calibration engine as its worker pool size and is also set as **`runtime.cpu: n_jobs`** in the WDL. The engine uses a **process pool** with at most **`n_jobs`** concurrent trait-level workers when **`n_jobs` > 1**; with **`n_jobs` = 1**, traits are processed sequentially. This controls **intra-task** parallelism only. It does **not** change how many Cromwell jobs are created.

**Disk sizing:** Task disk is **`ceil`** of the **total localised input size** of all path-listed files (both setups when comparing) plus a fixed padding (**20 GiB** in the WDL). Larger path lists or bigger sumstats therefore request proportionally larger SSDs on the same single worker.

**Sharding across many Cromwell jobs (optional pattern):** To obtain **true** scatter — many independent jobs, each with its own VM and smaller disk footprint — **split** your master path list into **N** text files (non-overlapping subsets of traits) and **submit N workflow runs** (or wrap **`gwas_calibration_qc`** in a parent WDL that uses Cromwell **`scatter`** over those lists). Each such **shard** is then one Cromwell job with its own **`n_jobs`**, **`memory_gb`**, and localised data subset. Keep **pairing order** consistent across **`paths_setup_a`** and **`paths_setup_b`** when using two setups so protein keys still align within each shard.

---

## Workflow outputs

Cromwell materialises **declared output files** (and may copy them to a bucket if **`final_workflow_outputs_dir`** is set in workflow options). Typical artefacts:


| Artefact                                | Role                                                            |
| --------------------------------------- | --------------------------------------------------------------- |
| `calibration_outputs.tgz`               | Archive of the entire `results/` tree.                          |
| `calibration_compare.metrics.long.csv`  | Long-format metrics per trait (and setup, if comparing).        |
| `calibration_compare.summary.json`      | Run-level summary.                                              |
| `calibration_compare.tiered_report.txt` | Human-readable tiered KPI report.                               |
| `gwas_calibration_qc.log`               | Execution log.                                                  |
| `results/lead_variants.json`            | Present when leads were auto-generated or written to that path. |


With **`diagnostic_plots: true`**, expect additional plot directories inside the tarball (e.g. QQ and tail-excess panels per trait).

---

## Building the container image

The **Dockerfile** expects a **build context** that contains **both**:

- this **`GWAS_calibration_engine`** tree (for workflow-specific requirements), and
- the **gwas-calibration-utils** source tree (sibling path **`Python_scripts/gwas_calibration_utils`** in the upstream layout).

### Option A — Local Docker

From that shared parent directory (call it **`REPO_ROOT`**):

```bash
cd REPO_ROOT

docker build \
  --progress=plain \
  -f Github_clones/GWAS_calibration_engine/docker/Dockerfile \
  -t LOCATION-docker.pkg.dev/PROJECT/REGISTRY/gwas-calibration-qc:TAG \
  .
```

Then push:

```bash
gcloud auth configure-docker LOCATION-docker.pkg.dev
docker push LOCATION-docker.pkg.dev/PROJECT/REGISTRY/gwas-calibration-qc:TAG
```

### Option B — Cloud Build (no local Docker)

When Docker is **not installed** on your machine (common on shared compute VMs), use **Google Cloud Build** to build and push the image remotely. The challenge is that the **`REPO_ROOT`** may be very large (hundreds of GiB of sumstats, analysis outputs, etc.), and `gcloud builds submit` tars the **entire** source directory before uploading — making it impractical to point at the monorepo root.

**Solution — lightweight staging context:** Create a small temporary directory that mirrors **only** the two `COPY` source paths the Dockerfile needs, add a `cloudbuild.yaml`, and submit that. Total upload is typically under **1 MiB**.

```bash
STAGE=$(mktemp -d)

# Mirror only the paths referenced in the Dockerfile
mkdir -p "$STAGE/Github_clones/GWAS_calibration_engine/docker" \
         "$STAGE/Python_scripts"

cp REPO_ROOT/Github_clones/GWAS_calibration_engine/docker/Dockerfile \
   "$STAGE/Github_clones/GWAS_calibration_engine/docker/"
cp REPO_ROOT/Github_clones/GWAS_calibration_engine/requirements.txt \
   "$STAGE/Github_clones/GWAS_calibration_engine/"
cp -r REPO_ROOT/Python_scripts/gwas_calibration_utils \
   "$STAGE/Python_scripts/"

# Write a minimal Cloud Build config
cat > "$STAGE/cloudbuild.yaml" <<'YAML'
steps:
  - name: 'gcr.io/cloud-builders/docker'
    args:
      - 'build'
      - '-f'
      - 'Github_clones/GWAS_calibration_engine/docker/Dockerfile'
      - '-t'
      - 'LOCATION-docker.pkg.dev/PROJECT/REGISTRY/gwas-calibration-qc:TAG'
      - '.'
images:
  - 'LOCATION-docker.pkg.dev/PROJECT/REGISTRY/gwas-calibration-qc:TAG'
YAML

# Submit — use a region-located staging bucket if your org policy restricts resources
gcloud builds submit "$STAGE" \
  --config "$STAGE/cloudbuild.yaml" \
  --project PROJECT \
  --region LOCATION \
  --gcs-source-staging-dir gs://YOUR_BUCKET/tmp/cloud_build_staging \
  --timeout=1200

# Clean up
rm -rf "$STAGE"
```

Replace **`LOCATION`**, **`PROJECT`**, **`REGISTRY`**, **`TAG`**, and **`YOUR_BUCKET`** with your values. The staging tarball is automatically removed from Cloud Build's source bucket after the build; you may also delete the `gs://` staging prefix manually.

**Why this works:** The Dockerfile's two `COPY` instructions reference **relative** paths (`Github_clones/GWAS_calibration_engine/requirements.txt` and `Python_scripts/gwas_calibration_utils`). By reconstructing **only** those subtrees in a temporary directory, we replicate the layout the Dockerfile expects without touching the hundreds of gigabytes of unrelated data in the monorepo.

### Troubleshooting

| Symptom | Likely cause |
|---------|-------------|
| `COPY` fails for `Python_scripts/…` | Build context does not mirror the Dockerfile's expected tree layout. |
| `HTTPError 412: 'us' violates constraint` | Org policy restricts resource locations. Pass **`--region`** and **`--gcs-source-staging-dir`** pointing to a bucket in the allowed region. |
| Image pull fails on Cromwell | Worker service account lacks Artifact Registry read permission, or URI/tag mismatch. |

### Smoke test

After build, the image exposes the **`gwas-calibration-qc`** console entry point:

```bash
docker run --rm LOCATION-docker.pkg.dev/PROJECT/REGISTRY/gwas-calibration-qc:TAG \
  gwas-calibration-qc --help
```

Exit code **0** confirms the environment is wired correctly. If Docker is not available locally, verify via Cloud Build logs or by pulling the image on a different machine.

Use the **same** URI (including tag) as **`gwas_calibration_qc.docker`** in your Cromwell inputs.

---

## Cromwell requirements

The workflow uses **conditional calls** (`if` on `two_setups`) and **`select_first`** over task outputs. Use a **recent Cromwell** (e.g. 50+ or your platform’s supported release) so optional outputs from branches resolve correctly.

Default runtime hints target **preemptible** VMs in **europe-west1** with **no public egress** from the worker (`noAddress: true`), matching common secure batch-QC deployments; adjust zones and flags to match your organisation’s policy.

### Optional `File?` inputs and path localisation

**`cis_json`** and optional **`lead_variants_json`** are **`File?`** workflow inputs. The WDL builds `--cis-json …` and `--lead-variants-json …` flags using **`~{…}`** interpolation **inside** the task `command` block, so Cromwell substitutes **localised disk paths** after download. Do **not** pre-compose those paths as `String` values in the task declaration section by concatenating a `File?` into a string: that can freeze the original **`gs://`** URI in the rendered script and cause the worker to fail when the engine expects a readable local file.

### `lead_window_bp` and the container image

The workflow passes **`--lead-window-bp ~{lead_window_bp}`** explicitly (default **1500000**). That keeps masking behaviour aligned with this repository even if an older image still carried a different Python default. To change the window, set **`gwas_calibration_qc.lead_window_bp`** in Cromwell inputs (use **0** to load leads but disable lead-window masking). Rebuild and push the image when you need an updated **gwas-calibration-utils** implementation inside the container.

### Workflow options and Google labels

Some Cromwell clients require **`google_labels`** (including **`product`**) in the **workflow options** JSON. If your tool treats **`--options`** and **`--l`** as mutually exclusive, add a **`google_labels`** object to the same JSON you pass as **`--options`** (see your platform’s submission examples).

---

## Attribution

- **Numerical methodology and metric definitions:** **gwas-calibration-utils** (Python; MIT-licensed in its source metadata).
- **Path-list and object-storage patterns** follow common practice for large-scale GWAS execution on GCS-backed Cromwell.

---

## Licence

- **This repository** (WDL, Docker build files, scripts under `scripts/`, and this README): **[MIT Licence](LICENSE)** — see **`LICENSE`** for the full text.
- **gwas-calibration-utils** (vendored or installed from source into the image): **MIT** under that package’s own copyright and licence notice; retain those notices when redistributing source or images that bundle the package.