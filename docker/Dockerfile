# Self-contained build: the calibration utils are a git submodule at
# ./gwas_calibration_utils (src-layout). Build from the repository root:
#
#   git submodule update --init --recursive
#   docker build -f docker/Dockerfile \
#     -t LOCATION-docker.pkg.dev/PROJECT/REGISTRY/gwas-calibration-qc:TAG \
#     .
#
# Or, without local Docker, via Cloud Build (see README.md § "Building the
# container image"): the staged directory is the repo root, and gcloud builds
# tars only files tracked by git plus the submodule tree.

FROM python:3.11-slim

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends bash gzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install heavy scientific deps first for better Docker layer cache reuse.
COPY requirements.txt /tmp/calibration_wdl_requirements.txt
RUN pip install --no-cache-dir -r /tmp/calibration_wdl_requirements.txt

# Self-contained: install the utils submodule as a Python package.
# `pip install '.[all]'` pulls matplotlib, gwaslab, intervaltree, tqdm so the
# WDL can run the optional diagnostic plots and the GWASlab mqq step.
COPY gwas_calibration_utils /opt/gwas_calibration_utils
RUN pip install --no-cache-dir '/opt/gwas_calibration_utils[all]'

# Entry points (gwas-calibration-qc, gwaslab-mqq-plots, build-lead-variants-json)
# are installed by the setuptools `[project.scripts]` block and resolve via $PATH.
