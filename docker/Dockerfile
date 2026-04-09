# Build from the long_gwas repository root so the calibration utils are copied unchanged:
#   docker build -f Github_clones/GWAS_calibration_engine/docker/Dockerfile \
#     -t YOUR_REGISTRY/gwas-calibration-qc:v1 \
#     /path/to/long_gwas
#
# Do not edit files under Python_scripts/gwas_calibration_utils here; the COPY is read-only at build time.

FROM python:3.11-slim

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends bash gzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY Github_clones/GWAS_calibration_engine/requirements.txt /tmp/calibration_wdl_requirements.txt
RUN pip install --no-cache-dir -r /tmp/calibration_wdl_requirements.txt

COPY Python_scripts/gwas_calibration_utils /opt/gwas_calibration_utils
RUN pip install --no-cache-dir -e /opt/gwas_calibration_utils

ENV PATH="/opt/gwas_calibration_utils:${PATH}"
