# MetaSplice Release 2.0: end-to-end scoring image (GRCh38; hg19 input is lifted).
#
#   variants -> SpliceTransformer + Pangolin (GPU if available, TF32 off) -> 38 features -> MetaSplice score
#
# Build:  ./build.sh            (stages the build context with checksums, then runs docker build)
# Run:    docker run --rm --gpus all --user $(id -u):$(id -g) -v /path/to/metasplice_v2_data:/data:ro -v "$PWD":/work \
#             metasplice:2.0.1 --input /work/variants.vcf --output /work/scores.tsv
#
# The annotation data bundle (219 GB; bundle/README.md) is mounted at /data. SpliceAI is not included (its licence does not
# permit redistribution). Third-party licences and modification notices: /opt/metasplice/THIRD_PARTY_NOTICES.

# python:3.12.12-slim-bookworm, pinned by digest (the same base layers as every tested build)
FROM python:3.12.12-slim-bookworm@sha256:593bd06efe90efa80dc4eee3948be7c0fde4134606dd40d8dd8dbcade98e669c

LABEL org.opencontainers.image.title="MetaSplice" \
      org.opencontainers.image.version="2.0.1" \
      org.opencontainers.image.description="Pathogenicity of intronic SNVs (GRCh38), with SpliceTransformer and Pangolin" \
      org.opencontainers.image.source="https://github.com/xiaoming-liu/MetaSplice" \
      org.opencontainers.image.licenses="LicenseRef-MetaSplice-NonCommercial"

ENV PIP_NO_CACHE_DIR=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 NVIDIA_TF32_OVERRIDE=0 \
    METASPLICE_APP=/opt/metasplice METASPLICE_DATA=/data XDG_CACHE_HOME=/tmp/.cache TINI_SUBREAPER=1

# tini runs as PID 1, so that docker stop reaches the pipeline at once (TINI_SUBREAPER: no warning where tini is not
# PID 1, as under Apptainer/Singularity)
RUN apt-get update && apt-get install -y --no-install-recommends patch tini && rm -rf /var/lib/apt/lists/*

# PyTorch 2.6.0 with CUDA 12.4 (as in training), then every other package pinned to the training environment.
RUN pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
COPY context/requirements.lock /opt/metasplice/requirements.lock
RUN pip install -r /opt/metasplice/requirements.lock

# SpliceTransformer (Apache-2.0), commit b67a51d, with the PyTorch 2.6 loading patch, and its published weights.
# (each tool's folder is bind-mounted for the build, so the archive does not stay in a layer, and the step is cached
# unless that folder changes)
RUN --mount=type=bind,source=context/spt,target=/ctx \
    mkdir -p /opt/SpliceTransformer && tar -xf /ctx/SpliceTransformer-b67a51d.tar -C /opt/SpliceTransformer \
 && sed '/^#/d' /ctx/spt_torch26.patch | patch -d /opt/SpliceTransformer -p1 \
 && cp /ctx/spt_torch26.patch /opt/SpliceTransformer/MODIFICATIONS.patch
COPY context/SpTransformer_pytorch.ckpt /opt/SpliceTransformer/model/weights/SpTransformer_pytorch.ckpt

# Pangolin 1.0.2 (GPL-3.0), commit 5cf94b8, with the soft-masked reference patch. Its source (including the change)
# stays in /opt/Pangolin.
RUN --mount=type=bind,source=context/pangolin,target=/ctx \
    mkdir -p /opt/Pangolin && tar -xf /ctx/Pangolin-5cf94b8.tar -C /opt/Pangolin \
 && sed '/^#/d' /ctx/pangolin_softmask.patch | patch -d /opt/Pangolin -p1 \
 && cp /ctx/pangolin_softmask.patch /opt/Pangolin/MODIFICATIONS.patch \
 && pip install --no-deps /opt/Pangolin

# GENCODE v50 MANE annotation: Pangolin's database, and SpliceTransformer's annotation files. The pipeline links these
# into a run directory inside its work directory, so the image writes nothing outside it (read-only root filesystems:
# docker --read-only, Apptainer/Singularity).
COPY context/annotation/ /opt/annotation/
COPY context/spt_run/ /opt/spt_run/

# MetaSplice: scorer, pipeline, model, gene-level extract, licences.
COPY context/metasplice/ /opt/metasplice/
RUN ln -s /opt/metasplice/metasplice_pipeline.py /usr/local/bin/metasplice && chmod +x /opt/metasplice/*.py \
 && python /opt/metasplice/list_licenses.py

WORKDIR /work
ENTRYPOINT ["/usr/bin/tini", "--", "python", "/opt/metasplice/metasplice_pipeline.py"]
CMD ["--help"]

# Under Apptainer/Singularity, which pass the host environment and mount $HOME, ignore the host's Python user site and
# PYTHONPATH, so that the pinned packages are used (release review of 2026-10-14). At the end: no layer changes.
ENV PYTHONNOUSERSITE=1 PYTHONPATH=

# The git commit of the sources, last, so that it changes no cached layer: a label-only rebuild has the tested layers.
ARG GIT_COMMIT=unknown
LABEL org.opencontainers.image.revision="${GIT_COMMIT}"
