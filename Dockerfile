# MetaSplice: Splice Variant Pathogenicity Predictor
#
# Build:  docker build -t metasplice .
# Run:    docker run --rm -v /path/to/data:/data metasplice --help
#
# Scoring (default entrypoint):
#   docker run --rm -v /path/to/data:/data metasplice \
#     --input /data/variants.vcf --output /data/scores.tsv ...
#
# One-time hg38 reverse-liftover precompute (alternate entrypoint):
#   docker run --rm -v /path/to/data:/data --entrypoint python metasplice \
#     /app/invert_liftover.py \
#     --input_dir  /data/hg19tohg38 \
#     --output_dir /data/hg38tohg19
#
# Two-stage build to minimize image size

FROM python:3.12-slim AS base

LABEL maintainer="[Author]"
LABEL description="MetaSplice: an ensemble predictor for splice variant pathogenicity"
LABEL version="1.0.1"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    tabix \
    libcurl4-openssl-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir \
    numpy==1.26.4 \
    pandas==2.2.2 \
    scikit-learn==1.5.1 \
    xgboost==2.0.3 \
    joblib==1.4.2 \
    pysam==0.22.1

# Copy application
COPY model/ /app/model/
COPY score.py /app/score.py
COPY invert_liftover.py /app/invert_liftover.py
COPY requirements.txt /app/requirements.txt

WORKDIR /app

ENTRYPOINT ["python", "/app/score.py"]
CMD ["--help"]
