#!/usr/bin/env bash
# Regression tests of the MetaSplice Release 2.0 image for defects fixed before release (check_regressions.py):
#   1. hg19 input REF checked against GRCh37 (verification of 2026-10-11);
#   2. full-precision scores (verification of 2026-10-11); this run also keeps its work directory under a relative
#      path, which SpliceTransformer used to resolve from the wrong directory (review of 2026-10-14);
#   3. an input with no SNVs gives an empty table, not a crash (review of 2026-10-13);
#   4. an hg38.fa.fai older than its FASTA, on the read-only /data, does not stop SpliceTransformer (review of 2026-10-13),
#      and the image runs with a read-only root filesystem, as under Apptainer/Singularity (re-review of 2026-10-13).
# A few minutes; GPU by default, --cpu for the CPU.
#
#   ./run_regressions.sh /path/to/metasplice_v2_data [image] [--cpu] [results_dir]
set -euo pipefail
DATA=${1:?usage: run_regressions.sh DATA_DIR [image] [--cpu] [results_dir]}
IMAGE=${2:-metasplice:2.0.1}
GPU=(--gpus all); OPTS=(); FLAG=--gpu
[[ "${3:-}" == "--cpu" ]] && { GPU=(); OPTS=(--cpu); FLAG=""; }
OUT=$(mkdir -p "${4:-metasplice_test_results}" && cd "${4:-metasplice_test_results}" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
U=(--user "$(id -u):$(id -g)" -v "$HERE":/tests:ro -v "$OUT":/work)
run() { docker run --rm "${GPU[@]}" "${U[@]}" -v "$DATA":/data:ro "$@" "${OPTS[@]}"; }
run "$IMAGE" --input /tests/regression_hg19.vcf --genome hg19 --output /work/regression_out_hg19.tsv --include-out-of-scope
rm -rf "$OUT/regression_work_relative"
run "$IMAGE" --input /tests/regression_grch38.vcf --output /work/regression_out_grch38.tsv \
    --features-out /work/regression_features_grch38.tsv --keep-work regression_work_relative
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > "$OUT/regression_empty.vcf"
run "$IMAGE" --input /work/regression_empty.vcf --output /work/regression_out_empty.tsv
cp "$DATA/grch38/hg38.fa.fai" "$OUT/stale_hg38.fa.fai" && touch -d 2000-01-01 "$OUT/stale_hg38.fa.fai"
run --read-only --tmpfs /tmp -v "$OUT/stale_hg38.fa.fai":/data/grch38/hg38.fa.fai:ro "$IMAGE" \
    --input /tests/regression_grch38.vcf --output /work/regression_out_stale_index.tsv
docker run --rm "${U[@]}" --entrypoint python "$IMAGE" /tests/check_regressions.py /work $FLAG
