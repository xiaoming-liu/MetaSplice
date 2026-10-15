#!/usr/bin/env bash
# Acceptance test of the MetaSplice Release 2.0 image: score ClinVar variants from scratch and compare them with the
# published scores (compare.py).
#
#   ./run_acceptance.sh /path/to/metasplice_v2_data [image] [gpu|hg19|cpu] [results_dir]
#
#   gpu   the 2,200 GRCh38 acceptance variants on the GPU (default)
#   hg19  the same variants in their GRCh37 coordinates (--genome hg19; make_hg19_set.py), on the GPU
#   cpu   the acceptance variants on chr21 and chr22, with SpliceTransformer and Pangolin on the CPU (--cpu)
# Results go to results_dir (default ./metasplice_test_results); the fixtures in this folder are only read.
set -euo pipefail
DATA=${1:?usage: run_acceptance.sh DATA_DIR [image] [gpu|hg19|cpu] [results_dir]}
IMAGE=${2:-metasplice:2.0.1}
MODE=${3:-gpu}
OUT=$(mkdir -p "${4:-metasplice_test_results}" && cd "${4:-metasplice_test_results}" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
GPU=(--gpus all); OPTS=()
case "$MODE" in
  gpu)  IN=/tests/acceptance_variants.vcf; EXP=/tests/acceptance_expected.tsv ;;
  hg19) IN=/tests/acceptance_variants_hg19.vcf; EXP=/tests/acceptance_expected_hg19.tsv; OPTS=(--genome hg19) ;;
  cpu)  IN=/work/acceptance_variants_cpu.vcf; EXP=/work/acceptance_expected_cpu.tsv; OPTS=(--cpu); GPU=()
        awk -F'\t' '/^#/ || $1 == "chr21" || $1 == "chr22"' "$HERE/acceptance_variants.vcf" > "$OUT/acceptance_variants_cpu.vcf"
        awk -F'\t' 'NR == 1 || $1 == "chr21" || $1 == "chr22"' "$HERE/acceptance_expected.tsv" > "$OUT/acceptance_expected_cpu.tsv" ;;
  *) echo "mode must be gpu, hg19 or cpu" >&2; exit 2 ;;
esac
U=(--user "$(id -u):$(id -g)" -v "$HERE":/tests:ro -v "$OUT":/work)
time docker run --rm "${GPU[@]}" "${U[@]}" -v "$DATA":/data:ro "$IMAGE" \
  --input "$IN" --output "/work/acceptance_out_$MODE.tsv" --keep-work "/work/acceptance_work_$MODE" "${OPTS[@]}"
docker run --rm "${U[@]}" --entrypoint python "$IMAGE" /tests/compare.py "$EXP" "/work/acceptance_out_$MODE.tsv"
