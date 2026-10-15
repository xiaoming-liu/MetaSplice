#!/usr/bin/env bash
# Stage the build context for the MetaSplice Release 2.0 image from the exact sources used for the published model,
# check every staged file against its recorded md5, then build the image.
#
#   ./build.sh [tag]        (default tag metasplice:2.0.1)
#
# Sources (all local; nothing is downloaded except Python packages during docker build):
#   SpliceTransformer  git archive of commit b67a51d (unmodified upstream; the patch is applied in the Dockerfile)
#   Pangolin           git archive of commit 5cf94b8 (unmodified upstream; the patch is applied in the Dockerfile)
#   weights, MANE annotation databases, model and gene-level extract: md5-checked copies
set -euo pipefail
cd "$(dirname "$0")"
P=/mnt/Work0/Splice_altering_variant
TAG=${1:-metasplice:2.0.1}
C=context
rm -rf "$C" && mkdir -p "$C"/{spt,pangolin,annotation,spt_run/data/data_package,spt_run/pyensembl_cache,metasplice/model,metasplice/data}

check() {  # check <file> <md5>
  local got; got=$(md5sum "$1" | cut -d' ' -f1)
  [[ "$got" == "$2" ]] || { echo "md5 mismatch: $1 is $got, expected $2" >&2; exit 1; }
  echo "  ok  $2  $1"
}

echo "Upstream sources (git archive of the pinned commits)"
# Each tool's archive and patch in its own folder: the Dockerfile bind-mounts only that folder, and BuildKit keys the
# step on the mounted content, so a change elsewhere in the context does not rebuild the tool's layer.
git -C "$P/SpliceTransformer" archive --format=tar b67a51d > "$C/spt/SpliceTransformer-b67a51d.tar"
git -C "$P/Pangolin" archive --format=tar 5cf94b8db938c658391b4305cd7ce33297d44ff7 > "$C/pangolin/Pangolin-5cf94b8.tar"
cp patches/spt_torch26.patch "$C/spt/"
cp patches/pangolin_softmask.patch "$C/pangolin/"
cp requirements.lock "$C/requirements.lock"

echo "Weights, annotation, model and gene-level extract"
cp "$P/SpliceTransformer/model/weights/SpTransformer_pytorch.ckpt" "$C/"
check "$C/SpTransformer_pytorch.ckpt" f9f9c491cae4e854eef644a7a97753de
cp "$P/annotation_v2/gencode.v50.mane.db" "$C/annotation/"
check "$C/annotation/gencode.v50.mane.db" 8796342f4a6d5e88506be6b051cc6467                      # V01: Pangolin MANE database
cp -p "$P/annotation_v2/spt_run/data/data_package/hg38.annotation.gtf.db" "$C/spt_run/data/data_package/"
cp -pL "$P/annotation_v2/spt_run/data/data_package/hg38.annotation.gtf.gz" "$C/spt_run/data/data_package/"
check "$C/spt_run/data/data_package/hg38.annotation.gtf.gz" 0fed8b55b40af40e25c8fdf184b6931d   # V01: MANE GTF
check "$C/spt_run/data/data_package/hg38.annotation.gtf.db" da949f5ea3b9ee57940aa8d650684e3c                         # V01b: SpliceTransformer database
cp "$P/MetaSplice/score_v2.py" "$C/metasplice/score_v2.py"
cp metasplice_pipeline.py "$C/metasplice/"
cp list_licenses.py "$C/metasplice/"
cp "$P/reproducibility_v2/models/metasplice_v2_full.joblib" "$P/reproducibility_v2/models/metasplice_v2_features.joblib" "$C/metasplice/model/"
check "$C/metasplice/model/metasplice_v2_full.joblib" 51c6e5d0db100b33b91907fa857cd9a7
check "$C/metasplice/model/metasplice_v2_features.joblib" 8be31360848f2cb43c2b02850bd03815
cp "$P/MetaSplice/data/metasplice_gene_level_dbNSFP5.4.tsv.gz" "$P/MetaSplice/data/README_gene_level.md" "$C/metasplice/data/"
check "$C/metasplice/data/metasplice_gene_level_dbNSFP5.4.tsv.gz" 088ba4d5101c01a8697fbbf503affd06
cp patches/*.patch "$C/metasplice/"
cp "$P/MetaSplice/LICENSE" "$C/metasplice/LICENSE"
[[ -f THIRD_PARTY_NOTICES ]] && cp THIRD_PARTY_NOTICES "$C/metasplice/" || cp "$P/MetaSplice/THIRD_PARTY_NOTICES" "$C/metasplice/"
[[ -f README.md ]] && cp README.md "$C/metasplice/" || true

echo "Data bundle tools and documents (bundle/, with its MANIFEST.tsv from tools/assemble_bundle.py)"
[[ -f bundle/MANIFEST.tsv ]] || { echo "bundle/MANIFEST.tsv missing: run bundle/tools/assemble_bundle.py, then copy its MANIFEST.tsv" >&2; exit 1; }
cp -rL bundle "$C/metasplice/bundle"
find "$C/metasplice/bundle" -name __pycache__ -prune -exec rm -rf {} +
awk -F'\t' -v d="$C/metasplice/bundle" 'NR > 1 && $5 == "local" {print $3 "  " d "/" $1}' bundle/MANIFEST.tsv | md5sum -c --quiet
echo "  ok  $(awk -F'\t' 'NR > 1 && $5 == "local"' bundle/MANIFEST.tsv | wc -l) local bundle files match MANIFEST.tsv"

echo "Code pins (update deliberately when a file changes)"
check "$C/metasplice/score_v2.py" 4a42ae0bb345ff20165bfa10b3cadb3c
check "$C/metasplice/metasplice_pipeline.py" 25bb866fad360e1363d4a8ca727fb258
check "$C/metasplice/list_licenses.py" 79ca55591b92ad5541966948c6321c18
check "$C/spt/spt_torch26.patch" 2090cea6338556fc88485969b19c5584
check "$C/pangolin/pangolin_softmask.patch" f2317c6ccee1050fc8c70f27f76f7d2c
check "$C/requirements.lock" 980adb4f9fb25c079e6d8166bf969700

# The image records the git commit; "-dirty" when the staged sources differ from it.
REV=$(git -C "$P" rev-parse --short HEAD)
[[ -n "$(git -C "$P" status --porcelain -- MetaSplice/score_v2.py MetaSplice/docker_v2 MetaSplice/data MetaSplice/LICENSE)" ]] && REV="$REV-dirty"
echo "Building $TAG (revision $REV)"
docker build --build-arg GIT_COMMIT="$REV" -t "$TAG" .
docker image inspect "$TAG" --format 'image {{.Id}}  size {{.Size}} bytes'
