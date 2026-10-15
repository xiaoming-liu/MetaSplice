# MetaSplice

MetaSplice predicts the pathogenicity of intronic single-nucleotide variants (SNVs) around and beyond the splice sites:
the donor and acceptor regions, the polypyrimidine tract, the branch-point window and deep-intronic positions. It
combines splicing predictions (SpliceTransformer, Pangolin), conservation, DANN, GPN-MSA and gene-level constraint in
an XGBoost model trained on ClinVar.

This is **MetaSplice 2.0** (scorer and model 2.0.0; Docker image 2.0.1): GRCh38, GENCODE v50 MANE transcripts, 38
features. MetaSplice 1 (GRCh37, 53 features, scorer v1.0.3) remains available at tag
[`v1.0.3`](https://github.com/xiaoming-liu/MetaSplice/tree/v1.0.3).

## Quick start: the Docker image

The image runs everything from a VCF to the score: scope, SpliceTransformer and Pangolin (on the GPU if one is
available), the 38 features, and the model. It gives the published model's scores: on 2,200 ClinVar variants scored
from scratch on the GPU, every score was bit-identical to the published score.

1. Pull the image:

       docker pull ghcr.io/xiaoming-liu/metasplice:2.0.1
       docker tag ghcr.io/xiaoming-liu/metasplice:2.0.1 metasplice:2.0.1

   The release is the image with manifest digest
   `sha256:c32fc8d49da016efd39bc0ac3738ee46d873476e775c0a01a6c714621914f105` (image ID `2dee3c2e…`).

2. Get the annotation data bundle (219 GB; 195 GB without the optional gnomAD tables). It is mirrored on USF Box
   (https://usf.box.com/s/44v8lpk6rxi8d6vdvw1032cfr5olmzdq); the image's tool downloads it and checks every file:

       docker run --rm --user $(id -u):$(id -g) -v "$PWD/metasplice_v2_data":/data --entrypoint python \
           metasplice:2.0.1 /opt/metasplice/bundle/tools/fetch_and_build.py --out /data --all-from-mirror

   The tool can also rebuild the bundle from the original sources; see
   [`docker_v2/bundle/README.md`](docker_v2/bundle/README.md).

3. Score:

       docker run --rm --gpus all --user $(id -u):$(id -g) \
           -v "$PWD/metasplice_v2_data":/data:ro -v "$PWD":/work \
           metasplice:2.0.1 --input /work/variants.vcf --output /work/scores.tsv

   Add `--genome hg19` for GRCh37 input (each REF is checked against GRCh37, then the variant is lifted to GRCh38), or
   `--cpu` without a GPU. On an HPC cluster, Apptainer/Singularity runs the image directly:

       apptainer build metasplice_2.0.1.sif docker://ghcr.io/xiaoming-liu/metasplice:2.0.1
       apptainer run --nv -B "$PWD/metasplice_v2_data":/data:ro -B "$PWD":/work metasplice_2.0.1.sif \
           --input /work/variants.vcf --output /work/scores.tsv

[`docker_v2/README.md`](docker_v2/README.md) describes the options, input handling, output columns, hardware, run
times and the acceptance and regression tests.

## Scope and features

A variant is scored if it lies at an intronic position of a GENCODE v50 transcript tagged MANE Select or MANE Plus
Clinical. Its region follows from its distance to the nearest exon boundary of its gene: canonical donor (+1 to +2) and
acceptor (−1 to −2) sites, donor region (+3 to +8), acceptor region (−3 to −10), polypyrimidine tract (PPT, −11 to
−25), branch-point window (−26 to −44) and deep intronic beyond these. Deep-intronic positions inside an exon of any
GENCODE v50 basic transcript are out of scope. Out-of-scope variants are reported with the reason
(`--include-out-of-scope`).

The 38 features:

| Group | Features |
|---|---|
| SpliceTransformer | the variant score and 15 tissue scores (16) |
| Pangolin | increase, decrease and maximum change in splice-site usage (3) |
| Conservation | phyloP and phastCons, 100-way vertebrate and 470-way mammalian (4) |
| DANN | GRCh37 score at the lifted position (1) |
| GPN-MSA | score (1) |
| Gene level | gnomAD pLI, gnomAD LOEUF, ClinGen haploinsufficiency score (3) |
| Exon context | exon length, exon length mod 3, distance to the boundary (3) |
| Region | one-hot region (7) |

Allele frequency is reported (`gnomad_af`, `gnomad_popmax_af`) but is not a model input.

## Interpreting the score

The score runs from 0 to 1; higher scores indicate stronger predicted pathogenicity. The paper derives a **candidate**
PP3/BP4 rule from likelihood ratios calibrated against ClinVar. The rule is offered for independent evaluation and is
not a clinical recommendation. It assigns:

- **PP3 at Supporting strength** to scores ≥ 0.5439710286928215 (rounded 0.544);
- **BP4 at Supporting strength** to scores ≤ 0.0017023810017853637 (rounded 0.0017).

Use the exact values: in the paper's calibration cohort, the rounded BP4 threshold changes 9 assignments. Scores are
written at full precision; read them exactly as well (for pandas, `pd.read_csv(..., float_precision="round_trip")`).

## The scorer without Docker

`score_v2.py` computes the 38 features and the score from SpliceTransformer and Pangolin outputs and the annotation
bundle. It is the scorer inside the image, and is useful when SpliceTransformer and Pangolin are run separately, for
example across many GPUs.

    pip install -r requirements.txt        # Python 3.12; the versions the model was trained with
    python test_v2.py                      # synthetic tests, no data needed

    D=metasplice_v2_data
    python score_v2.py --input variants.vcf --output scores.tsv \
        --spt-scores spt.csv --pangolin-scores pangolin.csv \
        --fasta $D/grch38/hg38.fa --hg19-fasta $D/grch37/hg19.fa --gtf $D/gencode/gencode.v50.annotation.gtf.gz \
        --reverse-liftover-dir $D/chain/hg38tohg19_ensembl --dann-dir $D/dann \
        --phylop-dir $D/conservation/phyloP_GRCh38 --phastcons-dir $D/conservation/phastCons_GRCh38 \
        --gpn-msa $D/gpn_msa/scores.tsv.bgz --gnomad-dir $D/gnomad \
        --dbnsfp-gene data/metasplice_gene_level_dbNSFP5.4.tsv.gz

Every option can also be set by an environment variable (`python score_v2.py --help`). `--features-out` writes the 38
features; `--features-only` writes them without scoring.

SpliceTransformer and Pangolin must be run exactly as the image runs them, or the scores change:

- SpliceTransformer at commit `b67a51d` with `docker_v2/patches/spt_torch26.patch`, with TF32 off
  (`NVIDIA_TF32_OVERRIDE=0`), and `--reference hg38 --vcf True --raw_score True`, from a run directory holding the
  GENCODE v50 MANE annotation;
- Pangolin 1.0.2 (commit `5cf94b8`) with `docker_v2/patches/pangolin_softmask.patch` and the GENCODE v50 MANE database,
  with `-c CHROM,POS,REF,ALT`.

Both need the GRCh38 variants: for GRCh38 input, run them on all input SNVs; for hg19 input, use the image, which lifts
and selects the variants first. The MANE annotation files are in the image (`/opt/spt_run/data/data_package/`,
`/opt/annotation/gencode.v50.mane.db`), and [`docker_v2/metasplice_pipeline.py`](docker_v2/metasplice_pipeline.py)
shows the exact commands and input formats. The Python environment for both tools is
[`docker_v2/requirements.lock`](docker_v2/requirements.lock) with PyTorch 2.6.0 (CUDA 12.4).

## Changes from MetaSplice 1

- GRCh38 throughout, with GENCODE v50 MANE Select and MANE Plus Clinical transcripts defining regions, genes and exon
  context. DANN, which exists only for GRCh37, is looked up by lifting each position back.
- 38 features instead of 53: no allele frequency, RVIS, Eigen, GERP++, exonic splicing enhancer or silencer motifs, or
  SMS scores.
- Gene-level features are looked up by Ensembl gene ID; a missing LOEUF stays missing.
- SpliceTransformer and Pangolin scores are matched to variants by chromosome, position and alleles.
- Nothing is skipped silently: out-of-scope variants and records that are not SNVs are reported with the reason.

## Files

| Path | Contents |
|---|---|
| `score_v2.py` | The scorer: scope, features, model |
| `test_v2.py` | Synthetic regression tests of the scorer |
| `model/` | The model (`metasplice_v2_full.joblib`) and its feature list (`metasplice_v2_features.joblib`) |
| `data/` | The gene-level values (gnomAD pLI and LOEUF, ClinGen haploinsufficiency score), their README, and `build_gene_level.py`, which builds them from gnomAD, ClinGen and HGNC (`data/sources/`) |
| `docker_v2/` | The Docker image: Dockerfile, pipeline, pinned environment, patches, tests, and the data bundle's README, manifest and tools |
| `requirements.txt` | Python packages for the scorer |
| `LICENSE`, `THIRD_PARTY_NOTICES` | Licences |
| `docs/FEATURE_LICENCES.md` | The licence terms of each input, quoted from its source |

Comments in the code refer to the study's planning documents (`docs/…`) and build scripts (`reproducibility_v2/…`).
The build scripts are in the reproducibility package (`scripts/build/`).

## Reproducibility

The analyses, data and models of the paper, with a script that reruns them, are in
[`xiaoming-liu/MetaSplice-reproducibility`](https://github.com/xiaoming-liu/MetaSplice-reproducibility). Its data
folder (291 MB) is an attachment of that repository's release `v2.0.1`, and its `get_data.py` downloads and checks it.

## Licences

MetaSplice's code, model and derived data are under the MetaSplice licence ([`LICENSE`](LICENSE)). It grants
permission "to use, copy, modify, and distribute MetaSplice and its associated model files for non-commercial
research, educational, and clinical (non-revenue-generating) purposes, subject to the following conditions", and
states: "For commercial use — including use within commercial products, fee-for-service offerings, and internal use by
for-profit entities — please contact xiaomingliu@usf.edu."

SpliceTransformer (Apache-2.0), Pangolin (GPL-3.0) and the annotation data keep their own terms: see
[`THIRD_PARTY_NOTICES`](THIRD_PARTY_NOTICES), and [`docs/FEATURE_LICENCES.md`](docs/FEATURE_LICENCES.md) for each
source's terms, quoted. SpliceAI is not included: its licence does not permit redistribution.

## Citation

A citation will be added here when the paper is published.
