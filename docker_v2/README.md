# MetaSplice 2.0 Docker image

MetaSplice predicts the pathogenicity of intronic single-nucleotide variants (SNVs) outside the canonical splice
sites. This image scores variants from start to finish:

    variants (VCF or table; GRCh38, or hg19 lifted to GRCh38)
      -> scope: intronic positions of GENCODE v50 MANE Select and MANE Plus Clinical transcripts
      -> SpliceTransformer and Pangolin (GPU if available)
      -> the 38 MetaSplice features
      -> MetaSplice score (0-1)

It gives the published model's scores. On 2,200 ClinVar variants scored from scratch on the GPU, every score was
bit-identical to the published score (`tests/`).

## Requirements

- Docker. For the GPU, an NVIDIA GPU with the NVIDIA Container Toolkit (`--gpus all`) and a CUDA 12 driver
  (525.60.13 or newer; tested with 580.173.02 on an RTX 3090). The run used about 1.25 GB of GPU memory and about 4 GB
  of RAM with the default 8 scoring processes. Without a GPU, add `--cpu`: SpliceTransformer and Pangolin then run more
  slowly, and scores can differ slightly from GPU scores (within 1.35e-5 in the CPU acceptance test). The run prints
  the device it uses; without `--gpus all` and without `--cpu` it falls back to the CPU and says so.
- The annotation data bundle `metasplice_v2_data` (219 GB; `bundle/README.md`), mounted at `/data`. It is mirrored on
  USF Box (https://usf.box.com/s/44v8lpk6rxi8d6vdvw1032cfr5olmzdq). Download it, with every file checked, using:

      mkdir -p "$PWD/metasplice_v2_data"
      docker run --rm --user $(id -u):$(id -g) -v "$PWD/metasplice_v2_data":/data --entrypoint python \
          metasplice:2.0.1 /opt/metasplice/bundle/tools/fetch_and_build.py --out /data --all-from-mirror

  Create the folder first, as above: a folder that Docker creates belongs to root, and the download then stops with
  "Permission denied". If an earlier attempt left an empty `metasplice_v2_data` owned by root, remove it
  (`sudo rmdir metasplice_v2_data`) or use another folder.

## Getting the image

Pull it from the GitHub Container Registry, and give it the short name the commands here use:

    docker pull ghcr.io/xiaoming-liu/metasplice:2.0.1
    docker tag ghcr.io/xiaoming-liu/metasplice:2.0.1 metasplice:2.0.1

**Versions.** Image 2.0.1 differs from image 2.0.0 only in documentation: the description of the gene-level values,
which are now built from gnomAD, ClinGen and HGNC (`data/README_gene_level.md`), the notices, and the data bundle's
README and manifest. The scorer and model are version 2.0.0, which a run reports, and both images give identical
scores.

`./build.sh` in this folder built the released image. Rebuilding needs the study's working tree: the pinned
SpliceTransformer and Pangolin repositories and the model files at the paths named in `build.sh`, each checked by md5.

## Running

    docker run --rm --gpus all --user $(id -u):$(id -g) \
        -v /path/to/metasplice_v2_data:/data:ro -v "$PWD":/work \
        metasplice:2.0.1 --input /work/variants.vcf --output /work/scores.tsv

| Option | Meaning |
|---|---|
| `--input` | VCF (`.vcf`, `.vcf.gz`) or table (`.tsv`; otherwise comma-separated) with columns `chr`, `pos`, `ref`, `alt` (or `CHROM`, `POS`, `REF`, `ALT`). Only SNVs are read, in either letter case; other records are skipped, and the run reports how many and why. Split multi-allelic VCF records first (for example `bcftools norm -m -`). An input with no SNVs gives a table with only the header. Tables are read uncompressed, as UTF-8 without a byte-order mark (Excel's "CSV UTF-8" adds one). A table without those four columns, a VCF line with fewer than five columns, a position that is not a whole number or is below 1, a compressed table, a byte-order mark or text that is not UTF-8 stops the run with a message naming the file (and the line, where there is one). |
| `--output` | Output table (tab-separated). |
| `--genome hg19` | The input is GRCh37/hg19. Each REF is checked against GRCh37, then the variant is lifted to GRCh38 with the bundle's Ensembl-built tables; the `status` column gives the reason for any variant that is not lifted. |
| `--cpu` | Run SpliceTransformer and Pangolin on the CPU. |
| `--processes N` | Scoring processes, each taking whole chromosomes (default: the smaller of 8 and the number of CPUs). |
| `--features-out FILE` | Also write the 38 features of each variant. |
| `--include-out-of-scope` | Also list input variants that are out of scope, with the reason. |
| `--keep-work DIR` | Keep the intermediate files (SpliceTransformer and Pangolin outputs, logs). Without it, the temporary work directory is removed when the run ends, whether or not it succeeds. |

`--user $(id -u):$(id -g)` makes the output files yours; the image runs under any user ID. `docker stop` ends a run
at once (`tini` is the container's init). The image writes only to its work directory (`/tmp`, or `--keep-work`), so it
also runs with a read-only root filesystem (`docker run --read-only --tmpfs /tmp`) and under Apptainer/Singularity:

    apptainer build metasplice_2.0.1.sif docker://ghcr.io/xiaoming-liu/metasplice:2.0.1    # no Docker needed
    apptainer run --nv -B /path/to/metasplice_v2_data:/data:ro -B "$PWD":/work metasplice_2.0.1.sif \
        --input /work/variants.vcf --output /work/scores.tsv

Apptainer passes the host's environment and mounts your home directory; the image sets `PYTHONNOUSERSITE=1` and an
empty `PYTHONPATH`, so Python packages installed on the host (`pip install --user`, `PYTHONPATH`) cannot replace the
image's pinned ones. On a machine with Docker and the image, `docker-daemon://metasplice:2.0.1` builds the same SIF.

**Output columns.** `chr`, `pos`, `ref`, `alt`, `region`, `gene`, `gene_id`, `gnomad_af`, `gnomad_popmax_af`,
`status`, `splice_score`. With `--genome hg19`, `chr`, `pos`, `ref` and `alt` are the input's GRCh37 coordinates, and
`chr38`, `pos38`, `ref38` and `alt38` give the GRCh38 variant that was scored.

- `gnomad_af` and `gnomad_popmax_af` are reported for information and are not model inputs. A variant absent from
  gnomAD v4.1 is reported as 0. The columns are empty if the bundle has no `gnomad/` folder.
- `status` says how the variant reached GRCh38: `GRCh38 input`, or lifted from hg19 (forward or inverted strand). An
  hg19 input is checked against GRCh37 first; a variant that is not lifted gives the reason, for example
  `not lifted: REF does not match hg19.fa` or `not lifted: REF differs between builds`.
- `splice_score` is written at full precision, so the exact thresholds below give the model's own assignments. Read
  it exactly as well. pandas' default CSV parser can change the last digits; use
  `pd.read_csv(..., float_precision="round_trip")`.

**Run time.** As a guide, the 2,200 acceptance variants took 38–43 minutes on an NVIDIA RTX 3090 with 16 CPUs:
about 7 minutes for SpliceTransformer and Pangolin, and 32–36 minutes for scoring with 8 processes. Most of the
scoring time is spent reading the DANN and conservation files, so time per variant falls as the number of variants
grows. hg19 input adds a liftover pre-filter, about 15 minutes for 2,200 variants (52–58 minutes in all).

## Scope

A variant is scored if it lies at an intronic position of a GENCODE v50 transcript tagged MANE Select or MANE Plus
Clinical (19,373 transcripts in 19,299 genes). Each variant is assigned a region by its distance from the nearest exon
boundary of its gene:

| Region | Positions |
|---|---|
| Canonical Donor Site | +1 to +2 |
| Canonical Acceptor Site | −1 to −2 |
| Donor Region | +3 to +8 |
| Acceptor Region | −3 to −10 |
| PPT Region | −11 to −25 |
| Branch Point Region | −26 to −44 |
| Deep Intronic | beyond these windows |

A position near more than one boundary takes the first applicable region, in the order canonical acceptor, canonical
donor, acceptor, donor, PPT, branch point, deep intronic. The PPT and branch-point regions are positional windows,
not mapped motifs. Deep-intronic positions inside an exon of any GENCODE v50 basic transcript are out of scope.
Canonical-site variants are scored too; the paper reports results for the five non-canonical regions separately.

## Interpreting the score

The score runs from 0 to 1; higher scores indicate stronger predicted pathogenicity. The model was trained on ClinVar
pathogenic/likely pathogenic and benign/likely benign variants.

The paper derives a **candidate** PP3/BP4 rule from likelihood ratios calibrated against ClinVar. The rule is offered
for independent evaluation and is not a clinical recommendation. It assigns:

- **PP3 at Supporting strength** to scores ≥ 0.5439710286928215 (rounded 0.544);
- **BP4 at Supporting strength** to scores ≤ 0.0017023810017853637 (rounded 0.0017).

These are the derived thresholds for PP3 Strong and BP4 Moderate, applied at Supporting strength. Use the exact
values: in the paper's calibration cohort, the rounded BP4 threshold changes 9 assignments. As
defined, the rule would replace, not add to, a SpliceAI-based PP3 or BP4 for the same variant.

## Exactness

- SpliceTransformer runs with TF32 off (`NVIDIA_TF32_OVERRIDE=0`, set in the image), as in training.
- SpliceTransformer (commit `b67a51d`) and Pangolin 1.0.2 (commit `5cf94b8`) are the versions used in training. Each
  carries one small change made for every Release 2.0 run, shipped as a patch with a modification notice
  (`THIRD_PARTY_NOTICES`).
- Python packages are pinned to the training environment (`requirements.lock`; PyTorch 2.6.0 with CUDA 12.4).
- The data bundle must be the published one. Some public copies differ in ways that change scores: UCSC's 2023
  phyloP 470-way files and DANN rounded to three decimals (`bundle/README.md`).

**Acceptance test.** The test set is 2,200 GRCh38 ClinVar variants (seed 42): 200 per non-canonical region and 50
per canonical site, from each of the training and test sets. It is scored from scratch and compared with the
published scores:

    tests/run_acceptance.sh /path/to/metasplice_v2_data [image] [gpu|hg19|cpu] [results_dir]

Results go to `results_dir` (default `./metasplice_test_results`); the fixtures in `tests/` are only read.

It passes when every variant is scored in the same region, at least 99% of scores agree within 1e-4 and none differs
by more than 0.01.

- `gpu`: all 2,200 scores were bit-identical to the published full-precision scores.
- `hg19`: the same variants in their GRCh37 coordinates (`tests/make_hg19_set.py`). The 2,198 that the Ensembl-built
  tables lift back to the same GRCh38 variant got scores identical to the GRCh38 run; the other 2 have no one-to-one
  position.
- `cpu`: the 57 variants on chr21 and chr22, with SpliceTransformer and Pangolin on the CPU, agreed within 1.35e-5.

**Regression tests.** `tests/run_regressions.sh /path/to/metasplice_v2_data` checks five defects that were fixed
before release:
- an hg19 input whose REF is wrong for GRCh37 but matches GRCh38 (chr19:11221454 C>G) is not lifted or scored;
- scores are exported at full precision. MEN1 chr11:64808119 T>C scores 0.0017024071421474218, 2.6e-8 above the BP4
  threshold; six decimals had put it below;
- a relative `--keep-work` directory works (SpliceTransformer used to look for its input in the wrong directory);
- an input with no SNVs gives an empty table instead of a crash;
- an `hg38.fa.fai` older than `hg38.fa` (as a download or copy can leave it) does not stop SpliceTransformer on the
  read-only `/data`; this run also has a read-only root filesystem (`docker --read-only`), as under Apptainer.

## Licences

MetaSplice's code, model and derived data are under the MetaSplice licence (`/opt/metasplice/LICENSE`): free for
non-commercial use.

SpliceTransformer (Apache-2.0) and Pangolin (GPL-3.0, with its full source in `/opt/Pangolin`) keep their own
licences. So do the Python packages and the annotation data. See `THIRD_PARTY_NOTICES`.

SpliceAI is not included: its licence does not permit redistribution.

## Files

| File | Purpose |
|---|---|
| `Dockerfile`, `build.sh` | Image definition; build script that stages and md5-checks every input |
| `metasplice_pipeline.py` | Entry point (`metasplice`) |
| `requirements.lock` | Pinned Python packages |
| `patches/` | The SpliceTransformer and Pangolin changes, with modification notices |
| `list_licenses.py` | Writes the licence list of the installed Python packages at build time |
| `THIRD_PARTY_NOTICES` | Third-party software and data notices |
| `bundle/` | The data bundle's README, manifest, licences and tools (`fetch_and_build.py`) |
| `tests/` | Acceptance and regression tests |
