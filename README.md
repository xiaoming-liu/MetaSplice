# MetaSplice

**An ensemble predictor for splice variant pathogenicity prediction**

MetaSplice integrates evolutionary constraint, deep learning splice predictions,
and exonic splice regulatory motifs to predict pathogenicity of intronic and
near-splice variants. It achieves auPRC 0.9952 on independent ClinVar test data,
outperforming SpliceAI, CADD, and other tools in non-canonical splice regions.

---

## Installation

```bash
git clone https://github.com/xiaoming-liu/MetaSplice.git
cd MetaSplice
```

### Docker (recommended)

```bash
docker build -t metasplice .
```

### Local (no Docker)

```bash
pip install -r requirements.txt
python score.py --help
# For hg38 input mode, also pre-compute the reverse liftover once:
python invert_liftover.py --help
```

### Package contents

```
MetaSplice/
├── Dockerfile            # multi-stage build
├── requirements.txt      # pinned Python deps
├── score.py              # main scorer (default ENTRYPOINT)
├── invert_liftover.py    # one-time hg38->hg19 precompute (alt entrypoint)
├── prepare_pangolin_input.py  # VCF-to-CSV converter for Pangolin
├── test.sh               # 3-variant smoke test
├── README.md             # this file
└── model/
    ├── metasplice.ubj      # trained XGBoost (53 features, native format)
    ├── features.json       # feature column definitions
    └── sms_7mer_scores.json
```

---

## Feature Tiers

| Tier | Features | Requires | GPU? |
|------|----------|----------|------|
| Minimal | 25 | hg19 FASTA + GTF + dbNSFP gene file | No |
| Standard | 44 | Minimal + conservation + DANN + Eigen + gnomAD + GPN-MSA + liftover | No |
| Full | 53 | Standard + SpliceTransformer + Pangolin | Yes |

The XGBoost model handles missing features natively; performance degrades
gracefully as tiers drop features.

---

## Expected Data Layout

For Full-tier scoring, organize your data directory as below. Files for tiers
you skip can be omitted.

```
/path/to/data/
├── variants.vcf                              # your input (hg19 or hg38)
│
├── hg19.fa                                   # ── hg19 input path ──
├── hg19.fa.fai
├── gencode.v19.annotation.gtf.gz
│
├── hg38.fa                                   # ── hg38 input path ──
├── hg38.fa.fai
├── gencode.v49.annotation.gtf.gz
│
├── dbNSFP5.3_gene.gz                         # gene-level (build-agnostic)
│
├── GERP_hg19/                                # hg19-direct databases
│   └── chr{1..22,X,Y}.maf.rates.gz           #   (used by both paths;
├── DANN/                                     #    hg38 input is mapped
│   └── DANN_whole_genome_SNVs.tsv.gz.chr{1..22,X,Y}.gz   # back to hg19)
├── Eigen/
│   └── Eigen_hg19_combined.tab.chr{1..22,X,Y}.gz
│
├── phyloP_GRCh38/                            # hg38-direct databases
│   ├── chr{1..22,X,Y}.phyloP100way.wigFix.gz #   (used by both paths;
│   └── chr{1..22,X,Y}.phyloP470way.wigFix.gz #    hg19 input is mapped
├── phastCons_GRCh38/                         #    forward to hg38)
│   ├── chr{1..22,X,Y}.phastCons100way.wigFix.gz
│   └── chr{1..22,X,Y}.phastCons470way.wigFix.gz
├── gnomAD/
│   └── gnomad.joint.v4.1.sites.chr{1..22,X,Y}.vcf.bgz.snp.tsv.gz
├── gpn_msa/
│   ├── scores.tsv.bgz
│   └── scores.tsv.bgz.tbi
│
├── hg19tohg38/                               # forward liftover (--genome hg19)
│   └── chr{1..22,X,Y}.hg19tohg38.gz
├── hg38tohg19/                               # reverse liftover (--genome hg38;
│   └── chr{1..22,X,Y}.hg38tohg19.gz          #   one-time precompute, see below)
│
├── spt_scores.tsv                            # ── Full tier (GPU) ──
└── pangolin_scores.csv
```

### Downloading the data

A pre-built data package containing all databases for both hg19 and hg38
input modes (including the pre-computed hg38→hg19 reverse liftover) is
available at:

> **https://usf.box.com/s/u583ubbntyjjim7i752lor2zf0e4gjoa**

Download and extract into a single directory (e.g., `/path/to/data/`).
The archive follows the layout shown above — no reorganization needed.

---

## Quick Start

### Minimal scoring (no GPU, no large databases)

```bash
docker run --rm \
  -v /path/to/data:/data \
  metasplice \
  --input /data/variants.vcf \
  --output /data/scores.tsv \
  --fasta /data/hg19.fa \
  --gtf /data/gencode.v19.annotation.gtf.gz \
  --dbnsfp-gene /data/dbNSFP5.3_gene.gz
```

### Standard scoring (no GPU)

```bash
docker run --rm \
  -v /path/to/data:/data \
  metasplice \
  --input /data/variants.vcf \
  --output /data/scores.tsv \
  --fasta /data/hg19.fa \
  --gtf /data/gencode.v19.annotation.gtf.gz \
  --dbnsfp-gene /data/dbNSFP5.3_gene.gz \
  --gerp-dir /data/GERP_hg19 \
  --dann-dir /data/DANN \
  --eigen-dir /data/Eigen \
  --liftover-dir /data/hg19tohg38 \
  --phylop-dir /data/phyloP_GRCh38 \
  --phastcons-dir /data/phastCons_GRCh38 \
  --gnomad-dir /data/gnomAD \
  --gpn-msa /data/gpn_msa/scores.tsv.bgz
```

### Full scoring (with SpliceTransformer + Pangolin)

Score SpT and Pangolin externally on GPU (see next section), then pass their
output files to MetaSplice:

```bash
docker run --rm \
  -v /path/to/data:/data \
  metasplice \
  --input /data/variants.vcf \
  --output /data/scores.tsv \
  --fasta /data/hg19.fa \
  --gtf /data/gencode.v19.annotation.gtf.gz \
  --dbnsfp-gene /data/dbNSFP5.3_gene.gz \
  --gerp-dir /data/GERP_hg19 \
  --dann-dir /data/DANN \
  --eigen-dir /data/Eigen \
  --liftover-dir /data/hg19tohg38 \
  --phylop-dir /data/phyloP_GRCh38 \
  --phastcons-dir /data/phastCons_GRCh38 \
  --gnomad-dir /data/gnomAD \
  --gpn-msa /data/gpn_msa/scores.tsv.bgz \
  --spt-scores /data/spt_scores.tsv \
  --pangolin-scores /data/pangolin_scores.csv
```

**IMPORTANT**: `spt_scores.tsv` and `pangolin_scores.csv` MUST be scored on the
SAME input VCF in the SAME variant order as `variants.vcf`. The Docker entry
point validates row alignment by chr:pos:ref:alt and aborts on mismatch.

---

## hg38 Input Mode

By default, MetaSplice expects hg19 input. To accept hg38 VCFs natively, pass
`--genome hg38` and substitute the hg38-compatible reference files. The model
was trained on hg19 features; on hg38 input, the pipeline routes hg38-bridged
features (phyloP, phastCons, gnomAD, GPN-MSA) through the hg38 databases
directly and routes hg19-direct features (GERP, DANN, Eigen) through an
hg38→hg19 reverse liftover.

### hg38→hg19 reverse liftover

The data package (see *Downloading the data* above) includes pre-computed
`hg38tohg19/` reverse liftover files. No additional setup is needed.

If you prefer to regenerate them from the `hg19tohg38/` files, the Docker
image ships `/app/invert_liftover.py`:

```bash
docker run --rm \
  -v /path/to/data:/data \
  --entrypoint python \
  metasplice \
  /app/invert_liftover.py \
  --input_dir  /data/hg19tohg38 \
  --output_dir /data/hg38tohg19
```

Runtime ~2 hours single-threaded, output ~13 GB.

### Run on hg38 VCF

```bash
docker run --rm \
  -v /path/to/data:/data \
  metasplice \
  --input /data/variants.hg38.vcf \
  --output /data/scores.tsv \
  --genome hg38 \
  --fasta /data/hg38.fa \
  --gtf /data/gencode.v49.annotation.gtf.gz \
  --dbnsfp-gene /data/dbNSFP5.3_gene.gz \
  --gerp-dir /data/GERP_hg19 \
  --dann-dir /data/DANN \
  --eigen-dir /data/Eigen \
  --reverse-liftover-dir /data/hg38tohg19 \
  --phylop-dir /data/phyloP_GRCh38 \
  --phastcons-dir /data/phastCons_GRCh38 \
  --gnomad-dir /data/gnomAD \
  --gpn-msa /data/gpn_msa/scores.tsv.bgz
```

Note that `--liftover-dir` is replaced with `--reverse-liftover-dir`, and FASTA
+ GTF are now hg38 + v49. GERP/DANN/Eigen directories themselves remain hg19
(those databases were not re-released for hg38); MetaSplice does the hg38→hg19
mapping internally.

### Caveats for the hg38 path

The hg38 path is best-effort and does NOT byte-reproduce the published
hg19-path predictions. On a 500-variant chr22 spot check:

- 94% of variants agree within ±0.05 score (471/500)
- 85% agree within ±0.01 (427/500)
- A small minority (29/500, ~6%) diverge by >0.05, up to 0.5 in extreme cases

The main sources of divergence:

1. **Gene renames v19→v49.** dbNSFP5.3 indexes by current HGNC names, so when
   v49 renames a gene (e.g. CECR1→ADA2), the hg38 path finds dbNSFP entries
   that the hg19 path misses. The hg38 path is arguably more accurate here,
   but predictions differ from the trained model's expectations for renamed
   genes (the model was trained with NaN for those variants and may behave
   non-monotonically when those features are now populated).
2. **Region reclassification.** A small fraction of variants near revised v49
   exon boundaries get a different splice region than v19, changing the region
   one-hot feature.
3. **Reverse liftover coverage.** Not every hg38 position successfully maps
   back to hg19; affected variants get NaN for DANN/Eigen/GERP. Typical loss
   is <0.1% in well-annotated regions.
4. **FASTA build.** Reference nucleotides at the same biological position are
   identical between hg19 and hg38, so ESE/ESS hexamer and SMS 7-mer features
   give the same answer with either FASTA (apart from rare indel-affected loci,
   which fall outside our SNV-only scope).

**For byte-identical reproduction of the published model evaluation, use the
default hg19 path.** Use the hg38 path when your VCF is already in hg38 and a
small score difference is acceptable; consider re-validating thresholds on a
held-out hg38 set if you make clinical decisions from raw scores on this path.

---

## Running SpliceTransformer (GPU)

SpliceTransformer requires PyTorch + CUDA and is not bundled in this Docker
image (it would add ~3 GB and require nvidia-docker). Score it in a separate
GPU environment:

```bash
# 1. Install SpliceTransformer
git clone https://github.com/ShenLab-Genomics/SpliceTransformer.git
cd SpliceTransformer
conda create -n spt python=3.10 -y
conda activate spt
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
pip install sinkhorn-transformer pyvcf3 pyensembl pyfaidx

# 2. Set up reference files (SpliceTransformer expects data/data_package/)
mkdir -p data/data_package
cp /path/to/data/hg19.fa data/data_package/hg19.fa
cp /path/to/data/hg19.fa.fai data/data_package/hg19.fa.fai
cp /path/to/data/gencode.v19.annotation.gtf.gz \
   data/data_package/hg19.annotation.gtf.gz

# 3. Build pyensembl database from the GTF (required for gene lookups)
pyensembl install --reference-name hg19 --annotation-name gencode.v19 \
   --gtf "data/data_package/hg19.annotation.gtf.gz"

# 4. Download pre-trained weights (see SpliceTransformer README for URL)
#    Place at: model/weights/SpTransformer_pytorch.ckpt

# 5. Fix PyTorch >=2.0 compatibility: the checkpoint uses underscore-separated
#    ParameterList keys (weights_0) but newer PyTorch expects dot-separated
#    (weights.0). Patch tasks_annotate_mutations.py around line 413:
#
#      save_dict = torch.load(weight_path, map_location="cpu")
#      # --- add these lines ---
#      sd = save_dict["state_dict"]
#      sd = {k.replace("weights_0", "weights.0").replace("weights_1", "weights.1"): v
#            for k, v in sd.items()}
#      save_dict["state_dict"] = sd
#      # --- end patch ---
#      model.load_state_dict(save_dict["state_dict"])

# 6. Score your variants (~9 variants/sec on a single GPU)
python sptransformer.py \
  --input /path/to/data/variants.vcf \
  --output /path/to/data/spt_scores.tsv \
  --reference hg19 \
  --vcf True \
  --raw_score True
```

Expected output columns: `index, #CHROM, POS, ID, REF, ALT, score,
Adipose Tissue, Blood, Blood Vessel, Brain, Colon, Heart, Kidney, Liver,
Lung, Muscle, Nerve, Small Intestine, Skin, Spleen, Stomach`.

---

## Running Pangolin (GPU)

```bash
# 1. Install Pangolin
git clone https://github.com/tkzeng/Pangolin.git
cd Pangolin
conda create -n pangolin python=3.10 -y
conda activate pangolin
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
pip install "setuptools<81" pyvcf3 gffutils pandas pyfastx
pip install .

# 2. Download the GENCODE annotation database (hg19/GRCh37, ~400 MB)
wget -O gencode.v38lift37.annotation.db \
  https://www.dropbox.com/sh/6zo0aegoalvgd9f/AAA9Q90Pi1UqSzX99R_NM803a/gencode.v38lift37.annotation.db

# 3. Patch for soft-masked genomes (required for hg19.fa):
#    In pangolin/pangolin.py, find the line:
#        if seq[5000+d:5000+d+len(ref)] != ref:
#    and add this line BEFORE it:
#        seq = seq.upper()
#    Then reinstall: pip install .
#    Without this patch, variants in soft-masked (repeat) regions are
#    skipped due to case-sensitive ref mismatch (e.g. 'g' != 'G').

# 4. Convert your VCF to CSV (Pangolin outputs CSV when given CSV input,
#    which score.py can read directly)
python /path/to/metasplice/prepare_pangolin_input.py \
  /path/to/data/variants.vcf \
  /path/to/data/pangolin_input.csv

# 5. Score your variants
pangolin \
  /path/to/data/pangolin_input.csv \
  /path/to/data/hg19.fa \
  gencode.v38lift37.annotation.db \
  /path/to/data/pangolin_scores
# Output: pangolin_scores.csv (directly usable by score.py --pangolin-scores)
```

---

## Input Format

- **VCF** (`.vcf`): Standard VCF format
- **TSV** (`.tsv`): Tab-delimited with columns `chr`, `pos`, `ref`, `alt`
- **CSV** (`.csv`): Comma-delimited with columns `chr`, `pos`, `ref`, `alt`
- **Coordinates**: hg19/GRCh37 (default) or hg38/GRCh38 (`--genome hg38`).
  - hg19 path: GENCODE v19 GTF; phyloP/phastCons/gnomAD/GPN-MSA via
    `hg19tohg38/` forward liftover; DANN/Eigen/GERP matched directly.
  - hg38 path: GENCODE v49 GTF; phyloP/phastCons/gnomAD/GPN-MSA matched
    directly; DANN/Eigen/GERP via `hg38tohg19/` reverse liftover
    (one-time precompute — see *hg38 Input Mode* below).
- **Variants**: SNVs only (single-base ref and alt)

---

## Output

Tab-delimited file with columns:

| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `pos` | Position |
| `ref` | Reference allele |
| `alt` | Alternative allele |
| `region` | Splice region (one of the 7 model regions) |
| `gene` | Assigned gene |
| `splice_score` | Pathogenicity score (0-1) |

**Only variants classified into one of the 7 model regions are reported by
default** (Canonical Donor/Acceptor Site, Donor/Acceptor Region, PPT Region,
Branch Point Region, Deep Intronic). Variants that fall in exonic or
unclassified regions are excluded from the output, because the model is not
trained on those positions and the resulting scores are not meaningful.

Pass `--include-out-of-scope` to keep every input variant in the output; an
additional `in_scope` boolean column is added in that case, and out-of-scope
rows carry `region=NA` and a `splice_score` produced by the model on
out-of-distribution input (treat those values with caution).

---

## Score Interpretation

Operating points calibrated on training cross-validation (381,226 variants,
5-fold gene-grouped CV):

| Threshold | Use case | Sensitivity | Specificity |
|-----------|----------|-------------|-------------|
| >= 0.1 | Screening  | 99.1% | 98.5% |
| >= 0.2 | Balanced   | 98.8% | 98.9% |
| >= 0.5 | Diagnostic | 98.1% | 99.4% |

These are post hoc operating points; users planning prospective application
should re-calibrate on an internal validation set. For exonic variants (which
the default output excludes), use CADD v1.7 instead.

---

## Splice Regions

| Region | Definition |
|--------|-----------|
| Canonical Donor Site | 1-2 bp downstream of exon end |
| Canonical Acceptor Site | 1-2 bp upstream of exon start |
| Donor Region | 3-8 bp into intron from donor |
| Acceptor Region | 3-10 bp into intron from acceptor |
| PPT Region (putative) | 11-25 bp upstream of acceptor |
| Branch Point Region (putative) | 26-44 bp upstream of acceptor |
| Deep Intronic | >44 bp from nearest exon boundary |

The polypyrimidine-tract and branch-point regions are defined by distance
heuristics; transcript- or motif-specific branch points are not inferred.

---

## Citation

```
[Citation to be added upon publication]
```
