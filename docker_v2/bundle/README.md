# MetaSplice Release 2.0 annotation data bundle (`metasplice_v2_data`)

The GRCh38 annotation files that the MetaSplice scorer (`score_v2.py`) and the MetaSplice Docker image
(`ghcr.io/xiaoming-liu/metasplice:2.0.1`, called `metasplice:2.0.1` below) read. Mount the folder at `/data`:

    docker run --rm --gpus all --user $(id -u):$(id -g) -v /path/to/metasplice_v2_data:/data:ro -v "$PWD":/work \
        metasplice:2.0.1 --input /work/variants.vcf --output /work/scores.tsv

219 GB in 201 data files (195 GB without the optional gnomAD tables). `MANIFEST.tsv` lists every file with its size,
md5, source and licence.

## Getting the bundle

There are two routes. Both give the same files, each checked against `MANIFEST.tsv`.

1. **Download the whole bundle from the MetaSplice mirror** on USF Box: https://usf.box.com/s/44v8lpk6rxi8d6vdvw1032cfr5olmzdq. `tools/fetch_and_build.py`
   reads each file's direct link from `MIRROR_LINKS.tsv` and checks every file:

       python tools/fetch_and_build.py --out metasplice_v2_data --all-from-mirror [--with-gnomad]

   Box takes files up to 15 GB, so the GPN-MSA scores (37 GB) are stored in three parts, which the script
   reassembles. `--only dann/,conservation/` limits a run to some folders. USF Box shows the Python scripts in `tools/`
   as previews rather than serving them by direct link. Take them from the MetaSplice repository or the image
   (`/opt/metasplice/bundle/tools/`); the script never downloads them.

2. **Fetch and build from the original sources.** Run the same script without `--all-from-mirror`. It:
   - downloads each resource from its original host (UCSC, GENCODE, Ensembl, Hugging Face, gnomAD);
   - builds the indexes and the chain tables;
   - takes from the mirror only the files that are no longer available from their original host (DANN and phyloP
     470-way; see below).

       python tools/fetch_and_build.py --out metasplice_v2_data [--with-gnomad]

The script needs Python 3.8+, pysam, pyfastx, bash, gzip and mawk. The image has all of them:

    mkdir -p "$PWD/metasplice_v2_data"
    docker run --rm --user $(id -u):$(id -g) -v "$PWD/metasplice_v2_data":/data --entrypoint python \
        metasplice:2.0.1 /opt/metasplice/bundle/tools/fetch_and_build.py --out /data --all-from-mirror

Create the folder first, as above: a folder that Docker creates belongs to root, and the download then stops with
"Permission denied". If an earlier attempt left an empty `metasplice_v2_data` owned by root, remove it
(`sudo rmdir metasplice_v2_data`) or use another folder.

It can be rerun after an interruption; files that already pass their check are kept. When it finishes, it makes each
index (`.fai`, `.fxi`, `.tbi`) newer than its data file. To check a bundle obtained any other way, run
`python tools/fetch_and_build.py --out metasplice_v2_data --verify`. It also warns about an index older than its data
file, as a copy that does not keep timestamps can leave one; `touch` that index.

## Contents

| Folder | Files | Used for | Source | Licence |
|---|---|---|---|---|
| `grch38/` | `hg38.fa`, `.fai`, `.fxi` | Reference allele check; SpliceTransformer and Pangolin input; conservation and GPN-MSA lookups | UCSC `hg38.fa.gz`, decompressed; indexes built with pysam and pyfastx | UCSC |
| `grch37/` | `hg19.fa`, `.fai` | GRCh37 allele check for the DANN lookup; hg19 input | UCSC `hg19.fa.gz`, decompressed | UCSC |
| `gencode/` | `gencode.v50.annotation.gtf.gz` | Scope rules: MANE Select and MANE Plus Clinical transcripts, regions, exon context | GENCODE release 50 | GENCODE/Ensembl, no restrictions |
| `chain/` | per-position tables, GRCh38→GRCh37 and GRCh37→GRCh38, 24 chromosomes each | DANN lookup (GRCh37 scores); hg19 input | Built by `tools/build_liftover_tables.py` from Ensembl's `GRCh37_to_GRCh38.chain.gz` | Ensembl, no restrictions |
| `dann/` | `DANN_whole_genome_SNVs.tsv.gz.chr{1..22,X,Y}.gz` | `DANN_score` | DANN's precomputed scores; host offline, so mirror only (`dann/README.md`) | MIT (`dann/LICENSE`) |
| `conservation/` | phyloP and phastCons, 100-way vertebrate and 470-way mammalian, per chromosome | Four conservation features | UCSC `phyloP100way`, `phyloP470way`, `phastCons100way`, `phastCons470way` | UCSC |
| `gpn_msa/` | `scores.tsv.bgz`, `.tbi` | `gpn_msa_score` | Hugging Face `songlab/gpn-msa-hg38-scores`, revision `cf1718a` | MIT (`gpn_msa/LICENSE`) |
| `gene_level/` | `metasplice_gene_level_dbNSFP5.4.tsv.gz` | `gnomAD_pLI`, `gnomAD_LOEUF`, `ClinGen_Haploinsufficiency_Score` | Built from gnomAD, ClinGen and HGNC for the MANE genes (`gene_level/README_gene_level.md`) | MetaSplice licence; gnomAD (open data; CC0 named for its primary data), ClinGen and HGNC (CC0) |
| `gnomad/` (optional) | `gnomad.joint.v4.1.sites.chr{1..22,X,Y}.vcf.bgz.snp.tsv.gz` | Reported allele frequencies (`gnomad_af`, `gnomad_popmax_af`); not model inputs | Built by `tools/build_gnomad_tables.py` from gnomAD v4.1 joint VCFs | CC0 |
| `tools/` | `fetch_and_build.py`, `build_liftover_tables.py`, `build_gnomad_tables.py`, `assemble_bundle.py`, `make_mirror_links.py` | Building, mirroring and checking the bundle | MetaSplice | MetaSplice licence |

chrM is left out: no MANE transcript lies on chrM, so no chrM variant is in scope.

## Versions that change scores

MetaSplice reproduces its published scores only with these exact files. Each check below was made on 2026-10-11.

- **phyloP 470-way.** The bundle carries UCSC's release of 2022-08-28, which the v2 model was trained on. UCSC
  replaced the 24 files on 2023-09-02 without updating its `md5sum.txt`, which still lists the 2022 files (and matches
  the bundle's). The 2023 files keep every 2022 value and add positions: on chr21, 989,255 positions (2.7%) were
  added, and none of the 36,131,530 existing values changed. A variant at an added position has no phyloP 470-way
  value in the 2022 files but has one in the 2023 files, so its score can change. These files therefore come from
  the mirror in both routes.
- **DANN.** Full-precision scores. Copies truncated to three decimals change scores (`dann/README.md`).
- **Reference sequences.** `hg38.fa` and `hg19.fa` are UCSC's `hg38.fa.gz` and `hg19.fa.gz`, decompressed.
  `fetch_and_build.py` checks the md5 of the decompressed file. UCSC's `hg38.fa.gz` differs by 8 bytes from the
  copy the analyses used, which carries a file name in its gzip header; the compressed data are identical.
- **Chain tables.** Built from Ensembl's chain (md5 `2a862e70fa1afc08f1f3e7a8cc5cda14`). The script checks the
  uncompressed content, because the gzip stream can differ between gzip versions.
- **gnomAD tables.** These are PASS SNVs only, with CR LF line endings. Popmax is taken over all nine ancestry
  groups (not gnomAD's grpmax). On chrY, the builder's output equals the bundle's table byte for byte.

## `MANIFEST.tsv`

Tab-separated, one row per file:

| Column | Meaning |
|---|---|
| `path` | Path within the bundle |
| `bytes` | Size |
| `md5` | md5 of the file |
| `content_md5` | md5 of the uncompressed content, for files the script rebuilds (`.` otherwise) |
| `route` | How `fetch_and_build.py` obtains the file: `download`, `gunzip`, `index`, `build`, `local` or `mirror` |
| `source` | Original URL, or the tool and input it is built from |
| `licence` | Licence of the file |

## Licences and citations

MetaSplice's own files (tools, README, the gene-level extract) are under the MetaSplice licence (`LICENSE`). The
third-party data keep their own terms:

- **UCSC** (reference sequences, phyloP, phastCons): free for public and commercial use (https://genome.ucsc.edu/license/).
  Pollard KS, et al. *Genome Res.* 2010;20:110-21. Siepel A, et al. *Genome Res.* 2005;15:1034-50.
- **GENCODE**: Harrow J, et al. *Genome Res.* 2012;22:1760-74. **MANE**: Morales J, et al. *Nature*
  2022;604:310-5.
- **Ensembl**: the chain file is from the Ensembl assembly-mapping release; Ensembl imposes no restrictions.
- **DANN**: MIT. Quang D, et al. *Bioinformatics* 2015;31:761-3.
- **GPN-MSA**: MIT. Benegas G, et al. *Nat Biotechnol.* 2025;43:1960-5.
- **gnomAD** (allele frequencies, `gnomad/`): CC0. Chen S, et al. *Nature* 2024;625:92-100.
- **Gene-level values**: from gnomAD (pLI, LOEUF; open data, with CC0 named for its primary data), ClinGen
  (haploinsufficiency; CC0) and HGNC (genes; CC0), with the rules of the dbNSFP 5.4 gene file (Liu X, et al. *Genome
  Med.* 2020;12:103), whose values they reproduce (`gene_level/README_gene_level.md`).
