# Third-party licences

MetaSplice's own code, model and derived data are under the MetaSplice licence ([`LICENSE`](../LICENSE)). It grants
permission "to use, copy, modify, and distribute MetaSplice and its associated model files for non-commercial
research, educational, and clinical (non-revenue-generating) purposes, subject to the following conditions", which are
in [`LICENSE`](../LICENSE), and states: "For commercial use — including use within commercial products, fee-for-service
offerings, and internal use by for-profit entities — please contact xiaomingliu@usf.edu."

MetaSplice 2.0 also uses software and data from other providers, which keep their own terms. This page lists them,
with each provider's licence statement and a link to it, so that you can check the terms for your own use. A summary of
the notices is in [`THIRD_PARTY_NOTICES`](../THIRD_PARTY_NOTICES); the licence texts are in the Docker image and the
data bundle, at the paths given below. The statements below were checked against the linked pages in October 2026.
This page is not legal advice.

## Summary

| Component | Used for | Where it is | Licence |
|---|---|---|---|
| SpliceTransformer | 16 features | Docker image | Apache-2.0 |
| Pangolin 1.0.2 | 3 features | Docker image, with its source | GPL-3.0 |
| DANN scores | 1 feature | Data bundle | MIT |
| GPN-MSA scores | 1 feature | Data bundle | MIT |
| phyloP and phastCons (100-way, 470-way) | 4 features | Data bundle | UCSC: free for public and commercial use |
| GRCh38 and GRCh37 reference sequences | Reference | Data bundle | UCSC: free for public and commercial use |
| GENCODE v50 annotation | 10 features (exon context 3, region 7) and scope | Data bundle (GTF); Docker image (MANE databases) | Open access, no restrictions |
| Ensembl GRCh37–GRCh38 chain | Liftover tables built from it | Data bundle (the tables) | No restrictions |
| gnomAD pLI and LOEUF, ClinGen haploinsufficiency score, for HGNC's genes | 3 features | This repository (`data/`), the Docker image and the data bundle | Sources CC0 (gnomAD, ClinGen, HGNC); the file is under the MetaSplice licence |
| gnomAD v4.1 allele frequencies | Reported only, not a model input | Data bundle (optional) | CC0 |

The Python packages in the Docker image keep their own licences; the image lists them in
`/opt/metasplice/PYTHON_PACKAGE_LICENSES.txt`.

## Licence statements

**SpliceTransformer.** "This project is covered under the Apache 2.0 License."
(https://github.com/ShenLab-Genomics/SpliceTransformer). The image carries one change to its code, described in
[`docker_v2/patches/spt_torch26.patch`](../docker_v2/patches/spt_torch26.patch) and, in the image,
`/opt/SpliceTransformer/MODIFICATIONS.patch`.

**Pangolin.** GNU General Public License v3 (https://github.com/tkzeng/Pangolin/blob/main/LICENSE). Its complete
source, including the one change the image makes
([`docker_v2/patches/pangolin_softmask.patch`](../docker_v2/patches/pangolin_softmask.patch)), is in the image at
`/opt/Pangolin`. On program output, GPLv3 section 2 says: "The output from running a covered work is covered by this
License only if the output, given its content, constitutes a covered work."

**DANN.** The DANN readme states: "Both DANN and pre-computer [sic] DANN scores are under the MIT license." Its
licence file reads "MIT License / Copyright (c) [2015] [Daniel Quang]". The original site is no longer online; both
are available from the Internet Archive
(https://web.archive.org/web/20210414191342/https://cbcl.ics.uci.edu/public_data/DANN/readme,
https://web.archive.org/web/20240707115443/https://cbcl.ics.uci.edu/public_data/DANN/LICENSE), and the data bundle
includes the licence (`dann/LICENSE`).

**GPN-MSA.** "GPN is developed in the Song Lab at UC Berkeley and distributed under the MIT License."
(https://github.com/songlab-cal/gpn; licence: https://github.com/songlab-cal/gpn/blob/main/LICENSE). The score
dataset is also marked MIT (https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores). The data bundle includes
the licence (`gpn_msa/LICENSE`).

**UCSC conservation scores and reference sequences.** "As far as the UCSC Genome Browser group is concerned, the raw
table data and binary files used to create the graphics by the browser is freely available for both public and
commercial use" (https://genome.ucsc.edu/license/). The same page adds: "Sometimes, the source databases or authors
place restrictions on data." The 100-way README states: "All the files in this directory are freely available for
public use. For data use restrictions regarding the genome assemblies used in this annotation, see
http://genome.ucsc.edu/goldenPath/credits.html."
(https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/README.txt).

**GENCODE and Ensembl.** "All the GENCODE project data is open access"
(https://www.gencodegenes.org/pages/data_access.html). "Ensembl imposes no restrictions on access to, or use of, the
data provided and the software used to analyse and present it"
(https://www.ensembl.org/info/about/legal/disclaimer.html). The liftover tables are built from Ensembl's
`GRCh37_to_GRCh38.chain.gz` by
[`docker_v2/bundle/tools/build_liftover_tables.py`](../docker_v2/bundle/tools/build_liftover_tables.py).

**gnomAD.** "The primary data from the gnomAD exomes and genomes are available free of restrictions under the
Creative Commons Zero Public Domain Dedication" (https://gnomad.broadinstitute.org/policies). The statement does not
name the constraint metrics (pLI, LOEUF) separately.

**ClinGen.** "All curated content published by ClinGen is available free of restriction under the CC0 1.0 Universal
(CC0 1.0) Public Domain Dedication" (https://clinicalgenome.org/docs/terms-of-use/).

**HGNC.** HGNC states that "all data is released under the Creative Commons Public Domain (CC0) License. This means
that any form of reuse of the content is permitted." (https://www.genenames.org/about/license/)

The gene-level values in [`data/`](../data/README_gene_level.md) are built from gnomAD (pLI from v2.1.1, LOEUF from
v4.1), ClinGen and HGNC by [`data/build_gene_level.py`](../data/build_gene_level.py). Its ClinGen and HGNC inputs,
from 28 July 2026, are in [`data/sources/`](../data/sources/); it downloads the gnomAD files. The file of gene-level
values is distributed under the MetaSplice licence; its sources are CC0. The file in Docker images 2.0.0 and 2.0.1 and
in the data bundle is this same file (md5 `088ba4d5101c01a8697fbbf503affd06`), which `build_gene_level.py` reproduces
byte for byte.

## Licence conditions

- DANN and GPN-MSA: the MIT licence states that "The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software." The bundle has both notices (`dann/LICENSE`,
  `gpn_msa/LICENSE`).
- SpliceTransformer: the Apache License 2.0 requires that redistributions include a copy of the licence and that
  modified files carry notices of the change. The image has both (`/opt/SpliceTransformer/LICENSE`,
  `/opt/SpliceTransformer/MODIFICATIONS.patch`).
- Pangolin: GPL-3.0 requires anyone who conveys Pangolin to provide its corresponding source. The image includes it at
  `/opt/Pangolin`, with its licence and the notice of the change. Pangolin runs as a separate program, which the
  MetaSplice pipeline calls.
- A redistributed Docker image carries each component under its own terms: the MetaSplice licence for MetaSplice, the
  licences above for SpliceTransformer and Pangolin, their own licences for the Python packages (including NVIDIA's
  proprietary CUDA libraries) and the Debian-based Python base image, and the terms above for the data in it (the
  GENCODE annotation and the gene-level values).
- gnomAD: "All users of gnomAD data agree to not attempt to reidentify participants."

## Requests and cautions

- gnomAD: "we request that you actively acknowledge and give attribution to the gnomAD project, and link back to the
  relevant page, wherever possible." It asks developers whose tools include gnomAD data to say so, for example "This
  tool includes data from the gnomAD v4.1 release.", and asks that tools not use "gnomAD" or "Genome Aggregation
  Database" in their names. MetaSplice includes data from gnomAD v2.1.1 (pLI) and v4.1 (LOEUF, and the optional
  allele-frequency tables).
- ClinGen "requests that you give attribution to ClinGen and provide the date accessed whenever possible and
  appropriate". The values in `data/` are from the ClinGen list of 28 July 2026.
- HGNC: "Although stating the source is not mandatory, the HGNC recommends attribution".
- UCSC: "When you use the UCSC Genome Browser in your work, please cite one of our publications."
  (https://genome.ucsc.edu/license/). Its conditions of use state: "Data on the Genome Browser is meant for research
  purposes only. All clinical data should be interpreted by qualified professionals."
  (https://genome.ucsc.edu/conditions.html)

## Not included

**SpliceAI** and **CADD** were used only as comparators in the MetaSplice paper and are not distributed with MetaSplice.
SpliceAI's code is under the PolyForm Strict License 1.0.0, and its precomputed scores "are free for academic and
not-for-profit use; other use requires a commercial license from Illumina, Inc."
(https://github.com/Illumina/SpliceAI). "CADD scores are freely available for all non-commercial applications"
(https://cadd.bihealth.org/).
