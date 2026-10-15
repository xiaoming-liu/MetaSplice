# DANN precomputed scores (GRCh37)

Quang D, Chen Y, Xie X. DANN: a deep learning approach for annotating the pathogenicity of genetic variants.
*Bioinformatics* 2015;31:761-3. doi:10.1093/bioinformatics/btu703

## Licence

MIT. The DANN readme says: "Both DANN and pre-computer DANN scores are under the MIT license." The licence text,
Copyright (c) 2015 Daniel Quang, is in `LICENSE` next to this file, unchanged.

## Source

These files are the full-precision scores of `DANN_whole_genome_SNVs.tsv.bgz`, which DANN's authors distributed at
`https://cbcl.ics.uci.edu/public_data/DANN/data/`. That host no longer resolves (checked 2026-09-27 and 2026-10-10).
Internet Archive copies of the readme and the licence:

- https://web.archive.org/web/20210414191342/https://cbcl.ics.uci.edu/public_data/DANN/readme
- https://web.archive.org/web/20240707115443/https://cbcl.ics.uci.edu/public_data/DANN/LICENSE

The MetaSplice bundle carries them so that MetaSplice's `DANN_score` feature can be computed exactly. Copies truncated to
three decimals (for example GREEN-VARAN, Zenodo record 3957486) change MetaSplice's scores and are not a substitute.

## Format

One gzip file per GRCh37 chromosome, `DANN_whole_genome_SNVs.tsv.gz.chr{1..22,X,Y}.gz`, split by chromosome from
the single distributed file. Each line is tab-separated, sorted by position, with no header and with
CR LF line endings:

    chromosome (no "chr")   position (GRCh37, 1-based)   reference allele   alternative allele   DANN score

The md5 of each file is in the bundle's `MANIFEST.tsv`.
