"""Assemble the MetaSplice Release 2.0 annotation data bundle from the local inputs of the v2 analyses (maintainers).

  python assemble_bundle.py [--out /mnt/Work0/MetaSplice_release_v2/metasplice_v2_data] [--workers 12]

Hard-links every data file from the inputs recorded in reproducibility_v2/results/V00_inputs.txt (same filesystem;
nothing is copied or rewritten), copies the bundle's own README, licence and tool files, and writes MANIFEST.tsv:

  path  bytes  md5  content_md5  route  source  licence

md5 is reused from V00 where V00 recorded the same file and computed otherwise. content_md5, the md5 of the
uncompressed content, is given for the files that fetch_and_build.py rebuilds (chain and gnomAD tables) and for gzip
downloads, whose gzip stream can differ from the bundle's while the content is identical. route says how fetch_and_build.py obtains a file:
download (verbatim from source), gunzip (download source, decompress), index (built from another bundle file),
build (built by a tool in tools/ from source), local (copied from this bundle folder), mirror (only from the
MetaSplice mirror).
"""
import argparse, gzip, hashlib, os, shutil
from concurrent.futures import ProcessPoolExecutor

P = "/mnt/Work0/Splice_altering_variant"
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # the bundle folder (MetaSplice/docker_v2/bundle)
RAW = f"{P}/reproducibility/raw_inputs"
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]          # chrM is never in scope (no MANE gene on chrM)

UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath"
GPN = "https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/cf1718a926db35d078bcb08b9dcebc803494bf37"
CHAIN = "https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz"
GNOMAD = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint"
LIC = {"ucsc": "UCSC Genome Browser data: free for public and commercial use",
       "gencode": "GENCODE/Ensembl: no restrictions",
       "ensembl": "Ensembl: no restrictions (tables built by MetaSplice)",
       "dann": "MIT (dann/LICENSE)",
       "gpn": "MIT (gpn_msa/LICENSE)",
       "gnomad": "CC0 (tables built by MetaSplice)",
       "gene": "MetaSplice licence; built from gnomAD, ClinGen and HGNC (CC0)",
       "metasplice": "MetaSplice licence (LICENSE)"}


def entries():
    """(bundle path, local source, route, source, licence key) for every data file."""
    e = [("grch38/hg38.fa", f"{P}/hg38/hg38.fa", "gunzip", f"{UCSC}/hg38/bigZips/hg38.fa.gz", "ucsc"),
         ("grch38/hg38.fa.fai", f"{P}/hg38/hg38.fa.fai", "index", "pysam.faidx(grch38/hg38.fa)", "ucsc"),
         ("grch38/hg38.fa.fxi", f"{P}/hg38/hg38.fa.fxi", "index", "pyfastx.Fasta(grch38/hg38.fa) [Pangolin]", "ucsc"),
         ("grch37/hg19.fa", f"{P}/hg19/hg19.fa", "gunzip", f"{UCSC}/hg19/bigZips/hg19.fa.gz", "ucsc"),
         ("grch37/hg19.fa.fai", f"{P}/hg19/hg19.fa.fai", "index", "pysam.faidx(grch37/hg19.fa)", "ucsc"),
         ("gencode/gencode.v50.annotation.gtf.gz", f"{P}/docker_data/gencode.v50.annotation.gtf.gz", "download",
          "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz", "gencode")]
    for c in CHROMS:
        e.append((f"chain/hg19tohg38_ensembl/chr{c}.hg19tohg38.gz", f"{P}/hg19tohg38_ensembl/chr{c}.hg19tohg38.gz",
                  "build", f"tools/build_liftover_tables.py from {CHAIN}", "ensembl"))
        e.append((f"chain/hg38tohg19_ensembl/chr{c}.hg38tohg19.gz", f"{P}/hg38tohg19_ensembl/chr{c}.hg38tohg19.gz",
                  "build", f"tools/build_liftover_tables.py from {CHAIN}", "ensembl"))
    for c in CHROMS:
        e.append((f"dann/DANN_whole_genome_SNVs.tsv.gz.chr{c}.gz", f"{RAW}/DANN/DANN_whole_genome_SNVs.tsv.gz.chr{c}.gz",
                  "mirror", "DANN_whole_genome_SNVs.tsv.bgz (cbcl.ics.uci.edu/public_data/DANN, offline), split by chromosome",
                  "dann"))
    cons = [("phyloP_GRCh38", "phyloP100way", f"{UCSC}/hg38/phyloP100way/hg38.100way.phyloP100way"),
            ("phyloP_GRCh38", "phyloP470way", None),
            ("phastCons_GRCh38", "phastCons100way", f"{UCSC}/hg38/phastCons100way/hg38.100way.phastCons"),
            ("phastCons_GRCh38", "phastCons470way", f"{UCSC}/hg38/phastCons470way/hg38.470way.phastCons")]
    for d, kind, url in cons:
        for c in CHROMS:
            f = f"chr{c}.{kind}.wigFix.gz"
            if url:
                e.append((f"conservation/{d}/{f}", f"{RAW}/{d}/{f}", "download", f"{url}/{f}", "ucsc"))
            else:                       # UCSC replaced these files on 2023-09-02 (README.md, "phyloP 470-way")
                e.append((f"conservation/{d}/{f}", f"{RAW}/{d}/{f}", "mirror",
                          f"{UCSC}/hg38/phyloP470way/hg38.470way.phyloP/{f} as released 2022-08-28 (md5 in UCSC md5sum.txt)",
                          "ucsc"))
    e += [("gpn_msa/scores.tsv.bgz", f"{RAW}/gpn_msa/scores.tsv.bgz", "download", f"{GPN}/scores.tsv.bgz", "gpn"),
          ("gpn_msa/scores.tsv.bgz.tbi", f"{RAW}/gpn_msa/scores.tsv.bgz.tbi", "download", f"{GPN}/scores.tsv.bgz.tbi", "gpn")]
    for c in CHROMS:
        f = f"gnomad.joint.v4.1.sites.chr{c}.vcf.bgz.snp.tsv.gz"
        e.append((f"gnomad/{f}", f"{RAW}/gnomAD/{f}", "build",
                  f"tools/build_gnomad_tables.py from {GNOMAD}/gnomad.joint.v4.1.sites.chr{c}.vcf.bgz", "gnomad"))
    e.append(("gene_level/metasplice_gene_level_dbNSFP5.4.tsv.gz", f"{P}/MetaSplice/data/metasplice_gene_level_dbNSFP5.4.tsv.gz",
              "local", "MetaSplice/data/build_gene_level.py from gnomAD, ClinGen and HGNC", "gene"))
    return e


OWN = [("README.md", "metasplice"), ("LICENSE", "metasplice"), ("dann/LICENSE", "dann"), ("dann/README.md", "dann"),
       ("gpn_msa/LICENSE", "gpn"), ("gene_level/README_gene_level.md", "gene"), ("tools/fetch_and_build.py", "metasplice"),
       ("tools/build_liftover_tables.py", "metasplice"), ("tools/build_gnomad_tables.py", "metasplice"),
       ("tools/assemble_bundle.py", "metasplice"), ("tools/make_mirror_links.py", "metasplice")]


def own_files():
    """(bundle path, local file, licence key) for the bundle's own documents and tools, kept in this bundle folder
    (some as symbolic links to MetaSplice/LICENSE, MetaSplice/data and the chain builder of the v1 deposit)."""
    return [(rel, os.path.join(HERE, rel), lic) for rel, lic in OWN]


def content_checked(rel, route):
    """Files checked by their uncompressed content: built tables, and gzip downloads (a host may re-compress a file
    without changing it, as UCSC did with phastCons 470-way on 2023-09-29)."""
    return route == "build" or (route == "download" and rel.endswith(".gz"))


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 23), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_content(path):
    h = hashlib.md5()
    with gzip.open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 23), b""):
            h.update(chunk)
    return h.hexdigest()


def v00():
    out = {}
    for line in open(f"{P}/reproducibility_v2/results/V00_inputs.txt"):
        p = line.rstrip("\n").split("\t")
        if len(p) == 5:
            out[p[4]] = (p[1], int(p[2]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/mnt/Work0/MetaSplice_release_v2/metasplice_v2_data")
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()
    data, own = entries(), own_files()
    for _, src, *_ in data + own:
        if not os.path.isfile(src):
            raise FileNotFoundError(src)
    for rel, src, *_ in data:
        dst = os.path.join(args.out, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if os.path.exists(dst):
            if os.path.samefile(dst, src):
                continue
            raise FileExistsError(f"{dst} exists and is not {src}")
        os.link(src, dst)
    for rel, src, _ in own:
        dst = os.path.join(args.out, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(src, dst)

    recorded = v00()
    previous = {}                       # content md5s of an earlier run, reused when path, size and md5 are unchanged
    old = os.path.join(args.out, "MANIFEST.tsv")
    if os.path.exists(old):
        for line in open(old).readlines()[1:]:
            p = line.rstrip("\n").split("\t")
            previous[(p[0], int(p[1]), p[2])] = p[3]
    need_md5 = [src for _, src, *_ in data if os.path.realpath(src) not in recorded] + [src for _, src, _ in own]
    need_content = [src for rel, src, route, *_ in data if content_checked(rel, route) and previous.get(
        (rel, os.path.getsize(src), recorded.get(os.path.realpath(src), ("", 0))[0]), ".") == "."]
    with ProcessPoolExecutor(args.workers) as ex:
        md5 = dict(zip(need_md5, ex.map(md5_file, need_md5)))
        content = dict(zip(need_content, ex.map(md5_content, need_content)))
    rows = []
    for rel, src, route, source, lic in data:
        m, size = recorded.get(os.path.realpath(src), (md5.get(src), os.path.getsize(src)))
        assert size == os.path.getsize(src), src
        rows.append((rel, size, m, content.get(src) or previous.get((rel, size, m), "."), route, source, LIC[lic]))
    for rel, src, lic in own:
        rows.append((rel, os.path.getsize(src), md5[src], ".", "local", "MetaSplice", LIC[lic]))
    with open(os.path.join(args.out, "MANIFEST.tsv"), "w") as f:
        f.write("path\tbytes\tmd5\tcontent_md5\troute\tsource\tlicence\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    total = sum(r[1] for r in rows)
    print(f"{len(rows)} files, {total / 1e9:,.1f} GB ({total:,} bytes); MANIFEST.tsv written")
    for route in ("download", "gunzip", "index", "build", "local", "mirror"):
        sel = [r for r in rows if r[4] == route]
        print(f"  {route:<8} {len(sel):>4} files  {sum(r[1] for r in sel) / 1e9:8.1f} GB")


if __name__ == "__main__":
    main()
