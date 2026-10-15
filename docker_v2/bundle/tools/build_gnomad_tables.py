"""Per-chromosome gnomAD v4.1 joint SNV tables of the MetaSplice Release 2.0 data bundle (gnomad/).

  python build_gnomad_tables.py --vcf-dir DIR --out-dir gnomad [--chroms 1,2,...] [--jobs 4]

Reads gnomad.joint.v4.1.sites.chr{N}.vcf.bgz (https://gnomad.broadinstitute.org/downloads, joint exome+genome sites,
CC0) and writes gnomad.joint.v4.1.sites.chr{N}.vcf.bgz.snp.tsv.gz, one row per PASS SNV in VCF order:

  chr pos ref alt gnomAD_joint_flag gnomAD_joint_{AC,AN,AF,nhomalt} gnomAD_joint_POPMAX_{AC,AN,AF,nhomalt}
  gnomAD_joint_{AFR,AMI,AMR,ASJ,EAS,FIN,MID,NFE,SAS}_{AC,AN,AF,nhomalt}

chr carries no "chr" prefix; AF values are written as %e; a field absent from the record is '.'; gnomAD_joint_flag is
always '.'; lines end in CR LF, as in the bundle. POPMAX is the first of the nine groups, in the order above, with the
highest AF among those with an AF; it is not gnomAD's grpmax, which leaves out AMI, ASJ, FIN and MID. The scorer reads
gnomAD_joint_AF and gnomAD_joint_POPMAX_AF and reports them; they are not model inputs.

The output's uncompressed content equals the bundle's tables (checked on chrY, 705,461 rows, 2026-10-11); the gzip
stream differs, so compare the uncompressed md5 (MANIFEST.tsv, content_md5).
"""
import argparse, gzip, os
from concurrent.futures import ProcessPoolExecutor

POPS = ["afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"]
KEYS = ["AC_joint", "AN_joint", "AF_joint", "nhomalt_joint"] + [f"{s}_joint_{p}" for p in POPS for s in ("AC", "AN", "AF", "nhomalt")]
HEADER = (["chr", "pos", "ref", "alt", "gnomAD_joint_flag"]
          + [f"gnomAD_joint_{s}" for s in ("AC", "AN", "AF", "nhomalt")]
          + [f"gnomAD_joint_POPMAX_{s}" for s in ("AC", "AN", "AF", "nhomalt")]
          + [f"gnomAD_joint_{p.upper()}_{s}" for p in POPS for s in ("AC", "AN", "AF", "nhomalt")])
WANTED = set(KEYS)
BASES = {"A", "C", "G", "T"}


def fmt(key, val):
    if val is None:
        return "."
    return f"{float(val):e}" if key.startswith("AF_") else val


def popmax(groups):
    """AC, AN, AF, nhomalt of the first group with the highest AF; '.' x 4 when no group has an AF."""
    best = None
    for i in range(0, len(groups), 4):
        if groups[i + 2] != "." and (best is None or float(groups[i + 2]) > float(groups[best + 2])):
            best = i
    return ["."] * 4 if best is None else groups[best:best + 4]


def build(vcf, out):
    n = 0
    with gzip.open(vcf, "rt") as f, gzip.open(out + ".part", "wt", compresslevel=6) as o:
        o.write("\t".join(HEADER) + "\r\n")
        for line in f:
            if line[0] == "#":
                continue
            chrom, pos, _, ref, alt, _, filt, info = line.split("\t", 8)[:8]
            if filt != "PASS" or ref not in BASES or alt not in BASES:
                continue
            vals = {}
            for kv in info.split(";"):
                k, _, v = kv.partition("=")
                if k in WANTED:
                    vals[k] = v.rstrip("\n")
            v = [fmt(k, vals.get(k)) for k in KEYS]
            o.write("\t".join([chrom[3:] if chrom.startswith("chr") else chrom, pos, ref, alt, "."] + v[:4] + popmax(v[4:]) + v[4:]) + "\r\n")
            n += 1
    os.replace(out + ".part", out)
    return out, n


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vcf-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--chroms", default=",".join([str(i) for i in range(1, 23)] + ["X", "Y"]))
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    jobs = [(os.path.join(args.vcf_dir, f"gnomad.joint.v4.1.sites.chr{c}.vcf.bgz"),
             os.path.join(args.out_dir, f"gnomad.joint.v4.1.sites.chr{c}.vcf.bgz.snp.tsv.gz"))
            for c in args.chroms.split(",")]
    for vcf, _ in jobs:
        if not os.path.exists(vcf):
            raise FileNotFoundError(vcf)
    with ProcessPoolExecutor(args.jobs) as ex:
        for out, n in ex.map(build, *zip(*jobs)):
            print(f"  {out}: {n:,} SNVs", flush=True)


if __name__ == "__main__":
    main()
