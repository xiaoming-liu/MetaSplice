"""
Step 0: Build the per-position hg19 <-> hg38 tables from a chain file.

Writes, in the formats read by score.py, 16_compute_test_features.py and
invert_liftover.py:
  <out19>/chr{N}.hg19tohg38.gz   #chr_hg19 pos_hg19 chr_hg38 pos_hg38, sorted by hg19 position
  <out38>/chr{M}.hg38tohg19.gz   #chr_hg38 pos_hg38 chr_hg19 pos_hg19, sorted by hg38 position

Positions are 1-based and chromosome names carry no "chr" prefix. Each chain block maps a run
of consecutive positions; a block on the reverse query strand maps them in descending order.
As before, the tables carry no strand column. For each position this gives what CrossMap
reports for a single-base interval with the same chain file.

The default chain is Ensembl's GRCh37_to_GRCh38 (no restrictions on use), replacing the UCSC
hg19ToHg38 chain behind the earlier tables (non-commercial use only). Ensembl's chain maps
each GRCh37 chromosome to the same GRCh38 chromosome only. chrM is not built: Ensembl's MT is
the rCRS sequence, whereas hg19's chrM is NC_001807, so the chain does not apply to it.

Usage:
  curl -O https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz
  python 00_build_liftover_tables.py --chain GRCh37_to_GRCh38.chain.gz \\
      --out19 raw_inputs/hg19tohg38 --out38 raw_inputs/hg38tohg19
"""
import argparse
import gzip
import os
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np

PRIMARY = [str(i) for i in range(1, 23)] + ["X", "Y"]

# Expand one block per line into one row per position; printf keeps integers exact.
EXPAND_19 = ('{for(i=0;i<$3;i++){p=($6=="+")?$5+i+1:$7-($5+i);'
             'printf "%s\\t%d\\t%s\\t%d\\n",$1,$2+i+1,$4,p}}')
EXPAND_38 = '{for(i=0;i<$2;i++) printf "%s\\t%d\\t%s\\t%d\\n",$5,$1+i,$3,$4+$6*i}'


def strip(name):
    return name[3:] if name.startswith("chr") else name


def parse_chain(path, chroms):
    """Blocks per target chromosome: (t_start0, size, q_chr, q_start0, q_strand, q_size).

    Ensembl's GRCh37_to_GRCh38 file carries every chain twice, under "1" and under "chr1"
    (identical content), so each chromosome is read under one name only: the bare name if
    the file has it, otherwise the "chr" name (UCSC chains)."""
    wanted = set(chroms) | {"chr" + c for c in chroms}
    blocks = defaultdict(list)
    with gzip.open(path, "rt") as f:
        cur = None
        for line in f:
            p = line.split()
            if not p:
                cur = None
                continue
            if p[0] == "chain":
                tname = p[2]
                meta = (strip(p[7]), p[9], int(p[8]))
                cur = [int(p[5]), int(p[10])] if tname in wanted else None
                continue
            if cur is None:
                continue
            size = int(p[0])
            blocks[tname].append((cur[0], size, meta[0], cur[1], meta[1], meta[2]))
            if len(p) == 3:
                cur[0] += size + int(p[1])
                cur[1] += size + int(p[2])
    return {c: sorted(blocks[c] if c in blocks else blocks.get("chr" + c, [])) for c in chroms}


def assert_no_overlap(intervals, what):
    a = np.array(sorted(intervals), dtype=np.int64)
    if len(a) > 1:
        ov = np.maximum(0, np.maximum.accumulate(a[:, 1])[:-1] - a[1:, 0]).sum()
        assert ov == 0, f"{what}: {ov} positions covered twice; a per-position table would be ambiguous"


def run(block_file, header, program, out_path):
    cmd = (f"set -o pipefail; {{ printf '{header}\\n'; mawk -F'\\t' '{program}' {block_file}; }} "
           f"| gzip -6 > {out_path}.part && mv {out_path}.part {out_path}")
    subprocess.run(["bash", "-c", cmd], check=True)
    return out_path


def main():
    ap = argparse.ArgumentParser(description="Build per-position hg19<->hg38 tables from a chain file")
    ap.add_argument("--chain", required=True)
    ap.add_argument("--out19", required=True)
    ap.add_argument("--out38", default=None, help="omit to build hg19->hg38 only")
    ap.add_argument("--chroms", default=",".join(PRIMARY), help="hg19 chromosomes to build")
    ap.add_argument("--jobs", type=int, default=12)
    args = ap.parse_args()
    chroms = args.chroms.split(",")
    blocks = parse_chain(args.chain, set(chroms))
    os.makedirs(args.out19, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="_blocks_", dir=args.out19)
    jobs = []

    for c in chroms:
        b = blocks.get(c, [])
        assert b, f"no chain blocks for {c}"
        assert_no_overlap([(t, t + s) for t, s, *_ in b], f"hg19 chr{c}")
        bf = os.path.join(tmp, f"19_{c}.tsv")
        with open(bf, "w") as f:
            for t, s, qc, qs, st, qsize in b:
                f.write(f"{c}\t{t}\t{s}\t{qc}\t{qs}\t{st}\t{qsize}\n")
        jobs.append((bf, "#chr_hg19\\tpos_hg19\\tchr_hg38\\tpos_hg38", EXPAND_19,
                     os.path.join(args.out19, f"chr{c}.hg19tohg38.gz")))

    if args.out38:
        os.makedirs(args.out38, exist_ok=True)
        by_q = defaultdict(list)
        for c in chroms:
            for t, s, qc, qs, st, qsize in blocks[c]:
                if st == "+":
                    by_q[qc].append((qs + 1, s, c, t + 1, 1))
                else:  # reverse strand: ascending hg38 positions pair with descending hg19 positions
                    by_q[qc].append((qsize - qs - s + 1, s, c, t + s, -1))
        for qc, b in by_q.items():
            b.sort()
            assert_no_overlap([(s38, s38 + s) for s38, s, *_ in b], f"hg38 {qc}")
            bf = os.path.join(tmp, f"38_{qc}.tsv")
            with open(bf, "w") as f:
                for s38, s, c19, s19, d in b:
                    f.write(f"{s38}\t{s}\t{c19}\t{s19}\t{qc}\t{d}\n")
            jobs.append((bf, "#chr_hg38\\tpos_hg38\\tchr_hg19\\tpos_hg19", EXPAND_38,
                         os.path.join(args.out38, f"chr{qc}.hg38tohg19.gz")))

    with ThreadPoolExecutor(args.jobs) as ex:
        for out in ex.map(lambda j: run(*j), jobs):
            print(f"  wrote {out}", flush=True)
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))
    os.rmdir(tmp)
    n = {c: sum(s for _, s, *_ in blocks[c]) for c in chroms}
    print(f"hg19 positions mapped: {sum(n.values()):,} across {len(chroms)} chromosomes")


if __name__ == "__main__":
    main()
