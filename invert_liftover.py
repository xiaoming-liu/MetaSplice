"""
Invert hg19->hg38 liftover files into hg38->hg19 mapping files.

Reads raw_inputs/hg19tohg38/chr{N}.hg19tohg38.gz (which has columns
chr_hg19, pos_hg19, chr_hg38, pos_hg38, sorted by hg19 position) and writes
raw_inputs/hg38tohg19/chr{M}.hg38tohg19.gz with columns chr_hg38, pos_hg38,
chr_hg19, pos_hg19 sorted by hg38 position. Output is grouped by hg38 chromosome
(M), which is not always the same as the source hg19 chromosome (rare but
possible at e.g. centromere or alt loci).

This is a one-time pre-compute step required to support hg38-input MetaSplice
scoring (Docker score.py --genome hg38).

Usage:
  python scripts/invert_liftover.py
  python scripts/invert_liftover.py --input_dir raw_inputs/hg19tohg38 \\
                                    --output_dir raw_inputs/hg38tohg19
"""
import argparse
import gzip
import os
import sys
import time
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description="Invert hg19->hg38 liftover")
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    parser.add_argument(
        "--input_dir",
        default=os.path.join(base, "raw_inputs/hg19tohg38"),
    )
    parser.add_argument(
        "--output_dir",
        default=os.path.join(base, "raw_inputs/hg38tohg19"),
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # First pass: collect all hg38 -> hg19 records grouped by hg38 chromosome
    # Memory note: each hg19 chromosome holds ~250M base mappings. Storing them
    # as Python lists is ~12 GB total. We stream chromosome-by-chromosome and
    # write per-hg38-chromosome temp files, then sort and merge.
    print(f"Inverting liftover from {args.input_dir} -> {args.output_dir}")

    # Open output handles (one per hg38 chromosome) - keep open for append
    # We accumulate (hg38_pos, hg19_chr, hg19_pos) tuples per hg38 chromosome,
    # process one input chromosome at a time, write its records, then proceed.
    # Final sort happens per output chromosome at the end.

    # Stage 1: per-hg38-chromosome temp files (text)
    temp_dir = os.path.join(args.output_dir, "_tmp")
    os.makedirs(temp_dir, exist_ok=True)
    temp_handles = {}

    def get_temp(hg38_chr):
        if hg38_chr not in temp_handles:
            temp_handles[hg38_chr] = open(
                os.path.join(temp_dir, f"chr{hg38_chr}.tmp"), "w"
            )
        return temp_handles[hg38_chr]

    chroms = list(range(1, 23)) + ["X", "Y"]
    total_records = 0
    for chrom_num in chroms:
        path = os.path.join(args.input_dir, f"chr{chrom_num}.hg19tohg38.gz")
        if not os.path.exists(path):
            print(f"  chr{chrom_num}: input not found, skipping")
            continue
        t0 = time.time()
        n = 0
        with gzip.open(path, "rt") as f:
            f.readline()  # skip header
            for line in f:
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                try:
                    hg19_pos = int(parts[1])
                    hg38_chr = parts[2].strip()
                    hg38_pos = int(parts[3])
                except ValueError:
                    continue
                get_temp(hg38_chr).write(
                    f"{hg38_pos}\t{chrom_num}\t{hg19_pos}\n"
                )
                n += 1
        total_records += n
        print(f"  chr{chrom_num}: {n:,} records ({time.time()-t0:.0f}s)")

    for h in temp_handles.values():
        h.close()
    print(f"\nTotal records collected: {total_records:,}")
    print(f"Output hg38 chromosomes: {sorted(temp_handles.keys(), key=lambda x: (len(x), x))}")

    # Stage 2: sort each temp file by hg38_pos and write compressed output
    print("\nSorting and compressing per hg38 chromosome:")
    for hg38_chr in sorted(temp_handles.keys(), key=lambda x: (len(x), x)):
        tmp_path = os.path.join(temp_dir, f"chr{hg38_chr}.tmp")
        out_path = os.path.join(args.output_dir, f"chr{hg38_chr}.hg38tohg19.gz")
        t0 = time.time()
        rows = []
        with open(tmp_path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    try:
                        rows.append((int(parts[0]), parts[1], int(parts[2])))
                    except ValueError:
                        continue
        rows.sort(key=lambda r: r[0])
        with gzip.open(out_path, "wt") as f:
            f.write("#chr_hg38\tpos_hg38\tchr_hg19\tpos_hg19\n")
            for hg38_pos, hg19_chr, hg19_pos in rows:
                f.write(f"{hg38_chr}\t{hg38_pos}\t{hg19_chr}\t{hg19_pos}\n")
        os.remove(tmp_path)
        print(f"  chr{hg38_chr}: {len(rows):,} records ({time.time()-t0:.0f}s)")

    os.rmdir(temp_dir)
    print("\nDONE")


if __name__ == "__main__":
    main()
