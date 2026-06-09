#!/usr/bin/env python3
"""Convert a VCF file to the CSV format Pangolin expects.

Pangolin outputs CSV when given CSV input, which MetaSplice score.py
can read directly. This avoids the VCF-to-CSV conversion step.

Usage:
  python prepare_pangolin_input.py input.vcf pangolin_input.csv

Then run Pangolin on the CSV:
  pangolin pangolin_input.csv hg19.fa annotation.db pangolin_scores
  # Output: pangolin_scores.csv  (directly usable by score.py --pangolin-scores)
"""
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Convert VCF to CSV for Pangolin input")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("output_csv", help="Output CSV for Pangolin")
    args = parser.parse_args()

    n = 0
    with open(args.input_vcf) as fin, open(args.output_csv, "w") as fout:
        fout.write("CHROM,POS,REF,ALT\n")
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, _, ref, alt = fields[:5]
            chrom = chrom.lstrip("chr")
            for a in alt.split(","):
                if len(ref) == 1 and len(a) == 1:
                    fout.write(f"{chrom},{pos},{ref},{a}\n")
                    n += 1
    print(f"Wrote {n} SNVs to {args.output_csv}")


if __name__ == "__main__":
    main()
