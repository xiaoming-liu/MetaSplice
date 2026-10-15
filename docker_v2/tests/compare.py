"""Compare the image's scores for the acceptance set with the published scores.

  python compare.py acceptance_expected.tsv acceptance_out.tsv

Pass: every expected variant is scored, >= 99% agree within 1e-4 and none differs by more than 0.01 (the V07 gate's
criterion). Exits 1 on failure.
"""
import sys
import numpy as np
import pandas as pd


def main(expected, observed):
    e = pd.read_csv(expected, sep="\t")
    o = pd.read_csv(observed, sep="\t")
    if "chr38" in o.columns:            # hg19 input: chr, pos, ref, alt are the input's; compare on GRCh38
        o = o.drop(columns=["chr", "pos", "ref", "alt"]).rename(
            columns={"chr38": "chr", "pos38": "pos", "ref38": "ref", "alt38": "alt"})
    m = e.merge(o[["chr", "pos", "ref", "alt", "region", "splice_score"]], on=["chr", "pos", "ref", "alt"],
                how="left", suffixes=("", "_image"), validate="one_to_one")
    missing = int(m.splice_score.isna().sum())
    region_diff = int((m.region_image.notna() & (m.region != m.region_image)).sum())
    d = (m.splice_score - m.expected).abs()
    within = float((d <= 1e-4).mean())
    print(f"variants {len(m):,}; scored {len(m) - missing:,}; region differs {region_diff}")
    print(f"|image - published|: within 1e-4 {within:.2%}; max {d.max():.2e}; over 0.01: {int((d > 0.01).sum())}")
    for (s, r), g in m.assign(d=d).groupby(["set", "region"]):
        print(f"  {s:<12} {r:<26} n={len(g):>4}  max {g.d.max():.2e}")
    ok = missing == 0 and region_diff == 0 and within >= 0.99 and (d > 0.01).sum() == 0
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
