#!/usr/bin/env python3
"""Regression test for duplicate handling in compute_gnomad.

v1.0.2 converted the annotation lookups (DANN, Eigen, dbNSFP) to map each
(pos, ref, alt) key to every input row needing it, so a VCF listing the same variant
twice is annotated twice. compute_gnomad was missed and still did:

    lookup[(pos, ref, alt)] = idx     # last writer wins
    ...
    idx = lookup[key]                 # one row annotated, the rest left missing

gnomAD_AF and gnomAD_POPMAX_AF are model features, so the two copies of a duplicated
variant received different feature vectors and therefore different predictions. Fixed in
v1.0.3.

Unlike test_eigen_regression.py this needs no annotation databases: it writes a small
gnomAD-format table to a temporary directory, so it runs anywhere in about a second.

Usage:  python test_gnomad_duplicates.py
"""
import gzip
import importlib.util
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("score", os.path.join(HERE, "score.py"))
score = importlib.util.module_from_spec(spec)
spec.loader.exec_module(score)

failures = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))
    if not cond:
        failures.append(name)


HEADER = ["chrom", "pos", "ref", "alt", "gnomAD_joint_AF", "gnomAD_joint_POPMAX_AF"]
ROWS = [
    ("chr1", 100000, "A", "G", "0.00123", "0.00456"),
    ("chr1", 200000, "C", "T", "0.02000", "0.03000"),
]


def write_gnomad(dirpath):
    path = os.path.join(dirpath, "gnomad.joint.v4.1.sites.chr1.vcf.bgz.snp.tsv.gz")
    with gzip.open(path, "wt") as f:
        f.write("\t".join(HEADER) + "\n")
        for r in ROWS:
            f.write("\t".join(str(x) for x in r) + "\n")
    return path


print("gnomAD duplicate-row regression test\n")

with tempfile.TemporaryDirectory() as tmp:
    write_gnomad(tmp)

    # Same variant listed three times, plus a distinct one and an unmatched one.
    # variants entries are (chrom, pos, ref, alt); hg38_map is {index: (chrom_num, pos)}.
    variants = [
        ("chr1", 100000, "A", "G"),   # 0  duplicate group
        ("chr1", 200000, "C", "T"),   # 1  distinct
        ("chr1", 100000, "A", "G"),   # 2  duplicate group
        ("chr1", 999999, "T", "C"),   # 3  absent from the table
        ("chr1", 100000, "A", "G"),   # 4  duplicate group
    ]
    hg38_map = {i: ("1", v[1]) for i, v in enumerate(variants)}

    af, pm = score.compute_gnomad(variants, hg38_map, tmp)

    dup_idx = [0, 2, 4]
    check("every copy of a duplicated variant is annotated (AF)",
          all(i in af for i in dup_idx),
          f"annotated {[i for i in dup_idx if i in af]} of {dup_idx}")
    check("every copy of a duplicated variant is annotated (POPMAX)",
          all(i in pm for i in dup_idx),
          f"annotated {[i for i in dup_idx if i in pm]} of {dup_idx}")
    check("duplicated copies receive identical AF",
          len({af.get(i) for i in dup_idx}) == 1,
          f"values {[af.get(i) for i in dup_idx]}")
    check("duplicated copies receive identical POPMAX",
          len({pm.get(i) for i in dup_idx}) == 1,
          f"values {[pm.get(i) for i in dup_idx]}")
    check("the value is the one in the table",
          af.get(0) == 0.00123 and pm.get(0) == 0.00456,
          f"AF={af.get(0)} POPMAX={pm.get(0)}")
    check("a distinct variant is still annotated",
          af.get(1) == 0.02 and pm.get(1) == 0.03,
          f"AF={af.get(1)} POPMAX={pm.get(1)}")
    check("a variant absent from the table stays unannotated",
          3 not in af and 3 not in pm)

    # Complement matching must behave the same way: the scorer accepts either strand.
    rc_variants = [("chr1", 100000, "T", "C"), ("chr1", 100000, "T", "C")]
    rc_map = {i: ("1", 100000) for i in range(2)}
    rc_af, _ = score.compute_gnomad(rc_variants, rc_map, tmp)
    check("complement-matched duplicates are both annotated",
          len(rc_af) == 2 and len(set(rc_af.values())) == 1,
          f"{len(rc_af)} annotated, values {sorted(set(rc_af.values()))}")

print()
if failures:
    sys.exit(f"{len(failures)} test(s) failed: {', '.join(failures)}")
print("all tests passed")
