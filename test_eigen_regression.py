#!/usr/bin/env python3
"""Regression tests for Eigen two-block scanning.

Eigen_hg19_combined.tab.chr*.gz stores all coding records first, then all noncoding
records, each block independently ascending. Intronic variants appear only in the second
block. Two failure modes are guarded here:

  v1.0.0 bug  -- the scan exited at the first position above max_pos, inside the coding
                 block, so intronic variants got no Eigen score at all.
  v1.0.1 bug  -- the exit was removed entirely for Eigen, so any unmatched variant forced
                 a read to end-of-file (~80 s on chr21, minutes on chr1).

Usage:  python test_eigen_regression.py /path/to/annotation/data
"""
import importlib.util, os, sys, time

if len(sys.argv) < 2:
    sys.exit("usage: test_eigen_regression.py <data_dir>   "
             "(must contain DANN/ and Eigen/)")
DATA = sys.argv[1]
DANN, EIGEN = os.path.join(DATA, "DANN"), os.path.join(DATA, "Eigen")
for d in (DANN, EIGEN):
    if not os.path.isdir(d):
        sys.exit(f"missing required directory: {d}")

spec = importlib.util.spec_from_file_location(
    "score", os.path.join(os.path.dirname(os.path.abspath(__file__)), "score.py"))
score = importlib.util.module_from_spec(spec)
spec.loader.exec_module(score)

# Deep intronic chr21 variants: these live in the noncoding block only.
# Verified present in the noncoding block of Eigen chr21.
INTRONIC = [("chr21", 15516934, "A", "T"), ("chr21", 15552721, "T", "C"),
            ("chr21", 15561344, "C", "G")]
UNMATCHED = ("chr21", 15516935, "N", "N")   # cannot match any record
SLOW_SECONDS = 40                            # v1.0.1 took ~80 s here; fixed takes ~5 s

failures = []

def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))
    if not cond:
        failures.append(name)

print("Eigen two-block regression tests\n")

# 1. narrow intronic input must reach the noncoding block (guards the v1.0.0 bug)
_, eigen = score.compute_dann_eigen(INTRONIC, DANN, EIGEN)
check("narrow intronic input reaches the noncoding block",
      len(eigen) == len(INTRONIC), f"{len(eigen)}/{len(INTRONIC)} scored")

# 2. an unmatched variant must not force a full-file read (guards the v1.0.1 bug)
t0 = time.time()
_, eigen = score.compute_dann_eigen(list(INTRONIC) + [UNMATCHED], DANN, EIGEN)
elapsed = time.time() - t0
check("unmatched variant terminates promptly",
      elapsed < SLOW_SECONDS, f"{elapsed:.1f}s (limit {SLOW_SECONDS}s)")
check("matched variants still scored alongside an unmatched one",
      len(eigen) == len(INTRONIC), f"{len(eigen)}/{len(INTRONIC)} scored")

# 3. duplicate input rows must all be annotated
dup = list(INTRONIC) + list(INTRONIC)
_, eigen = score.compute_dann_eigen(dup, DANN, EIGEN)
check("duplicate input rows all annotated",
      len(eigen) == len(dup), f"{len(eigen)}/{len(dup)} scored")

# 4. DANN must keep its early-exit optimisation (single sorted block)
t0 = time.time()
dann, _ = score.compute_dann_eigen(INTRONIC, DANN, "")
check("DANN still scored", len(dann) == len(INTRONIC), f"{len(dann)}/{len(INTRONIC)}")
check("DANN retains its early exit", time.time() - t0 < SLOW_SECONDS,
      f"{time.time()-t0:.1f}s")

print()
if failures:
    sys.exit(f"{len(failures)} test(s) failed: {', '.join(failures)}")
print("all tests passed")
