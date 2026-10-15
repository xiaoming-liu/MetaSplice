#!/usr/bin/env python3
"""Regression tests for the v2.0 scorer (score_v2.py; docs/METASPLICE_V2_PLAN.md, Task 2).

v2.0 scores GRCh38 natively with GENCODE v50 MANE Select + MANE Plus Clinical transcripts. These
tests pin down the behaviour the plan requires:
  1. only MANE transcripts define exons: a variant inside a non-MANE exon is in scope and
     classified relative to the MANE exons;
  2. a variant outside every MANE transcript of its gene is out of scope (decision O1);
  3. gene ties follow GTF order, as in training (v1.0.4);
  4. exon context uses the assigned gene's MANE exons with no distance limit;
  5. liftover keeps a forward segment's alleles, complements an inverted segment's, and refuses
     a REF that differs between builds, a position with two targets and a target with two sources;
     with the source check (hg19 input), it refuses a REF that is wrong for hg19 but matches hg38;
  6. DANN, the one hg19-only feature, is looked up at the lifted hg19 position with those alleles;
  7. gene-level features are looked up by Ensembl gene ID; a missing LOEUF is NaN;
  8. a variant not in gnomAD gets gnomad_af 0.0, and allele frequency is not a model input;
  9. the feature list is the 38 features, in model order (and equals the saved model's, once
     the model exists);
 10. a GTF without MANE tags is refused;
 11. a missing source raises instead of leaving features missing or allele frequency 0, and a
     variant that was not looked up in gnomAD gets NaN, not 0.0;
 12. scoring needs every source of the 38 features and the model; --features-only writes
     features without them;
 13. a Deep Intronic position inside an exon of a GENCODE basic transcript is out of scope (decision
     O5), while a MANE splice-region variant inside such an exon stays in scope; a GTF without basic
     transcripts outside MANE (the MANE subset) is refused by scope().
Like test_training_alignment.py, it needs no annotation databases: every input is synthetic.

Usage:  python test_v2.py
"""
import gzip
import importlib.util
import math
import os
import random
import sys
import tempfile

import subprocess
from argparse import Namespace

import numpy as np
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("score_v2", os.path.join(HERE, "score_v2.py"))
ms = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ms)

failures = []
COMP = str.maketrans("ACGT", "TGCA")


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))
    if not cond:
        failures.append(name)


TAGS = {"mane": ' tag "basic"; tag "MANE_Select";', "basic": ' tag "basic";', None: ""}


def gtf_line(chrom, feature, start, end, strand, gene, gene_id, tx, kind="mane"):
    attrs = f'gene_id "{gene_id}"; transcript_id "{tx}"; gene_name "{gene}";' + TAGS[kind]
    return f"{chrom}\tTEST\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n"


def write_gtf(path, genes, tags=True, mane_only=False):
    """genes: [(chrom, strand, name, gene_id, [(tx, kind, [(s, e), ...]), ...])], written in order;
    kind is "mane", "basic" or None. tags=False writes no tags; mane_only writes MANE transcripts only."""
    with gzip.open(path, "wt") as f:
        for chrom, strand, name, gid, txs in genes:
            for tx, kind, exons in txs:
                if mane_only and kind != "mane":
                    continue
                k = kind if tags else None
                f.write(gtf_line(chrom, "transcript", exons[0][0], exons[-1][1], strand, name, gid, tx, k))
                for s, e in exons:
                    f.write(gtf_line(chrom, "exon", s, e, strand, name, gid, tx, k))


GENES = [
    # GENEA: MANE exons 1000-1100, 2000-2100, 3000-3500; an untagged transcript adds exon 1500-1550;
    # a basic (not MANE) transcript adds exons 2500-2550 (deep in intron 2) and 2980-2995 (acceptor region of exon 3)
    ("chr1", "+", "GENEA", "ENSG00000000001.5", [("TA1", "mane", [(1000, 1100), (2000, 2100), (3000, 3500)]),
                                                 ("TA2", None, [(1000, 1100), (1500, 1550), (2000, 2100)]),
                                                 ("TA3", "basic", [(2500, 2550), (2980, 2995)])]),
    # ZZZ1 and AAA1 share both exons; ZZZ1 comes first in the GTF
    ("chr2", "+", "ZZZ1", "ENSG00000000026.2", [("TZ1", "mane", [(5000, 5100), (6000, 6100)])]),
    ("chr2", "+", "AAA1", "ENSG00000000077.1", [("TB1", "mane", [(5000, 5100), (6000, 6100)])]),
    # GENEX has exons far from the test position; GENEY's exon is nearer but belongs to another gene
    ("chr3", "+", "GENEX", "ENSG00000000031.1", [("TX1", "mane", [(1000, 1100), (50000, 50100)])]),
    ("chr3", "+", "GENEY", "ENSG00000000032.1", [("TY1", "mane", [(30000, 30100)])]),
    # GENEZ: positions before its first exon lie outside every MANE transcript
    ("chr4", "+", "GENEZ", "ENSG00000000041.3", [("TZZ", "mane", [(10000, 10100), (12000, 12100)])]),
]


def make_fastas(d):
    """hg19 and hg38 chr9 of 300 bp. hg19 11-40 -> hg38 111-140 (forward); hg19 61-90 -> hg38 190-161
    (inverted); hg19 101 and 102 -> hg38 50 (one target, two sources); hg19 150 -> two targets.
    hg38 120 carries a base that is neither hg19 20's base nor its complement."""
    rng = random.Random(3)
    h19 = [rng.choice("ACGT") for _ in range(300)]
    h38 = [rng.choice("ACGT") for _ in range(300)]
    table = []
    for p in range(11, 41):
        h38[p + 100 - 1] = h19[p - 1]; table.append((p, p + 100))
    for p in range(61, 91):
        h38[251 - p - 1] = h19[p - 1].translate(COMP); table.append((p, 251 - p))
    other = [b for b in "ACGT" if b not in (h19[19], h19[19].translate(COMP))]
    h38[119] = other[0]
    table += [(101, 50), (102, 50), (150, 200), (150, 201)]
    for name, seq in [("hg19", h19), ("hg38", h38)]:
        with open(f"{d}/{name}.fa", "w") as f:
            f.write(">chr9\n" + "".join(seq) + "\n")
        pysam.faidx(f"{d}/{name}.fa")
    os.makedirs(f"{d}/hg19tohg38"); os.makedirs(f"{d}/hg38tohg19")
    with gzip.open(f"{d}/hg19tohg38/chr9.hg19tohg38.gz", "wt") as f:
        f.write("#chr_hg19\tpos_hg19\tchr_hg38\tpos_hg38\n")
        for a, b in sorted(table):
            f.write(f"9\t{a}\t9\t{b}\n")
    with gzip.open(f"{d}/hg38tohg19/chr9.hg38tohg19.gz", "wt") as f:
        f.write("#chr_hg38\tpos_hg38\tchr_hg19\tpos_hg19\n")
        for a, b in sorted((b, a) for a, b in table):
            f.write(f"9\t{a}\t9\t{b}\n")
    return "".join(h19), "".join(h38)


with tempfile.TemporaryDirectory() as d:
    gtf = f"{d}/mane.gtf.gz"
    write_gtf(gtf, GENES)
    variants = [("chr1", 1525, "A", "G"),    # 0 inside the non-MANE exon 1500-1550
                ("chr1", 3200, "A", "G"),    # 1 inside MANE exon 3
                ("chr1", 1995, "A", "G"),    # 2 5 bp upstream of MANE exon 2: Acceptor Region
                ("chr2", 5995, "A", "G"),    # 3 tie between ZZZ1 and AAA1
                ("chr4", 9000, "A", "G"),    # 4 1 kb before GENEZ: outside every MANE transcript
                ("chr4", 9999, "A", "G"),    # 5 1 bp before GENEZ's first exon: outside too
                ("chr1", 2525, "A", "G"),    # 6 Deep Intronic, inside the basic exon 2500-2550 (O5)
                ("chr1", 2990, "A", "G")]    # 7 Acceptor Region of MANE exon 3, inside the basic exon 2980-2995
    regions, genes = ms.assign_regions(variants, gtf)
    inside = ms.in_mane_transcript(variants, genes, gtf)
    print("1. Only MANE transcripts define exons")
    check("variant inside a non-MANE exon is in scope, Deep Intronic relative to MANE exons",
          regions.get(0) == "Deep Intronic" and inside.get(0) is True, f"got {regions.get(0)}, inside {inside.get(0)}")
    check("variant inside a MANE exon is out of scope", 1 not in regions, f"got {regions.get(1)}")
    check("acceptor-region variant is classified", regions.get(2) == "Acceptor Region", f"got {regions.get(2)}")
    print("2. Outside every MANE transcript: out of scope (O1)")
    check("1 kb before the first exon is outside the gene's MANE transcripts", inside.get(4) is False, f"got {inside.get(4)}")
    check("1 bp before the first exon is outside too", inside.get(5) is False, f"got {inside.get(5)}")
    sc, gn, why = ms.scope(variants, gtf)
    check("scope() gives the reason 'outside MANE transcripts'", why.get(4) == why.get(5) == "outside MANE transcripts",
          f"got {why.get(4)!r}, {why.get(5)!r}")
    check("scope() keeps the in-scope variants", set(sc) == {0, 2, 3, 7}, f"got {sorted(sc)}")
    print("3. Gene ties follow GTF order")
    check("tie resolved to the first gene in the GTF, with its Ensembl ID",
          genes.get(3) == ("ZZZ1", "ENSG00000000026.2"), f"got {genes.get(3)}")

    print("4. Exon context: assigned gene's MANE exons, no distance limit")
    parsed = ms.parse_mane_gtf(gtf)
    ctx = ms.compute_exon_context([("chr3", 25000, "A", "G")], {0: "GENEX"}, parsed["exons_by_chr"])
    check("uses the assigned gene's exons, not a nearer exon of another gene", ctx.get(0, {}).get("exon_length") == 101.0,
          f"got {ctx.get(0)}")
    check("no distance limit (23,900 bp)", ctx.get(0, {}).get("dist_to_boundary") == 23900.0, f"got {ctx.get(0)}")

    print("5. Liftover with orientation from neighbouring positions")
    h19, h38 = make_fastas(d)
    b = lambda seq, p: seq[p - 1]
    alt_of = lambda r: "A" if r != "A" else "C"
    lv = [("chr9", 15, b(h19, 15), alt_of(b(h19, 15))),     # forward
          ("chr9", 70, b(h19, 70), alt_of(b(h19, 70))),     # inverted
          ("chr9", 20, b(h19, 20), alt_of(b(h19, 20))),     # REF differs between builds
          ("chr9", 101, b(h19, 101), alt_of(b(h19, 101))),  # target with two sources
          ("chr9", 150, b(h19, 150), alt_of(b(h19, 150)))]  # two targets
    lifted, failed = ms.lift_variants(lv, f"{d}/hg19tohg38", f"{d}/hg38.fa")
    check("forward segment keeps REF and ALT", lifted.get(0) == ("chr9", 115, lv[0][2], lv[0][3], 1), f"got {lifted.get(0)}")
    check("inverted segment complements REF and ALT",
          lifted.get(1) == ("chr9", 181, lv[1][2].translate(COMP), lv[1][3].translate(COMP), -1), f"got {lifted.get(1)}")
    check("REF that differs between builds is not lifted", failed.get(2) == "REF differs between builds", f"got {failed.get(2)}")
    check("a target reached from two sources is not lifted", 3 in failed and 3 not in lifted, f"got {failed.get(3)}")
    check("a position with two targets is not lifted", 4 in failed and 4 not in lifted, f"got {failed.get(4)}")
    # release 2.0.0: user input (--genome hg19) is checked against the source build first. hg19 20 maps forward to
    # hg38 120, whose base differs; an input REF equal to the hg38 base is wrong for hg19 and must not be lifted.
    wrong = ("chr9", 20, b(h38, 120), alt_of(b(h38, 120)))
    lifted2, failed2 = ms.lift_variants(lv + [wrong], f"{d}/hg19tohg38", f"{d}/hg38.fa", source_fasta=f"{d}/hg19.fa")
    check("with the source check, valid inputs lift as before", {k: lifted2.get(k) for k in lifted} == lifted,
          f"got {lifted2}")
    check("with the source check, a REF wrong for hg19 is not lifted although it matches hg38",
          failed2.get(5) == "REF does not match hg19.fa" and 5 not in lifted2, f"got {failed2.get(5)}")

    print("6. DANN at the lifted hg19 position")
    os.makedirs(f"{d}/DANN")
    with gzip.open(f"{d}/DANN/DANN_whole_genome_SNVs.tsv.gz.chr9.gz", "wt") as f:
        f.write("#Chrom\tPos\tRef\tAlt\tDANN_score\n")
        for p in range(1, 301):
            for k, a in enumerate(x for x in "ACGT" if x != b(h19, p)):
                f.write(f"9\t{p}\t{b(h19, p)}\t{a}\t{p / 1000 + k / 10:.4f}\n")
    v38 = [("chr9", 115, lv[0][2], lv[0][3]), ("chr9", 181, lv[1][2].translate(COMP), lv[1][3].translate(COMP))]
    dann, dfail = ms.compute_dann(v38, f"{d}/DANN", f"{d}/hg38tohg19", f"{d}/hg19.fa")
    want = lambda p, a: float(f"{p / 1000 + [x for x in 'ACGT' if x != b(h19, p)].index(a) / 10:.4f}")
    check("forward: DANN of the same alleles at hg19 15", dann.get(0) == want(15, lv[0][3]), f"got {dann.get(0)}")
    check("inverted: DANN of the complemented alleles at hg19 70", dann.get(1) == want(70, lv[1][3]), f"got {dann.get(1)}")

    print("7. Gene level by Ensembl gene ID")
    gene_file = f"{d}/gene.gz"
    with gzip.open(gene_file, "wt") as f:
        f.write("Gene_name\tEnsembl_gene\tgnomAD_pLI\tgnomAD_LOEUF\tClinGen_Haploinsufficiency_Score\n")
        f.write("GENEA_OLDNAME\tENSG00000000001\t0.10\t0.85\t.\n")                    # symbol differs; ID matches
        f.write("ZZZ1\tENSG00000000099;ENSG00000000026\t0.91\t.\t3\n")               # two IDs in one row; ZZZ1's second
    feats = ms.compute_gene_level({0: ("GENEA", "ENSG00000000001.5"), 1: ("ZZZ1", "ENSG00000000026.2"),
                                   2: ("AAA1", "ENSG00000000077.1")}, gene_file)
    check("found by Ensembl ID although the symbol differs", feats.get(0, {}).get("gnomAD_pLI") == 0.10, f"got {feats.get(0)}")
    check("present LOEUF is kept", feats.get(0, {}).get("gnomAD_LOEUF") == 0.85)
    check("missing ClinGen HI is NaN", math.isnan(feats.get(0, {}).get("ClinGen_Haploinsufficiency_Score", 0.0)))
    check("missing LOEUF is NaN (not 0.0)", math.isnan(feats.get(1, {}).get("gnomAD_LOEUF", 0.0)), f"got {feats.get(1)}")
    check("an ID listed second in a row is found", feats.get(1, {}).get("ClinGen_Haploinsufficiency_Score") == 3.0)
    check("an ID absent from the file gets no gene-level values", 2 not in feats, f"got {feats.get(2)}")

    print("10. A GTF without MANE tags is refused")
    plain = f"{d}/plain.gtf.gz"
    write_gtf(plain, GENES, tags=False)
    try:
        ms.assign_regions(variants, plain)
        check("GTF without MANE tags raises", False, "no error")
    except ValueError as e:
        check("GTF without MANE tags raises", "MANE" in str(e), str(e))

print("8. gnomAD allele frequency is an annotation, not a model input")
af, pm = ms.gnomad_values(3, {1: 0.001}, {1: 0.002})
check("a variant not observed in gnomAD gets 0.0", af[0] == 0.0 and pm[2] == 0.0, f"got {af}, {pm}")
check("allele frequency is not a model input", not {"gnomAD_AF", "gnomAD_POPMAX_AF"} & set(ms.FEATURES))

print("9. Feature list")
EXPECTED = (["phyloP100way_vertebrate", "phyloP470way_mammalian", "phastCons100way_vertebrate", "phastCons470way_mammalian",
             "gnomAD_pLI", "gnomAD_LOEUF", "ClinGen_Haploinsufficiency_Score", "DANN_score", "exon_length",
             "exon_length_mod3", "dist_to_boundary", "spt_orig"]
            + [f"spt_o_{t}" for t in ["Adipose Tissue", "Blood", "Blood Vessel", "Brain", "Colon", "Heart", "Kidney", "Liver",
                                      "Lung", "Muscle", "Nerve", "Small Intestine", "Skin", "Spleen", "Stomach"]]
            + ["pangolin_inc", "pangolin_dec", "pangolin_max", "gpn_msa_score"]
            + [f"r_{r}" for r in ["Acceptor Region", "Branch Point Region", "Canonical Acceptor Site", "Canonical Donor Site",
                                  "Deep Intronic", "Donor Region", "PPT Region"]])
check("the 38 features in model order", ms.FEATURES == EXPECTED, f"{len(ms.FEATURES)} features")
saved = os.path.join(HERE, "model", "metasplice_v2_features.joblib")
if os.path.exists(saved):
    import joblib
    check("the saved model's feature list equals it", list(joblib.load(saved)) == EXPECTED)
else:
    print("  (model not trained yet: saved feature list not checked)")

print("13. Deep Intronic inside a GENCODE basic exon is out of scope (O5)")
with tempfile.TemporaryDirectory() as d:
    gtf = f"{d}/full.gtf.gz"
    write_gtf(gtf, GENES)
    sc, gn, why = ms.scope(variants, gtf)
    check("Deep Intronic inside a basic exon: out of scope", 6 not in sc and why.get(6) == "inside a GENCODE basic exon",
          f"got {sc.get(6)}, {why.get(6)!r}")
    check("MANE splice-region variant inside a basic exon: in scope", sc.get(7) == "Acceptor Region", f"got {sc.get(7)}")
    check("Deep Intronic inside an exon that is not basic: in scope", sc.get(0) == "Deep Intronic", f"got {sc.get(0)}")
    mane_only = f"{d}/mane_only.gtf.gz"
    write_gtf(mane_only, GENES, mane_only=True)
    regions_m, _ = ms.assign_regions(variants, mane_only)
    check("the MANE subset still serves region assignment", regions_m.get(2) == "Acceptor Region")
    try:
        ms.scope(variants, mane_only)
        check("scope() refuses the MANE subset", False, "no error")
    except ValueError as e:
        check("scope() refuses the MANE subset", "full" in str(e) and "MANE subset" in str(e), str(e)[:90])

print("11. A missing source raises")
with tempfile.TemporaryDirectory() as d:
    empty = f"{d}/empty"; os.makedirs(empty)
    v = [("chr9", 20, "A", "C")]

    def raises(name, fn, exc):
        try:
            fn()
            check(name, False, "no error")
        except exc as e:
            check(name, True, str(e)[:90])

    raises("gnomAD file missing for a needed chromosome", lambda: ms.compute_gnomad(v, {0: ("9", 20)}, empty), FileNotFoundError)
    raises("conservation file missing", lambda: ms.compute_conservation(v, {0: ("9", 20)}, empty, empty), FileNotFoundError)
    raises("DANN file missing", lambda: ms.scan_dann({0: ("chr9", 20, "A", "C")}, empty), FileNotFoundError)
    raises("GPN-MSA file missing", lambda: ms.compute_gpn_msa(v, {0: ("9", 20)}, f"{d}/none.bgz"), FileNotFoundError)
    with open(f"{d}/gpn.tsv", "w") as f:
        f.write("9\t20\tA\tC\t1.5\n")
    pysam.tabix_compress(f"{d}/gpn.tsv", f"{d}/gpn.tsv.bgz")
    pysam.tabix_index(f"{d}/gpn.tsv.bgz", seq_col=0, start_col=1, end_col=1, force=True)
    check("GPN-MSA found for a covered contig", ms.compute_gpn_msa(v, {0: ("9", 20)}, f"{d}/gpn.tsv.bgz") == {0: 1.5})
    raises("GPN-MSA contig missing", lambda: ms.compute_gpn_msa([("chr21", 20, "A", "C")], {0: ("21", 20)},
                                                                 f"{d}/gpn.tsv.bgz"), ValueError)
    af, pm = ms.gnomad_values(3, {1: 0.001}, {1: 0.002}, looked_up={0, 1})
    check("gnomAD: looked up and not found 0.0; found kept; not looked up NaN",
          af[0] == 0.0 and af[1] == 0.001 and np.isnan(af[2]) and np.isnan(pm[2]), f"got {af}")

print("12. Scoring needs every source and the model; --features-only does not")
with tempfile.TemporaryDirectory() as d:
    full = {k: "x" for k in ms.SCORING_INPUTS}
    check("nothing missing when every source is given", ms.missing_scoring_inputs(Namespace(**full)) == [])
    part = dict(full, spt_scores="", pangolin_scores="")
    miss = ms.missing_scoring_inputs(Namespace(**part))
    check("SpliceTransformer and Pangolin reported missing", len(miss) == 2 and "SpliceTransformer" in miss[0]
          and "Pangolin" in miss[1], f"got {miss}")

    def run(*extra):
        r = subprocess.run([sys.executable, os.path.join(HERE, "score_v2.py"), "--input", f"{d}/in.vcf",
                            "--output", f"{d}/out.tsv", "--fasta", "x.fa", "--gtf", "x.gtf.gz", *extra],
                           capture_output=True, text=True)
        return r.returncode, r.stderr
    given = [a for k in ms.SCORING_INPUTS if k not in ("spt_scores", "pangolin_scores")
             for a in (f"--{k.replace('_', '-')}", "x")]
    code, err = run(*given)
    check("scoring without SpliceTransformer or Pangolin is refused", code == 2 and "SpliceTransformer" in err
          and "Pangolin" in err and "--features-only" in err, f"exit {code}")
    code, err = run(*given, "--spt-scores", "x", "--pangolin-scores", "x", "--model-dir", d)
    check("scoring without the model is refused", code == 2 and "metasplice_v2_full.joblib" in err, f"exit {code}")
    code, err = run("--features-only")
    check("--features-only needs --features-out", code == 2 and "--features-out" in err, f"exit {code}")

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("All checks passed.")
