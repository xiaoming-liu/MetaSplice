"""hg19-input acceptance set: the acceptance variants in their GRCh37 coordinates.

For each variant of acceptance_expected.tsv, the GRCh37 coordinates stored with its features (chr19, pos19, ref19,
alt19: the GRCh37 ClinVar record with the same ClinVar ID) are written to acceptance_variants_hg19.vcf. The scorer's own
liftover (score_v2.lift_variants on the Ensembl-built hg19->hg38 tables) is applied here as well. Variants it returns to
their own GRCh38 position and alleles are expected to get their published score from an --genome hg19 run; they are
written to acceptance_expected_hg19.tsv. The others are counted by reason.

  python make_hg19_set.py /path/to/metasplice_v2_data

Maintainers only: it reads the study's working tree (the paths below). The fixtures it wrote are committed, and
the tests need only them.
"""
import importlib.util, os, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = "/mnt/Work0/Splice_altering_variant/reproducibility_v2"
KEY = ["chr", "pos", "ref", "alt"]
spec = importlib.util.spec_from_file_location("score_v2", "/mnt/Work0/Splice_altering_variant/MetaSplice/score_v2.py")
ms = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ms)


def main(data):
    e = pd.read_csv(os.path.join(HERE, "acceptance_expected.tsv"), sep="\t")
    k19 = pd.concat([pd.read_csv(f"{V2}/data/{f}", usecols=KEY + ["chr19", "pos19", "ref19", "alt19"], low_memory=False)
                     for f in ("v2_features_training.csv", "v2_features_pooled_test.csv")]).drop_duplicates(KEY)
    m = e.merge(k19, on=KEY, how="left", validate="one_to_one")
    has = m[m.pos19.notna()].copy()
    has["pos19"] = has.pos19.astype(int)
    raw = list(zip(has.chr19, has.pos19, has.ref19, has.alt19))
    lifted, failed = ms.lift_variants(raw, f"{data}/chain/hg19tohg38_ensembl", f"{data}/grch38/hg38.fa",
                                      source_fasta=f"{data}/grch37/hg19.fa")
    back = [lifted[i][:4] == (c, p, r, a) if i in lifted else False
            for i, (c, p, r, a) in enumerate(zip(has.chr, has.pos, has.ref, has.alt))]
    has["back"] = back
    print(f"{len(e):,} acceptance variants; GRCh37 coordinates for {len(has):,}; lifted back to the same GRCh38 variant "
          f"{sum(back):,}; not lifted {len(failed):,}; lifted elsewhere {len(lifted) - sum(back):,}")
    for reason, n in pd.Series(list(failed.values())).value_counts().items():
        print(f"  not lifted: {reason}: {n}")
    o = has.assign(o=has.chr19.map({f"chr{c}": k for k, c in enumerate([str(x) for x in range(1, 23)] + ["X", "Y"])}))
    o = o.sort_values(["o", "pos19"])
    with open(os.path.join(HERE, "acceptance_variants_hg19.vcf"), "w") as f:
        f.write("##fileformat=VCFv4.2\n##reference=GRCh37\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in o.itertuples(index=False):
            f.write(f"{r.chr19}\t{r.pos19}\t.\t{r.ref19}\t{r.alt19}\t.\tPASS\t.\n")
    has.loc[has.back, e.columns].to_csv(os.path.join(HERE, "acceptance_expected_hg19.tsv"), sep="\t", index=False,
                                        float_format="%.8f")


if __name__ == "__main__":
    main(sys.argv[1])
