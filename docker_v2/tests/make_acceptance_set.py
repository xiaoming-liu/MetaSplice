"""Acceptance set for the MetaSplice Release 2.0 image: seeded ClinVar variants with their published scores.

Sample (seed 42): up to 200 variants per non-canonical region and 50 per canonical site, from the training set and from
the pooled temporal test set (deduplicated). The expected score of each is the full model's prediction from the stored
38-feature tables, i.e. the published score. The image recomputes everything from the variant alone (fresh
SpliceTransformer and Pangolin runs, every lookup), and must agree within 1e-4 for >= 99% of variants with none over 0.01,
the criterion of the v2.0 equivalence gate (V07).

Writes acceptance_variants.vcf (GRCh38) and acceptance_expected.tsv next to this script. ClinVar variants are public
domain; the scores are MetaSplice's.

Maintainers only: it reads the study's working tree (the paths below). The fixtures it wrote are committed, and
the tests need only them.
"""
import os
import joblib
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = "/mnt/Work0/Splice_altering_variant/reproducibility_v2"
NONCANON = ["Donor Region", "Acceptor Region", "PPT Region", "Branch Point Region", "Deep Intronic"]
CANON = ["Canonical Donor Site", "Canonical Acceptor Site"]


def main():
    model = joblib.load(f"{V2}/models/metasplice_v2_full.joblib")
    feats = list(joblib.load(f"{V2}/models/metasplice_v2_features.joblib"))
    parts = []
    for name, path in [("training", "v2_features_training.csv"), ("pooled test", "v2_features_pooled_test.csv")]:
        d = pd.read_csv(f"{V2}/data/{path}", low_memory=False)
        d = d.sample(frac=1, random_state=42)
        take = pd.concat([d[d.region == r].head(200) for r in NONCANON] + [d[d.region == r].head(50) for r in CANON])
        take = take.assign(set=name, expected=model.predict_proba(take[feats].astype(np.float32).values)[:, 1])
        parts.append(take[["chr", "pos", "ref", "alt", "region", "gene", "set", "expected"]])
    a = pd.concat(parts).drop_duplicates(["chr", "pos", "ref", "alt"])
    order = {f"chr{c}": k for k, c in enumerate([str(x) for x in range(1, 23)] + ["X", "Y"])}
    a = a.assign(o=a.chr.map(order)).sort_values(["o", "pos", "ref", "alt"]).drop(columns="o")
    with open(os.path.join(HERE, "acceptance_variants.vcf"), "w") as f:
        f.write("##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in a.itertuples(index=False):
            f.write(f"{r.chr}\t{r.pos}\t.\t{r.ref}\t{r.alt}\t.\tPASS\t.\n")
    a.to_csv(os.path.join(HERE, "acceptance_expected.tsv"), sep="\t", index=False, float_format="%.8f")
    print(f"{len(a):,} variants:", a.groupby(["set", "region"]).size().to_dict())


if __name__ == "__main__":
    main()
