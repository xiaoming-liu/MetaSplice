"""Check the regression runs of run_regressions.sh (fixed in release 2.0.0 before publication; verification report of
2026-10-11).

  python check_regressions.py RESULTS_DIR [--gpu]       (the folder run_regressions.sh writes to)

1. hg19 input REF: an input whose REF is wrong for GRCh37 but matches GRCh38 (chr19:11221454 C>G; hg19 has T, hg38 C)
   is not lifted and not scored; the correct REF (T>G) is not lifted because the base differs between builds; a valid
   hg19 variant is lifted and scored.
2. Score precision: splice_score is written at full precision. It equals the model's prediction from the exported
   features exactly, so the exact PP3/BP4 thresholds give the model's own assignments. For MEN1 chr11:64808119 T>C
   the model's score, 0.0017024071421474218, lies 2.6e-8 above the BP4 threshold; six decimals (0.001702) fell below
   it. With --gpu the exported score must also equal that value exactly. On the CPU, SpliceTransformer can differ
   slightly (README), so only the equality with the model's prediction is required.
   The GRCh38 run keeps its work directory under a relative path (--keep-work regression_work_relative, from
   /work): SpliceTransformer must read its input and write its scores there (review of 2026-10-14).
3. An input with no SNVs gives a table with the header and no rows, and exit 0 (review of 2026-10-13).
4. With an hg38.fa.fai older than hg38.fa on the read-only /data, the run still completes and scores MEN1
   (review of 2026-10-13; pyfaidx used to try to rewrite the index and SpliceTransformer failed). The run also has a
   read-only root filesystem (docker --read-only), as under Apptainer/Singularity: nothing is written outside the
   work directory (re-review of 2026-10-13).
Exits 1 on any failure.
"""
import os
import sys
import joblib
import numpy as np
import pandas as pd

BP4 = 0.0017023810017853637
MEN1 = 0.0017024071421474218
FEATURES = "/opt/metasplice/model/metasplice_v2_features.joblib"
MODEL = "/opt/metasplice/model/metasplice_v2_full.joblib"


def main(results, gpu=False):
    hg19_out, grch38_out, grch38_features = (os.path.join(results, f) for f in (
        "regression_out_hg19.tsv", "regression_out_grch38.tsv", "regression_features_grch38.tsv"))
    ok = True

    def check(name, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))

    print("1. hg19 input: REF checked against GRCh37 before liftover")
    h = pd.read_csv(hg19_out, sep="\t", dtype=str, keep_default_na=False)
    row = lambda p, r, a: h[(h.chr == "chr19") & (h.pos == p) & (h.ref == r) & (h.alt == a)]
    wrong, differs = row("11221454", "C", "G"), row("11221454", "T", "G")
    valid = h[(h.chr == "chr1") & (h.pos == "879191")]
    check("REF wrong for GRCh37 (C; hg19 has T) is not lifted", len(wrong) == 1 and
          wrong.status.iloc[0] == "not lifted: REF does not match hg19.fa" and wrong.splice_score.iloc[0] == "",
          f"{wrong[['status', 'splice_score']].to_dict('records')}")
    check("correct REF where the base differs between builds is not lifted", len(differs) == 1 and
          differs.status.iloc[0] == "not lifted: REF differs between builds" and differs.splice_score.iloc[0] == "",
          f"{differs[['status', 'splice_score']].to_dict('records')}")
    check("a valid hg19 variant is lifted and scored", len(valid) == 1 and valid.status.iloc[0] == "lifted (forward)"
          and valid.splice_score.iloc[0] != "", f"{valid[['chr38', 'pos38', 'status', 'splice_score']].to_dict('records')}")

    print("2. Score precision")
    o = pd.read_csv(grch38_out, sep="\t", float_precision="round_trip")    # the default parser is not exact
    f = pd.read_csv(grch38_features, sep="\t", float_precision="round_trip")
    feats = list(joblib.load(FEATURES))
    pred = joblib.load(MODEL).predict_proba(f[feats].astype(np.float32).values)[:, 1]
    m = o.merge(f[["chr", "pos", "ref", "alt"]].assign(model=[float(x) for x in pred]), on=["chr", "pos", "ref", "alt"])
    s = float(m.splice_score.iloc[0])
    check("exported score equals the model's prediction from the exported features exactly",
          len(m) == 1 and s == m.model.iloc[0], f"exported {s!r}, model {m.model.iloc[0]!r}")
    check("exported score has more than six decimals", len(repr(s).split(".")[1]) > 6, repr(s))
    check("BP4 by the exact threshold agrees with the model's prediction", (s <= BP4) == (m.model.iloc[0] <= BP4))
    if gpu:
        check("GPU: MEN1 score equals the published model score, above the BP4 threshold", s == MEN1 and s > BP4,
              f"{s!r}")
    kw = os.path.join(results, "regression_work_relative", "inscope.spt.csv")
    rows = open(kw).read().splitlines() if os.path.exists(kw) else []
    check("relative --keep-work: SpliceTransformer's scores are in the kept work directory", len(rows) == 2,
          f"{len(rows)} lines in {kw}")
    print("3. Input with no SNVs")
    e = os.path.join(results, "regression_out_empty.tsv")
    rows = open(e).read().splitlines() if os.path.exists(e) else []
    check("a header-only table is written", len(rows) == 1 and rows[0].endswith("splice_score"), f"{len(rows)} lines")

    print("4. hg38.fa.fai older than hg38.fa on the read-only bundle, read-only root filesystem")
    st = os.path.join(results, "regression_out_stale_index.tsv")
    t = pd.read_csv(st, sep="\t", float_precision="round_trip") if os.path.exists(st) else pd.DataFrame()
    check("the run completes and MEN1 is scored", len(t) == 1 and pd.notna(t.splice_score.iloc[0]),
          f"{t.to_dict('records') if len(t) else 'no output'}")
    if gpu and len(t):
        check("GPU: the same score as with a fresh index", float(t.splice_score.iloc[0]) == MEN1)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if x != "--gpu"]
    sys.exit(main(a[0], gpu="--gpu" in sys.argv))
