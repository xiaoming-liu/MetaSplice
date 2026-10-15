#!/usr/bin/env python3
"""MetaSplice Release 2.0, end to end: variants -> SpliceTransformer + Pangolin -> 38 features -> score.

The image's entry point. It runs the two splicing models exactly as the Release 2.0 model's training data were produced
(reproducibility_v2/scripts/build/V03_upstream_scores.py):
  - SpliceTransformer, from a run directory holding the GENCODE v50 MANE annotation database, with TF32 off
    (NVIDIA_TF32_OVERRIDE=0) and --reference hg38 --vcf True --raw_score True;
  - Pangolin 1.0.2, with the GENCODE v50 MANE database and -c CHROM,POS,REF,ALT.
Then it calls the unchanged scorer (score_v2.py) with their outputs and the annotation resources under --data.

Only in-scope variants are sent to the two models. Scope, liftover (--genome hg19) and the REF check use the scorer's own
functions, so the variants selected here are exactly those the scorer scores.

Usage (inside the image):
  metasplice --input variants.vcf --output scores.tsv [--genome hg19] [--data /data] [--cpu] [--processes N]
             [--features-out FILE] [--include-out-of-scope] [--keep-work DIR]
"""
import argparse, csv, importlib.util, os, shutil, subprocess, sys, tempfile, time
import pandas as pd
import pysam

APP = os.environ.get("METASPLICE_APP", "/opt/metasplice")
SPT_DIR = os.environ.get("METASPLICE_SPT_DIR", "/opt/SpliceTransformer")
SPT_RUN = os.environ.get("METASPLICE_SPT_RUN", "/opt/spt_run")
PANGOLIN_DB = os.environ.get("METASPLICE_PANGOLIN_DB", "/opt/annotation/gencode.v50.mane.db")

spec = importlib.util.spec_from_file_location("score_v2", os.path.join(APP, "score_v2.py"))
ms = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ms)

# Default locations of the annotation resources inside the data bundle (docs/V2_RELEASE_PLAN.md, Task B); each can be
# overridden by the scorer's environment variable.
BUNDLE = {
    "METASPLICE_FASTA": "grch38/hg38.fa",
    "METASPLICE_HG19_FASTA": "grch37/hg19.fa",
    "METASPLICE_GTF": "gencode/gencode.v50.annotation.gtf.gz",
    "METASPLICE_LIFTOVER": "chain/hg19tohg38_ensembl",
    "METASPLICE_RLIFTOVER": "chain/hg38tohg19_ensembl",
    "METASPLICE_DANN": "dann",
    "METASPLICE_PHYLOP": "conservation/phyloP_GRCh38",
    "METASPLICE_PHASTCONS": "conservation/phastCons_GRCh38",
    "METASPLICE_GPNMSA": "gpn_msa/scores.tsv.bgz",
    "METASPLICE_GNOMAD": "gnomad",
}
SCORER_ARG = {"METASPLICE_FASTA": "--fasta", "METASPLICE_HG19_FASTA": "--hg19-fasta", "METASPLICE_GTF": "--gtf",
              "METASPLICE_LIFTOVER": "--liftover-dir", "METASPLICE_RLIFTOVER": "--reverse-liftover-dir",
              "METASPLICE_DANN": "--dann-dir", "METASPLICE_PHYLOP": "--phylop-dir",
              "METASPLICE_PHASTCONS": "--phastcons-dir", "METASPLICE_GPNMSA": "--gpn-msa", "METASPLICE_GNOMAD": "--gnomad-dir",
              "METASPLICE_DBNSFP": "--dbnsfp-gene"}
OPTIONAL = {"METASPLICE_GNOMAD", "METASPLICE_LIFTOVER"}     # reported AF only; hg19 input only
KEEP = "--keep-work DIR keeps the full logs"


def resources(data_root):
    """{env var: path}: an explicit environment variable wins, else the bundle location under data_root."""
    out = {}
    for var, rel in BUNDLE.items():
        out[var] = os.environ.get(var) or os.path.join(data_root, rel)
    out["METASPLICE_DBNSFP"] = os.environ.get("METASPLICE_DBNSFP") or os.path.join(
        APP, "data", "metasplice_gene_level_dbNSFP5.4.tsv.gz")
    return out


def device(cpu):
    """The device SpliceTransformer and Pangolin will use (both run on PyTorch in this container)."""
    if cpu:
        return "CPU (--cpu)"
    import torch
    if torch.cuda.is_available():
        return f"GPU ({torch.cuda.get_device_name(0)})"
    return "CPU: no GPU is visible in the container (run docker with --gpus all to use one)"


def check_input(path):
    """Inputs the scorer's parser (score_v2.parse_input, pinned) would misread: a compressed file other than .vcf.gz, a
    byte-order mark, and a table (not a VCF) without the four columns it reads, in either case."""
    if path.endswith(".vcf.gz"):
        return
    with open(path, "rb") as f:
        head = f.read(3)
    if head.startswith(b"\x1f\x8b"):
        raise ValueError(f"{path}: compressed; only a VCF is read compressed (.vcf.gz). Use a .vcf.gz, or an "
                         "uncompressed .vcf, .tsv or .csv")
    if head == b"\xef\xbb\xbf":
        raise ValueError(f"{path}: starts with a byte-order mark (as Excel's \"CSV UTF-8\" writes), which would be read "
                         "as part of the first column's name. Save it as plain CSV, or remove the mark: "
                         "sed -i '1s/^\\xEF\\xBB\\xBF//' FILE")
    if path.endswith(".vcf"):
        return
    with open(path) as f:
        header = next(csv.reader(f, delimiter="\t" if path.endswith(".tsv") else ","), [])
    want = {"chr": ("chr", "CHROM", "#CHROM"), "pos": ("pos", "POS"), "ref": ("ref", "REF"), "alt": ("alt", "ALT")}
    missing = [k for k, names in want.items() if not set(names) & set(header)]
    if missing:
        sep = "tab" if path.endswith(".tsv") else "comma"
        raise ValueError(f"{path}: no {', '.join(missing)} column in the header ({sep}-separated; a table needs "
                         "columns chr, pos, ref and alt, or CHROM, POS, REF and ALT; or use a .vcf)")


def in_scope_grch38(args, res):
    """The unique in-scope GRCh38 (chr, pos, ref, alt), selected with the scorer's own functions."""
    skipped = {}
    check_input(args.input)
    raw = ms.parse_input(args.input, skipped)
    low = next((v for v in raw if v[1] < 1), None)
    if low:
        raise ValueError(f"{args.input}: position {low[1]} ({low[0]} {low[2]}>{low[3]}) is below 1; positions are 1-based")
    if skipped:
        print(f"  skipped {sum(skipped.values()):,} records that are not SNVs: "
              + ", ".join(f"{k} {v:,}" for k, v in sorted(skipped.items())), flush=True)
    if args.genome == "hg19":
        lifted, _ = ms.lift_variants(raw, res["METASPLICE_LIFTOVER"], res["METASPLICE_FASTA"],
                                     source_fasta=res["METASPLICE_HG19_FASTA"])
        variants = [lifted[i][:4] if i in lifted else None for i in range(len(raw))]
    else:
        fa = pysam.FastaFile(res["METASPLICE_FASTA"])
        refs = set(fa.references)
        variants = [v if v[0] in refs and fa.fetch(v[0], v[1] - 1, v[1]).upper() == v[2] else None for v in raw]
        fa.close()
    regions, _, _ = ms.scope(variants, res["METASPLICE_GTF"])
    order = {f"chr{c}": k for k, c in enumerate([str(x) for x in range(1, 23)] + ["X", "Y"])}
    uniq = sorted({tuple(variants[i]) for i in regions}, key=lambda v: (order.get(v[0], 99), v[1], v[2], v[3]))
    return raw, uniq


def chromosome_groups(raw, fasta, n):
    """Split input rows into at most n groups of whole chromosomes, balanced by chromosome length (the scorer's lookups
    read each chromosome's annotation files once, so their cost follows chromosome length, not the number of variants)."""
    length = {}
    if os.path.exists(fasta + ".fai"):
        for line in open(fasta + ".fai"):
            f = line.split("\t")
            length[f[0]] = int(f[1])
    rows = {}
    for i, v in enumerate(raw):
        rows.setdefault(v[0], []).append(i)
    groups, load = [[] for _ in range(n)], [0] * n
    for c in sorted(rows, key=lambda c: -length.get(c, 0)):
        k = load.index(min(load))
        groups[k] += rows[c]
        load[k] += length.get(c, 1)
    return [sorted(g) for g in groups if g]


def score_groups(args, raw, base_cmd, work):
    """Run the scorer on chromosome groups in parallel and reassemble its output in input order, as text, exactly as a
    single run would write it. Each group runs with --include-out-of-scope so its rows match its input one to one."""
    groups = chromosome_groups(raw, base_cmd[base_cmd.index("--fasta") + 1], args.processes)
    procs = []
    for k, rows in enumerate(groups):
        stem = os.path.join(work, f"group{k}")
        with open(f"{stem}.vcf", "w") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for i in rows:
                c, p, r, a = raw[i]
                f.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\n")
        cmd = base_cmd + ["--input", f"{stem}.vcf", "--output", f"{stem}.out.tsv", "--include-out-of-scope"]
        if args.features_out:
            cmd += ["--features-out", f"{stem}.features.tsv"]
        procs.append((k, subprocess.Popen(cmd, stdout=open(f"{stem}.log", "w"), stderr=subprocess.STDOUT)))
    for k, p in procs:
        if p.wait():
            raise SystemExit(f"scoring group {k} failed; the end of its log ({KEEP}):\n"
                             + open(os.path.join(work, f"group{k}.log")).read()[-2000:])
    for kind, target in [("out", args.output), ("features", args.features_out)]:
        if not target:
            continue
        parts = []
        for k, rows in enumerate(groups):
            t = pd.read_csv(os.path.join(work, f"group{k}.{kind}.tsv"), sep="\t", dtype=str, keep_default_na=False)
            if len(t) != len(rows):
                raise SystemExit(f"group {k}: {len(t)} output rows for {len(rows)} input rows")
            parts.append(t.assign(_row=rows))
        t = pd.concat(parts).sort_values("_row", kind="stable").drop(columns="_row")
        if not args.include_out_of_scope:
            t = t[t.in_scope == "True"].drop(columns=["in_scope", "out_of_scope_reason"])
        t.to_csv(target, sep="\t", index=False)
    n_out = sum(1 for _ in open(args.output)) - 1
    print(f"Saved: {args.output} ({n_out:,} rows; {len(groups)} scoring processes)")


def write_inputs(variants, stem, pangolin=True):
    """stem.vcf for SpliceTransformer (and the scorer) and stem.pangolin.csv for Pangolin, in V03's formats."""
    with open(f"{stem}.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, p, r, a in variants:
            f.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\n")
    if not pangolin:
        return
    with open(f"{stem}.pangolin.csv", "w") as f:
        f.write("CHROM,POS,REF,ALT\n")
        for c, p, r, a in variants:
            f.write(f"{c.replace('chr', '')},{p},{r},{a}\n")


def spt_run_dir(work, fasta):
    """SpliceTransformer's run directory for this run, inside the work directory, so that nothing is written outside it
    (the image then also runs with a read-only root filesystem: docker --read-only, Apptainer/Singularity; re-review of
    2026-10-13). SpliceTransformer reads data/data_package/ and model/ relative to it. The image's annotation files and
    model are linked, as is the GRCh38 FASTA in use. Its .fai is a fresh copy: pyfaidx rewrites an index older than its
    FASTA, which fails on the read-only bundle, and a bundle copied without timestamps can have such an index (review of
    2026-10-13). The copy's content is unchanged."""
    run = os.path.join(work, "spt_run")
    pkg = os.path.join(run, "data", "data_package")
    os.makedirs(pkg, exist_ok=True)
    os.makedirs(os.path.join(run, "pyensembl_cache"), exist_ok=True)
    links = {os.path.join(run, "model"): os.path.join(SPT_DIR, "model"), os.path.join(pkg, "hg38.fa"): fasta}
    for name in ("hg38.annotation.gtf.gz", "hg38.annotation.gtf.db"):
        links[os.path.join(pkg, name)] = os.path.join(SPT_RUN, "data", "data_package", name)
    links[os.path.join(pkg, "hg38.fa.fai")] = None
    for link, target in links.items():
        if os.path.lexists(link):
            os.remove(link)
        if target:
            os.symlink(target, link)
    shutil.copyfile(fasta + ".fai", os.path.join(pkg, "hg38.fa.fai"))
    return run


def run_tools(stem, fasta, cpu):
    """SpliceTransformer and Pangolin at the same time; TF32 off; raise if either fails."""
    run = spt_run_dir(os.path.dirname(stem), fasta)
    env = dict(os.environ, PYENSEMBL_CACHE_DIR=os.path.join(run, "pyensembl_cache"), NVIDIA_TF32_OVERRIDE="0")
    if cpu:
        env["CUDA_VISIBLE_DEVICES"] = ""
    jobs = {
        "SpliceTransformer": subprocess.Popen(
            [sys.executable, os.path.join(SPT_DIR, "sptransformer.py"), "--input", f"{stem}.vcf", "--output",
             f"{stem}.spt.csv", "--reference", "hg38", "--vcf", "True", "--raw_score", "True"],
            cwd=run, env=env, stdout=open(f"{stem}.spt.log", "w"), stderr=subprocess.STDOUT),
        "Pangolin": subprocess.Popen(
            ["pangolin", f"{stem}.pangolin.csv", fasta, PANGOLIN_DB, f"{stem}.pangolin_out", "-c", "CHROM,POS,REF,ALT"],
            env=env, stdout=open(f"{stem}.pangolin.log", "w"), stderr=subprocess.STDOUT),
    }
    for name, p in jobs.items():
        if p.wait():
            log = f"{stem}.{'spt' if name == 'SpliceTransformer' else 'pangolin'}.log"
            raise SystemExit(f"{name} failed (exit {p.returncode}); the end of its log ({KEEP}):\n"
                             + open(log).read()[-2000:])


def main():
    ap = argparse.ArgumentParser(description=f"MetaSplice {ms.VERSION}, end to end (GRCh38; hg19 input is lifted)")
    ap.add_argument("--input", required=True, help="VCF or TSV/CSV of SNVs (chr, pos, ref, alt)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--genome", default="hg38", choices=["hg38", "hg19"])
    ap.add_argument("--data", default=os.environ.get("METASPLICE_DATA", "/data"), help="annotation bundle root")
    ap.add_argument("--cpu", action="store_true", help="run SpliceTransformer and Pangolin on the CPU (slower)")
    ap.add_argument("--keep-work", default="", help="keep intermediate files in this directory")
    ap.add_argument("--features-out", default="")
    ap.add_argument("--include-out-of-scope", action="store_true")
    ap.add_argument("--processes", type=int, default=int(os.environ.get("METASPLICE_PROCESSES", min(8, os.cpu_count() or 1))),
                    help="scoring processes, each taking whole chromosomes (default min(8, CPUs))")
    args = ap.parse_args()

    res = resources(args.data)
    need = [v for v in res if v not in OPTIONAL and not os.path.exists(res[v])]
    for idx in (".fai", ".fxi"):      # without them, pyfaidx and pyfastx would try to build them in the bundle
        if os.path.exists(res["METASPLICE_FASTA"]) and not os.path.exists(res["METASPLICE_FASTA"] + idx):
            res[f"METASPLICE_FASTA{idx}"] = res["METASPLICE_FASTA"] + idx
            need.append(f"METASPLICE_FASTA{idx}")
    if args.genome == "hg19" and not os.path.exists(res["METASPLICE_LIFTOVER"]):
        need.append("METASPLICE_LIFTOVER")
    if need:
        raise SystemExit("missing annotation resources (mount the data bundle at --data, or set the variables):\n"
                         + "\n".join(f"  {v}: {res[v]}" for v in need))

    # Absolute, because SpliceTransformer runs from its run directory inside it (review of 2026-10-14). A temporary
    # work directory is removed however the run ends (under Apptainer /tmp is the host's).
    work = os.path.abspath(args.keep_work) if args.keep_work else tempfile.mkdtemp(prefix="metasplice_")
    os.makedirs(work, exist_ok=True)
    try:
        rc = run(args, res, work)
    finally:
        if not args.keep_work:
            shutil.rmtree(work, ignore_errors=True)
    sys.exit(rc)


def run(args, res, work):
    """Select, run the two models, score; returns the scorer's exit status."""
    stem = os.path.join(work, "inscope")
    t0 = time.time()
    print("Reading the input and selecting in-scope variants"
          + (" (checking each REF against GRCh37, then lifting to GRCh38)" if args.genome == "hg19" else "") + "...",
          flush=True)
    try:
        raw, variants = in_scope_grch38(args, res)
    except UnicodeDecodeError as e:
        raise SystemExit(f"error: {args.input}: not UTF-8 text ({e.reason} at byte {e.start})")
    except ValueError as e:                     # malformed input (score_v2.parse_input, check_input)
        raise SystemExit(f"error: {e}")
    n = len(raw)
    print(f"MetaSplice {ms.VERSION}: {n:,} input SNVs; {len(variants):,} unique in-scope GRCh38 variants")
    write_inputs(variants, stem)
    if variants:
        print(f"SpliceTransformer and Pangolin on {device(args.cpu)}; TF32 off...", flush=True)
        run_tools(stem, res["METASPLICE_FASTA"], args.cpu)
    else:
        with open(f"{stem}.spt.csv", "w") as f:
            f.write(",".join(["CHROM", "POS", "REF", "ALT", "score"] + ms.TISSUES) + "\n")
        with open(f"{stem}.pangolin_out.csv", "w") as f:
            f.write("CHROM,POS,REF,ALT,Pangolin\n")
    print(f"  done in {time.time() - t0:.0f} s")

    base = [sys.executable, os.path.join(APP, "score_v2.py"), "--genome", args.genome,
            "--spt-scores", f"{stem}.spt.csv", "--pangolin-scores", f"{stem}.pangolin_out.csv"]
    for var, flag in SCORER_ARG.items():
        if var in OPTIONAL and not os.path.exists(res[var]):
            continue
        base += [flag, res[var]]
    t1 = time.time()
    if args.processes > 1 and len({v[0] for v in raw}) > 1:
        print(f"Scoring in up to {args.processes} parallel processes (by chromosome)...", flush=True)
        score_groups(args, raw, base, work)
        rc = 0
    else:
        # The scorer reads the SNVs parsed here, as in score_groups, so the skip report above is not repeated; it
        # parses them to the same (chr, pos, ref, alt) in the same order.
        write_inputs(raw, os.path.join(work, "input"), pangolin=False)
        cmd = base + ["--input", os.path.join(work, "input.vcf"), "--output", args.output]
        if args.features_out:
            cmd += ["--features-out", args.features_out]
        if args.include_out_of_scope:
            cmd += ["--include-out-of-scope"]
        rc = subprocess.call(cmd)
    print(f"  scoring done in {(time.time() - t1) / 60:.1f} min")
    return rc


if __name__ == "__main__":
    main()
