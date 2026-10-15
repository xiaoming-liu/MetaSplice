"""Fetch and build the MetaSplice Release 2.0 annotation data bundle, and verify it against MANIFEST.tsv.

  python fetch_and_build.py --out metasplice_v2_data [--all-from-mirror] [--with-gnomad] [--jobs 4]
  python fetch_and_build.py --out metasplice_v2_data --verify

Each file in MANIFEST.tsv is obtained by its route and then checked:

  download  verbatim from its original host                                  md5 (gzip files: or of the content)
  gunzip    from its original host, decompressed (the reference FASTAs)       md5
  index     built from the FASTA: .fai with pysam, .fxi with pyfastx          md5
  build     built from the original data by a tool in tools/                  md5 of the uncompressed content
            (chain tables from Ensembl's chain; gnomAD tables, only with --with-gnomad)
  local     copied from this bundle folder (README, licences, tools, gene-level extract)          md5
  mirror    from the MetaSplice mirror: DANN, whose host is offline, and phyloP 470-way
            as released 2022-08-28, which UCSC replaced on 2023-09-02 (README.md)                  md5

The mirror holds the whole bundle (USF Box). MIRROR_LINKS.tsv, next to MANIFEST.tsv, gives each file's direct link
and is read by default; --mirror-links names another table, and --mirror a mirror that keeps the folder layout under
one base URL. A file stored in parts (<path>.part00, .part01, ...) is reassembled before its check. With
--all-from-mirror every file is taken from the mirror, so nothing is downloaded from the original hosts or built.

Files already present that pass their check are kept, so the script can be rerun after an interruption; --only
limits a run (or --verify) to some folders. A file that
fails its check is reported and left as <name>.failed. gnomAD is optional (the scorer reports its allele frequency;
it is not a model input); --with-gnomad downloads gnomAD's joint VCFs one chromosome at a time.

Needs Python 3.8+, pysam and pyfastx (indexes), and bash, gzip and mawk (chain tables). The MetaSplice image has all
of them:

  docker run --rm --user $(id -u):$(id -g) -v "$PWD/metasplice_v2_data":/data --entrypoint python \\
      metasplice:2.0.1 /opt/metasplice/bundle/tools/fetch_and_build.py --out /data --all-from-mirror
"""
import argparse, gzip, hashlib, os, shutil, subprocess, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

TOOLS = os.path.dirname(os.path.abspath(__file__))
HERE = os.path.dirname(TOOLS)                          # the bundle folder this script belongs to
CHAIN_URL = "https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz"
CHAIN_MD5 = "2a862e70fa1afc08f1f3e7a8cc5cda14"
GNOMAD_URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr{}.vcf.bgz"


def read_manifest(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in f if line.strip()]


def md5_of(path, uncompressed=False):
    h = hashlib.md5()
    with (gzip.open(path, "rb") if uncompressed else open(path, "rb")) as f:
        for chunk in iter(lambda: f.read(1 << 23), b""):
            h.update(chunk)
    return h.hexdigest()


def passes(row, path):
    """True if the file at path is the manifest's file: same md5, or for a built file the same uncompressed content."""
    if not os.path.isfile(path) or (row["content_md5"] == "." and os.path.getsize(path) != int(row["bytes"])):
        return False
    if md5_of(path) == row["md5"]:
        return True
    return row["content_md5"] != "." and md5_of(path, uncompressed=True) == row["content_md5"]


def fetch(url, dst, tries=5):
    """Download url to dst, resuming a partial dst.part."""
    part = dst + ".part"
    os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
    for attempt in range(1, tries + 1):
        done = os.path.getsize(part) if os.path.exists(part) else 0
        req = urllib.request.Request(url, headers={"Range": f"bytes={done}-"} if done else {})
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                if done and r.status != 206:           # server ignored the range: start again
                    done = 0
                with open(part, "ab" if done else "wb") as f:
                    shutil.copyfileobj(r, f, 1 << 23)
            os.replace(part, dst)
            return dst
        except Exception as e:                         # network errors: retry with resume
            if attempt == tries:
                raise RuntimeError(f"{url}: {e}") from e
            time.sleep(10 * attempt)


def settle(row, out, path):
    """Check a file just obtained; on failure rename it to .failed and say so."""
    if passes(row, path):
        print(f"  ok      {row['path']}", flush=True)
        return True
    os.replace(path, path + ".failed")
    print(f"  FAILED  {row['path']} (kept as {row['path']}.failed)", flush=True)
    return False


class Mirror:
    """Where the mirror keeps each file: a table of direct links (path, url), or one base URL with the bundle's layout."""

    def __init__(self, base=None, links=None):
        self.base, self.links = (base or "").rstrip("/"), {}
        if links:
            text = (urllib.request.urlopen(links, timeout=120).read().decode() if links.startswith(("http://", "https://"))
                    else open(links).read())
            for line in text.splitlines()[1:]:
                if line.strip():
                    path, url = line.split("\t")[:2]
                    self.links[path] = url

    def urls(self, path):
        """The URLs of path, in order (several for a file stored in parts); None if the mirror does not have it."""
        if path in self.links:
            return [self.links[path]]
        parts = sorted(k for k in self.links if k.startswith(path + ".part"))
        if parts:
            return [self.links[k] for k in parts]
        return [f"{self.base}/{path}"] if self.base else None


def from_mirror(urls, dst):
    if len(urls) == 1:
        return fetch(urls[0], dst)
    pieces = [fetch(u, f"{dst}.p{k:02d}") for k, u in enumerate(urls)]
    with open(dst + ".part", "wb") as fo:
        for piece in pieces:
            with open(piece, "rb") as fi:
                shutil.copyfileobj(fi, fo, 1 << 23)
    os.replace(dst + ".part", dst)
    for piece in pieces:
        os.remove(piece)
    return dst


def get_download(row, out, mirror, via_mirror=False):
    dst = os.path.join(out, row["path"])
    if via_mirror or row["route"] == "mirror":
        urls = mirror.urls(row["path"])
        if not urls:
            print(f"  skipped {row['path']} (mirror only; no mirror link for it)")
            return None
        from_mirror(urls, dst)
    elif row["route"] == "download":
        fetch(row["source"], dst)
    elif row["route"] == "gunzip":
        gz = fetch(row["source"], dst + ".download.gz")
        with gzip.open(gz, "rb") as fi, open(dst + ".part", "wb") as fo:
            shutil.copyfileobj(fi, fo, 1 << 23)
        os.replace(dst + ".part", dst)
        os.remove(gz)
    elif row["route"] == "local":
        os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
        shutil.copyfile(os.path.join(HERE, row["path"]), dst)
    return settle(row, out, dst)


INDEXES = (".fai", ".fxi", ".tbi")


def stale_indexes(rows, out):
    """Index files present but older than their data file. pyfaidx rewrites such an index, which fails on a
    read-only mount; downloads in parallel, or copies that do not keep timestamps, can leave one."""
    stale = []
    for r in rows:
        p = os.path.join(out, r["path"])
        if r["path"].endswith(INDEXES) and os.path.exists(p):
            data = p.rsplit(".", 1)[0]
            if os.path.exists(data) and os.path.getmtime(p) < os.path.getmtime(data):
                stale.append(p)
    return stale


def get_index(row, out):
    dst = os.path.join(out, row["path"])
    fasta = dst.rsplit(".", 1)[0]
    if dst.endswith(".fai"):
        import pysam
        pysam.faidx(fasta)
    elif dst.endswith(".fxi"):
        import pyfastx
        pyfastx.Fasta(fasta)
    return settle(row, out, dst)


def build_chain(rows, out, work, jobs):
    chain = fetch(CHAIN_URL, os.path.join(work, "GRCh37_to_GRCh38.chain.gz"))
    if md5_of(chain) != CHAIN_MD5:
        raise RuntimeError(f"{CHAIN_URL}: md5 differs from {CHAIN_MD5}; Ensembl has changed the file")
    subprocess.run([sys.executable, os.path.join(TOOLS, "build_liftover_tables.py"), "--chain", chain,
                    "--out19", os.path.join(out, "chain/hg19tohg38_ensembl"),
                    "--out38", os.path.join(out, "chain/hg38tohg19_ensembl"), "--jobs", str(jobs)], check=True)
    return [settle(r, out, os.path.join(out, r["path"])) for r in rows]


def build_gnomad(rows, out, work):
    sys.path.insert(0, TOOLS)
    import build_gnomad_tables
    ok = []
    for r in rows:
        chrom = r["path"].split(".chr")[1].split(".")[0]
        vcf = fetch(GNOMAD_URL.format(chrom), os.path.join(work, f"gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz"))
        build_gnomad_tables.build(vcf, os.path.join(out, r["path"]))
        os.remove(vcf)
        ok.append(settle(r, out, os.path.join(out, r["path"])))
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help="bundle folder to fill (mounted at /data by the image)")
    ap.add_argument("--manifest", default=os.path.join(HERE, "MANIFEST.tsv"))
    ap.add_argument("--mirror-links", default=os.path.join(HERE, "MIRROR_LINKS.tsv"),
                    help="table of the mirror's direct links (path, url): a file or URL (default: MIRROR_LINKS.tsv)")
    ap.add_argument("--mirror", help="base URL of a mirror that keeps the bundle's folder layout (instead of links)")
    ap.add_argument("--all-from-mirror", action="store_true",
                    help="take every file from the mirror (no downloads from the original hosts, no builds)")
    ap.add_argument("--with-gnomad", action="store_true", help="also build the optional gnomAD tables")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--verify", action="store_true", help="only check the files present; fetch nothing")
    ap.add_argument("--only", default="", help="comma-separated path prefixes, e.g. dann/,conservation/ (default: all)")
    args = ap.parse_args()
    rows = read_manifest(args.manifest)
    if args.only:
        rows = [r for r in rows if r["path"].startswith(tuple(args.only.split(",")))]
    gnomad = [r for r in rows if r["path"].startswith("gnomad/")]
    if not args.with_gnomad:
        rows = [r for r in rows if r not in gnomad]

    present = {r["path"]: passes(r, os.path.join(args.out, r["path"])) for r in rows}
    if args.verify:
        bad = [p for p, ok in present.items() if not ok]
        for p in bad:
            print(f"  {'missing' if not os.path.exists(os.path.join(args.out, p)) else 'DIFFERS'}  {p}")
        print(f"{len(rows) - len(bad)} of {len(rows)} files verified" + ("" if args.with_gnomad else " (gnomAD not checked)"))
        for p in stale_indexes(rows, args.out):
            print(f"  warning: {p} is older than its data file; run `touch {p}` (pyfaidx and htslib expect indexes "
                  "newer than their data file)")
        return 1 if bad else 0

    todo = [r for r in rows if not present[r["path"]]]
    print(f"{len(rows) - len(todo)} of {len(rows)} files already verified; {len(todo)} to obtain")
    links = args.mirror_links if (args.mirror_links.startswith(("http://", "https://")) or os.path.exists(args.mirror_links)) else None
    mirror = Mirror(args.mirror, links)
    work = os.path.join(args.out, "_work")
    os.makedirs(work, exist_ok=True)
    results = []
    viam = [r for r in todo if r["route"] != "local" and args.all_from_mirror and mirror.urls(r["path"])]
    todo = [r for r in todo if r not in viam]
    simple = [r for r in todo if r["route"] in ("download", "gunzip", "local", "mirror")]
    with ThreadPoolExecutor(args.jobs) as ex:
        results += list(ex.map(lambda r: get_download(r, args.out, mirror, via_mirror=True), viam))
        results += list(ex.map(lambda r: get_download(r, args.out, mirror), simple))
    results += [get_index(r, args.out) for r in todo if r["route"] == "index"]
    chain = [r for r in todo if r["route"] == "build" and r["path"].startswith("chain/")]
    if chain:
        results += build_chain(chain, args.out, work, args.jobs)
    gn = [r for r in todo if r["route"] == "build" and r["path"].startswith("gnomad/")]
    if gn:
        results += build_gnomad(gn, args.out, work)
    shutil.rmtree(work, ignore_errors=True)
    for p in stale_indexes(rows, args.out):          # downloads end in any order: make each index newer than its data
        os.utime(p)
    failed, skipped = results.count(False), results.count(None)
    print(f"done: {results.count(True)} obtained, {failed} failed, {skipped} skipped (mirror only)")
    return 1 if failed or skipped else 0


if __name__ == "__main__":
    sys.exit(main())
