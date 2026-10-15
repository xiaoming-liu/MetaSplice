"""Write MIRROR_LINKS.tsv: a public direct link for every file of the bundle on the Box mirror (maintainers).

  python make_mirror_links.py [--remote box:MetaSplice_v2_data] [--out ../MIRROR_LINKS.tsv]

Uses rclone (a remote configured for USF Box). For each file of the mirror, `rclone link` creates or returns the
file's shared link ("people with the link"), https://usf.box.com/s/<id>; the direct download link is
https://usf.box.com/shared/static/<id>. Each link is checked anonymously against the local bundle: its first 64
bytes and its total size (Content-Range) must be the file's. A link that fails is left out of the table and listed.
USF Box serves Python scripts as a preview page rather than the file, so tools/*.py are browse-only on Box; they are
local files for fetch_and_build.py and never downloaded from the mirror. The exit status is 1 if any other file is left
out. A file stored in parts appears as its parts (<path>.part00, ...), which fetch_and_build.py reassembles.
MANIFEST.tsv and MIRROR_LINKS.tsv are not listed.
"""
import argparse, os, subprocess, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor

RCLONE = os.environ.get("RCLONE", os.path.expanduser("~/bin/rclone"))
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def direct(remote, path):
    shared = subprocess.run([RCLONE, "link", f"{remote}/{path}"], capture_output=True, text=True, check=True).stdout.strip()
    host, _, sid = shared.rpartition("/s/")
    if not sid:
        raise RuntimeError(f"{path}: unexpected shared link {shared}")
    return f"{host}/shared/static/{sid}"


PART = 13107200000          # split -b 12500M: the size of every part but the last


def serves(url, path, bundle):
    """True if url serves the bundle's file (or part): the same first 64 bytes and the same total size."""
    base, offset, size = path, 0, None
    if ".part" in path.rsplit("/", 1)[-1]:
        base, k = path.rsplit(".part", 1)
        whole = os.path.getsize(os.path.join(bundle, base))
        offset = int(k) * PART
        size = min(PART, whole - offset)
    local = os.path.join(bundle, base)
    size = os.path.getsize(local) if size is None else size
    with open(local, "rb") as f:
        f.seek(offset)
        want = f.read(min(64, size))
    req = urllib.request.Request(url, headers={"Range": f"bytes=0-{len(want) - 1}"})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            got, total = r.read(), r.headers.get("Content-Range", "").rsplit("/", 1)[-1]
    except Exception:
        return False
    return got == want and total == str(size)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--remote", default="box:MetaSplice_v2_data")
    ap.add_argument("--out", default=os.path.join(HERE, "MIRROR_LINKS.tsv"))
    ap.add_argument("--bundle", default="/mnt/Work0/MetaSplice_release_v2/metasplice_v2_data",
                    help="the local bundle the links are checked against")
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()
    files = subprocess.run([RCLONE, "lsf", "-R", "--files-only", args.remote], capture_output=True, text=True,
                           check=True).stdout.split()
    files = sorted(f for f in files if f not in ("MANIFEST.tsv", "MIRROR_LINKS.tsv"))
    with ThreadPoolExecutor(args.jobs) as ex:
        urls = list(ex.map(lambda f: direct(args.remote, f), files))
        ok = list(ex.map(lambda pu: serves(pu[1], pu[0], args.bundle), zip(files, urls)))
    with open(args.out, "w") as f:
        f.write("path\turl\n")
        for p, u, k in zip(files, urls, ok):
            if k:
                f.write(f"{p}\t{u}\n")
    bad = [p for p, k in zip(files, ok) if not k]
    print(f"{len(files) - len(bad)} of {len(files)} links serve their file and are written to {args.out}")
    for p in bad:
        print(f"  left out (does not serve the file): {p}")
    data = [p for p in bad if not (p.startswith("tools/") and p.endswith(".py"))]     # tools/*.py: Box previews
    if data:
        print(f"error: {len(data)} files other than tools/*.py have no working link", file=sys.stderr)
    return 1 if data else 0


if __name__ == "__main__":
    sys.exit(main())
