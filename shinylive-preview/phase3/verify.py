"""Verify documentation preservation and checksum the combined site."""
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote, urljoin, urlsplit


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.paths = set()

    def handle_starttag(self, tag, attrs):
        self.paths.update(v for k, v in attrs if k in ("href", "src") and v)


def checksum(path):
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()


def inventory(site):
    return {f.relative_to(site).as_posix(): f for f in sorted(site.rglob("*")) if f.is_file()}


def verify(site, previous, source, output, docs_snapshot):
    expected = set(subprocess.check_output(
        ["git", "ls-tree", "-rz", "--name-only", previous], text=True
    ).rstrip("\0").split("\0"))
    paths = inventory(site)
    required = {"index.html", "articles/index.html", "reference/index.html",
                ".nojekyll", "app/index.html", "app/app.json", "app/shinylive-sw.js"}
    missing = (expected | required) - paths.keys()
    if missing:
        raise ValueError(f"Required routes missing: {sorted(missing)}")

    articles = sorted(site.glob("articles/paweh-*.html"))
    index = (site / "articles/index.html").read_text()
    if len(articles) != 9 or any(a.name not in index for a in articles):
        raise ValueError("Expected nine linked articles")
    exports = set(re.findall(r"^export\(([^)]+)\)", (source / "NAMESPACE").read_text(), re.M))
    index = (site / "reference/index.html").read_text()
    if len(exports) != 61 or any(f"{f}()</a>" not in index for f in exports):
        raise ValueError("Expected all 61 exported functions in the reference index")

    broken = set()
    for path, file in paths.items():
        if not path.endswith(".html") or path.startswith("app/"):
            continue
        links = Links()
        links.feed(file.read_text(encoding="utf-8"))
        for value in links.paths:
            link = urlsplit(urljoin(f"https://akilanthony.github.io/paweh/{path}", value))
            if link.netloc != "akilanthony.github.io" or not link.path.startswith("/paweh/"):
                continue
            relative = unquote(link.path[len("/paweh/"):])
            if relative not in paths and str(Path(relative) / "index.html") not in paths:
                broken.add((path, relative))
    if broken:
        raise ValueError(f"Broken internal routes: {sorted(broken)}")

    docs = json.loads(docs_snapshot.read_text())
    if set(docs) != {p for p in paths if not p.startswith("app/")} or any(
        checksum(paths[path]) != digest for path, digest in docs.items()
    ):
        raise ValueError("Documentation changed after the app export")

    manifest = json.loads((site / "app/asset-manifest.json").read_text())
    for path, record in manifest["assets"].items():
        if checksum(paths[f"app/{path}"]) != record["sha256"]:
            raise ValueError(f"Runtime checksum mismatch for {path}")
    lock = json.loads(Path("shinylive-preview/phase2/lock.json").read_text())
    if manifest["source_commit"] != lock["source_commit"]:
        raise ValueError("Wrong source revision in the app")
    for package in lock["packages"]:
        artifact = paths["app/shinylive/webr/" + package["path"]]
        if checksum(artifact) != package["sha256"]:
            raise ValueError(f"Package checksum mismatch for {package['name']}")

    records = {path: {"bytes": file.stat().st_size, "sha256": checksum(file)}
               for path, file in paths.items()}
    total = sum(record["bytes"] for record in records.values())
    if total >= 10**9 or any(record["bytes"] >= 100 * 1024**2 for record in records.values()):
        raise ValueError("Combined output exceeds the Pages site or Git file size limit")
    print(f"{len(records)} files and {total:,} bytes; {len(expected)} previous paths preserved")
    output.write_text(json.dumps({
        "files": records, "bytes": total, "prior_paths": len(expected),
        "articles": len(articles), "exports": len(exports),
    }, indent=2) + "\n")


if __name__ == "__main__":
    if len(sys.argv) == 3:
        site, output = map(Path, sys.argv[1:])
        output.write_text(json.dumps({
            path: checksum(f) for path, f in inventory(site).items()
        }, indent=2) + "\n")
    else:
        site, previous, source, output, docs = sys.argv[1:]
        verify(Path(site), previous, Path(source), Path(output), Path(docs))
