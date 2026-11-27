"""Postprocess a fresh Shinylive 0.10.12 app export in place."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
MIME = {
    ".wasm": "application/wasm",
    ".so": "application/octet-stream",
    ".js": "text/javascript",
    ".mjs": "text/javascript",
    ".css": "text/css",
}

def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"Unexpected pinned export structure around {old!r}")
    return text.replace(old, new, 1)

def prepare(app, commit):
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("Provide the full paweh source commit")
    worker = app / "shinylive-sw.js"
    html = app / "index.html"
    text = worker.read_text()
    if "pawehAssetRequest" in text or "Shinylive 0.10.12" not in text:
        raise ValueError("Use a fresh Shinylive 0.10.12 export")
    assets = {}
    records = {}
    for file in sorted((app / "shinylive").rglob("*")):
        if not file.is_file():
            continue
        path = file.relative_to(app).as_posix()
        data = file.read_bytes()
        assets[path] = None
        records[path] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
        if file.suffix in MIME and len(data) >= 128 * 1024:
            compressed = gzip.compress(data, compresslevel=9, mtime=0)
            if len(compressed) < .85 * len(data):
                Path(str(file) + ".gz").write_bytes(compressed)
                assets[path] = MIME[file.suffix]
                records[path]["gzip_bytes"] = len(compressed)
    record_data = json.dumps(records, sort_keys=True, separators=(",", ":"))
    revision = hashlib.sha256(record_data.encode()).hexdigest()[:20]

    helper = (HERE / "worker.js").read_text()
    text = replace_once(text, 'if (useCaching) {',
                        'if (pawehAssetRequest(request)) {\n'
                        '    event.respondWith(pawehAsset(event));\n'
                        '    return;\n  }\n  if (useCaching) {')
    text = replace_once(text, "return key.indexOf(version + cacheName) !== 0;",
                       "return key.startsWith(pawehCachePrefix) && key !== pawehCacheName;")
    prefix = f"const PAWEH_ASSET_REVISION = {json.dumps(revision)};\nconst PAWEH_ASSETS = {json.dumps(assets)};\n"
    worker.write_text(prefix + helper + "\n" + text)
    for filename in ("shell.css", "shell.js"):
        (app / filename).write_bytes((HERE / filename).read_bytes())
    page = replace_once(html.read_text(), "</head>",
                        '<link rel="stylesheet" href="./shell.css" />\n'
                        '<script src="./shell.js" defer></script>\n</head>')
    page = replace_once(page, '<div style="height: 100vh; width: 100vw" id="root"></div>',
                        '<header id="paweh-build"><strong>paweh 0.99.0</strong> '
                        f'<a href="https://github.com/akilanthony/paweh/commit/{commit}">'
                        f'Source {commit[:7]}</a></header>\n'
                        '<div id="paweh-loading"><section aria-labelledby="paweh-status">'
                        '<h1 id="paweh-status" role="status">Starting the paweh dashboard…</h1>'
                        '<p>The R runtime and packages run in your browser. The first visit can take a few minutes.</p>'
                        '<small id="paweh-requirements">Use an up-to-date browser with service workers enabled. '
                        'The optional 3D view also needs WebGL.</small></section></div>\n'
                        '<div id="root" aria-busy="true"></div>')
    html.write_text(page)
    (app / "asset-manifest.json").write_text(json.dumps({
        "source_commit": commit, "paweh_version": "0.99.0", "shinylive_assets": "0.10.12",
        "cache_revision": revision, "assets": records
    }, indent=2) + "\n")
    print(f"Static cache generation {revision}")
    print(f"Compressed {sum(v is not None for v in assets.values())} runtime assets")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("app", type=Path)
    parser.add_argument("source_commit")
    args = parser.parse_args()
    prepare(args.app, args.source_commit)
