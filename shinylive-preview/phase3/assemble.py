"""Finish the pinned export with a bounded app cache and loading interface."""
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "phase2"))
from optimize import prepare, replace_once


def assemble(app, commit):
    prepare(app, commit)
    manifest_path = app / "asset-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["cache_policy"] = "at most 1 MiB per resource and 8 MiB total"

    # Use ordinary browser HTTP caching for the large runtime binaries.
    # Avoid teeing downloads into Cache Storage while webR is initializing.
    cacheable = []
    budget = 8 * 1024 * 1024
    for path, record in manifest["assets"].items():
        size = record["bytes"] + record.get("gzip_bytes", 0)
        if size <= min(budget, 1024 * 1024):
            cacheable.append(path)
            budget -= size

    worker = app / "shinylive-sw.js"
    original = worker.read_text()
    worker.write_text("const PAWEH_SMALL_ASSETS = new Set(" + json.dumps(cacheable) + ");\n"
        + replace_once(
            replace_once(original, "pawehCachePrefix + PAWEH_ASSET_REVISION;",
                         'pawehCachePrefix + "bounded-v1:" + PAWEH_ASSET_REVISION;'),
            "cache = await caches.open(pawehCacheName);",
            "if (PAWEH_SMALL_ASSETS.has(relative)) cache = await caches.open(pawehCacheName);"))

    page = app / "index.html"
    page.write_text(replace_once(
        page.read_text(), 'Starting the paweh dashboard…</h1>',
        'Preparing the R runtime…</h1>').replace(
        "</small></section></div>",
        '</small><p><button id="paweh-retry" type="button" hidden>Retry loading</button></p>'
        '<noscript>Enable JavaScript to run the dashboard.</noscript></section></div>', 1))
    (app / "shell.js").write_bytes((HERE / "shell.js").read_bytes())
    with (app / "shell.css").open("a") as css:
        css.write("\n#paweh-loading { top: var(--paweh-banner-height, 3rem); }\n"
                  "#paweh-retry { font: inherit; color: #145761; background: #fff; "
                  "padding: .45rem .65rem; border: 1px solid currentColor; border-radius: .35rem; }\n")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    assemble(Path(sys.argv[1]), sys.argv[2])
