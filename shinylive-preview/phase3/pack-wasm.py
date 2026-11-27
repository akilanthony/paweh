"""Validate a pure-R artifact from the pinned Wasm toolchain and wrap it for Shinylive."""
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
import tarfile


def pack(directory, commit):
    lock = json.loads(Path("shinylive-preview/phase2/lock.json").read_text())
    if commit != lock["source_commit"]:
        raise ValueError("The build must use the locked paweh source")
    binary = directory / f"paweh_{lock['package_version']}.tgz"
    output = directory / "paweh-library.tar.gz"
    with tarfile.open(binary) as package:
        files = [m for m in package.getmembers() if m.name.partition("/")[0] == "paweh"]
        desc = package.extractfile("paweh/DESCRIPTION").read().decode()
        built = re.search(r"^Built: (.+)$", desc, re.M).group(1)
        version, platform, *_ = (x.strip() for x in built.split(";"))
        if version != f"R {lock['wasm_r']}" or platform not in ("", "wasm32-unknown-emscripten"):
            raise ValueError(f"Unexpected R version or target: {built}")
        if f"Version: {lock['package_version']}\n" not in desc or any(
            not (m.isfile() or m.isdir()) or ".." in Path(m.name).parts
            or m.name.startswith("paweh/libs/") for m in files
        ):
            raise ValueError("Expected the locked, pure-R paweh package")
        with gzip.GzipFile(output, "wb", mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w", format=tarfile.USTAR_FORMAT) as portable:
                for entry in sorted(files, key=lambda m: m.name):
                    data = package.extractfile(entry) if entry.isfile() else None
                    portable.addfile(entry, data)
    provenance = {
        "paweh_source": commit,
        "wasm_binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "wrapper_archive_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "wasm_r": lock["wasm_r"],
        "built": built,
    }
    (directory / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    pack(Path(sys.argv[1]), sys.argv[2])
