"""Make the installed-library archive consumed by the local app wrapper."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import tarfile

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("source_commit")
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.source_commit):
        raise ValueError("The source commit must be a full Git SHA")
    (binary,) = args.directory.glob("paweh_*.tgz")
    output = args.directory / "paweh-library.tar.gz"
    with tarfile.open(binary) as package:
        files = [member for member in package.getmembers()
                 if member.name.partition("/")[0] == "paweh"]
        desc = package.extractfile("paweh/DESCRIPTION").read().decode()
        if not re.search(r"^Built: R 4\.6\.0; wasm32-unknown-emscripten;", desc, re.M):
            raise ValueError("The published artifact must target Wasm R 4.6.0")
        if any(not (member.isfile() or member.isdir()) or ".." in Path(member.name).parts
               or member.name.startswith("paweh/libs/") for member in files):
            raise ValueError("Expected an installed pure-R paweh package")
        with gzip.GzipFile(output, "wb", mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w", format=tarfile.USTAR_FORMAT) as portable:
                for entry in sorted(files, key=lambda member: member.name):
                    portable.addfile(entry, package.extractfile(entry) if entry.isfile() else None)
    provenance = {
        "paweh_source": args.source_commit,
        "wasm_binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "wrapper_archive_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "wasm_r": "4.6.0",
    }
    (args.directory / "provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")

if __name__ == "__main__":
    main()
