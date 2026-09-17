#!/usr/bin/env python3
"""Create a checksum manifest for a local PDB/CCD cache."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description="Hash all files in a PhosphoFill structure cache")
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir).resolve()
    output = Path(args.output).resolve()
    files = sorted(path for path in cache_dir.rglob("*") if path.is_file())
    records = []
    for path in files:
        relative = path.relative_to(cache_dir).as_posix()
        records.append(
            {
                "path": relative,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "kind": "ccd_template" if path.stem.upper() in {"TPO", "SEP", "PTR"} else "structure",
            }
        )

    payload = {
        "cache_label": "Benchmark_combined_final/benchmark_cache",
        "file_count": len(records),
        "total_bytes": sum(record["bytes"] for record in records),
        "files": records,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(records)} cache records to {output}")


if __name__ == "__main__":
    main()
