#!/usr/bin/env python3
"""Portable runner for the frozen PhosphoFill strip-and-regraft benchmark.

The original scripts are retained without modification in ``scripts/frozen``.
This wrapper only adapts their paths to the portable package directory layout.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def run(command: list[str], cwd: Path) -> None:
    print("\n$ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=str(cwd), check=True)


def resolve_from_root(root: Path, value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the frozen TPO/SEP/PTR strip-and-regraft benchmark and "
            "regenerate Figure 3 and Tables 2-3."
        )
    )
    parser.add_argument(
        "--outdir",
        default="outputs/reproduced_run",
        help="Output directory relative to the Benchmark package root",
    )
    parser.add_argument(
        "--cache-dir",
        default="inputs/structure_cache",
        help="mmCIF/CCD cache; missing public files are downloaded",
    )
    parser.add_argument("--max-rows", type=int, default=None, help="Rows per residue type for a smoke test")
    parser.add_argument("--archive-root", default=None, help="Optional per-case structure archive")
    parser.add_argument("--no-resume", action="store_true", help="Start without checkpoint resume")
    args = parser.parse_args()

    package_root = Path(__file__).resolve().parents[1]
    frozen = package_root / "scripts" / "frozen"
    definitions = package_root / "inputs" / "benchmark_definitions"
    templates = package_root / "inputs" / "ccd_templates"
    outdir = resolve_from_root(package_root, args.outdir)
    cache_dir = resolve_from_root(package_root, args.cache_dir)
    results_dir = outdir / "results"
    analysis_dir = outdir / "analysis"
    results_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Seed the cache with the exact CCD templates used by the frozen run.
    for code in ("TPO", "SEP", "PTR"):
        source = templates / f"{code}.cif"
        target = cache_dir / source.name
        if not target.exists():
            shutil.copy2(source, target)

    benchmark_script = frozen / "phosphofill_benchmark_production.py"
    graft_script = frozen / "graft_phospho_openmm_production.py"
    report_script = frozen / "phosphofill_report_v4_modes_geom.py"
    analysis_script = frozen / "analyse_benchmark_postmin.py"

    output_tsvs: list[str] = []
    for restype in ("TPO", "SEP", "PTR"):
        output_prefix = results_dir / f"{restype}_postmin"
        command = [
            sys.executable,
            str(benchmark_script),
            str(definitions / f"{restype}_benchmark.tsv"),
            "--graft-script",
            str(graft_script),
            "--report-script",
            str(report_script),
            "--output-prefix",
            str(output_prefix),
            "--cache-dir",
            str(cache_dir),
            "--prescan-contact-weight",
            "0.0",
            "--prescan-pack-weight",
            "0.2",
            "--graft-mode",
            "joint",
            "--checkpoint-every",
            "1",
        ]
        if not args.no_resume:
            command.append("--resume")
        if args.max_rows is not None:
            command.extend(["--max-rows", str(args.max_rows)])
        if args.archive_root:
            archive_root = resolve_from_root(package_root, args.archive_root)
            command.extend(["--archive-root", str(archive_root / restype)])
        run(command, package_root)
        output_tsvs.append(str(output_prefix) + ".tsv")

    run(
        [
            sys.executable,
            str(analysis_script),
            "--phosphofill-tsv",
            *output_tsvs,
            "--outdir",
            str(analysis_dir),
        ],
        package_root,
    )
    print(f"\nReproduced outputs: {outdir}", flush=True)


if __name__ == "__main__":
    main()
