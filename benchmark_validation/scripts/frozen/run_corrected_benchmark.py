#!/usr/bin/env python3
"""Run the corrected TPO/SEP/PTR benchmark and regenerate Figure 3/Tables 2-3."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run(command: list[str], cwd: Path) -> None:
    print("\n$ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=str(cwd), check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run all three corrected strip-and-regraft benchmarks, then create the validation figure and tables."
    )
    parser.add_argument("--outdir", default="postmin_benchmark", help="Directory for benchmark TSVs, summaries, and plots")
    parser.add_argument("--cache-dir", default="benchmark_cache", help="Existing mmCIF/CCD cache")
    parser.add_argument("--max-rows", type=int, default=None, help="Limit each residue benchmark for a smoke test")
    parser.add_argument("--archive-root", default=None, help="Optionally retain all per-case structures and reports")
    parser.add_argument("--no-resume", action="store_true", help="Do not resume from checkpointed output TSVs")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    outdir = (root / args.outdir).resolve() if not Path(args.outdir).is_absolute() else Path(args.outdir)
    results_dir = outdir / "results"
    analysis_dir = outdir / "analysis"
    results_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir.mkdir(parents=True, exist_ok=True)

    benchmark_script = root / "phosphofill_benchmark_production.py"
    graft_script = root / "graft_phospho_openmm_production.py"
    report_script = root / "phosphofill_report_v4_modes_geom.py"
    analysis_script = root / "analyse_benchmark_postmin.py"
    cache_dir = (root / args.cache_dir).resolve() if not Path(args.cache_dir).is_absolute() else Path(args.cache_dir)

    output_tsvs: list[str] = []
    for restype in ("TPO", "SEP", "PTR"):
        benchmark_tsv = root / f"{restype}_benchmark.tsv"
        output_prefix = results_dir / f"{restype}_postmin"
        command = [
            sys.executable,
            str(benchmark_script),
            str(benchmark_tsv),
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
            archive_base = Path(args.archive_root)
            if not archive_base.is_absolute():
                archive_base = (root / archive_base).resolve()
            command.extend(["--archive-root", str(archive_base / restype)])
        run(command, root)
        output_tsvs.append(str(output_prefix) + ".tsv")

    analysis_command = [
        sys.executable,
        str(analysis_script),
        "--phosphofill-tsv",
        *output_tsvs,
        "--outdir",
        str(analysis_dir),
    ]
    run(analysis_command, root)
    print(f"\nCorrected benchmark and manuscript outputs are in: {outdir}", flush=True)


if __name__ == "__main__":
    main()
