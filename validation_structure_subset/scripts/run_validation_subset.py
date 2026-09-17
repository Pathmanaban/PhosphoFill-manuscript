#!/usr/bin/env python3
"""Prepare and run the archived PhosphoFill structural validation subset."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "inputs" / "subset_cases.json"
SOURCE_DIR = ROOT / "inputs" / "source_benchmarks"
GENERATED_DIR = ROOT / "outputs" / "generated_inputs"
RESULTS_DIR = ROOT / "outputs" / "results"
ARCHIVE_DIR = ROOT / "outputs" / "archives"
WORK_DIR = ROOT / "outputs" / "work"
VERIFICATION_DIR = ROOT / "verification"
FROZEN_DIR = ROOT / "scripts" / "frozen"
CACHE_DIR = ROOT / "inputs" / "cache"

PRIMARY_TYPES = ("TPO", "SEP", "PTR")


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_tsv(path: Path) -> tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fieldnames: List[str], rows: Iterable[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def locate_source_row(rows: List[Dict[str, str]], match: Dict[str, str], case_id: str) -> Dict[str, str]:
    hits = [row for row in rows if all(str(row.get(key, "")).strip() == str(value) for key, value in match.items())]
    if len(hits) != 1:
        raise RuntimeError(f"{case_id}: expected exactly one source row, found {len(hits)} for {match}")
    return dict(hits[0])


def prepare_inputs(selected_types: set[str]) -> Dict[str, Path]:
    manifest = read_json(MANIFEST_PATH)
    prepared: Dict[str, Path] = {}
    total_cases = 0
    total_sites = 0

    for primary_type in PRIMARY_TYPES:
        if primary_type not in selected_types:
            continue
        source_path = SOURCE_DIR / f"{primary_type}_benchmark.tsv"
        fieldnames, source_rows = read_tsv(source_path)
        output_rows: List[Dict[str, str]] = []

        for case in manifest["cases"]:
            if case["primary_type"] != primary_type:
                continue
            row = locate_source_row(source_rows, case["match"], case["case_id"])
            row.update({str(key): str(value) for key, value in case.get("overrides", {}).items()})
            row["validation_case_id"] = case["case_id"]
            row["validation_case_study"] = "1" if case.get("case_study") else "0"
            row["validation_purpose"] = case["purpose"]
            output_rows.append(row)
            total_cases += 1
            total_sites += len(case["expected_ok_positions"]) + len(case["expected_no_coords_positions"])

        if not output_rows:
            continue
        out_fields = list(fieldnames)
        for extra in ("validation_case_id", "validation_case_study", "validation_purpose"):
            if extra not in out_fields:
                out_fields.append(extra)
        output_path = GENERATED_DIR / f"{primary_type}_validation_subset.tsv"
        write_tsv(output_path, out_fields, output_rows)
        prepared[primary_type] = output_path

    expected_cases = sum(1 for case in manifest["cases"] if case["primary_type"] in selected_types)
    expected_sites = sum(
        len(case["expected_ok_positions"]) + len(case["expected_no_coords_positions"])
        for case in manifest["cases"]
        if case["primary_type"] in selected_types
    )
    if total_cases != expected_cases or total_sites != expected_sites:
        raise RuntimeError(
            f"Prepared subset count mismatch: cases {total_cases}/{expected_cases}, sites {total_sites}/{expected_sites}"
        )
    return prepared


def get_versions() -> dict:
    packages = ("numpy", "Bio", "pandas", "scipy", "matplotlib", "openmm")
    versions: Dict[str, str] = {}
    for package in packages:
        try:
            module = __import__(package)
            versions[package] = str(getattr(module, "__version__", "unknown"))
        except Exception as exc:  # version capture must not prevent the benchmark
            versions[package] = f"unavailable: {type(exc).__name__}"
    return {
        "captured_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_executable": sys.executable,
        "python_version": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "packages": versions,
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def capture_frozen_files() -> dict:
    roots = [MANIFEST_PATH, SOURCE_DIR, CACHE_DIR, FROZEN_DIR, ROOT / "environment" / "requirements.txt"]
    files: List[Path] = []
    for item in roots:
        if item.is_file():
            files.append(item)
        elif item.is_dir():
            files.extend(path for path in item.rglob("*") if path.is_file())
    records = [
        {"path": relative(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(set(files))
    ]
    return {"algorithm": "SHA-256", "file_count": len(records), "files": records}


def relative(path: Path) -> str:
    return os.path.relpath(path, ROOT)


def build_command(primary_type: str, input_path: Path, resume: bool) -> List[str]:
    command = [
        sys.executable,
        relative(FROZEN_DIR / "phosphofill_benchmark_production.py"),
        relative(input_path),
        "--graft-script",
        relative(FROZEN_DIR / "graft_phospho_openmm_production.py"),
        "--report-script",
        relative(FROZEN_DIR / "phosphofill_report_v4_modes_geom.py"),
        "--output-prefix",
        relative(RESULTS_DIR / f"{primary_type}_validation_subset"),
        "--cache-dir",
        relative(CACHE_DIR),
        "--work-dir",
        relative(WORK_DIR / primary_type),
        "--prescan-contact-weight",
        "0.0",
        "--prescan-pack-weight",
        "0.2",
        "--graft-mode",
        "joint",
        "--archive-root",
        relative(ARCHIVE_DIR / primary_type),
        "--checkpoint-every",
        "1",
    ]
    if resume:
        command.append("--resume")
    return command


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare and run the self-contained archived structural validation subset."
    )
    parser.add_argument(
        "--types",
        nargs="+",
        choices=PRIMARY_TYPES,
        default=list(PRIMARY_TYPES),
        help="Residue groups to run (default: TPO SEP PTR).",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Generate the subset TSVs and provenance records without launching PhosphoFill.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from checkpointed output TSVs rather than starting a new result file.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    selected_types = set(args.types)
    for directory in (GENERATED_DIR, RESULTS_DIR, ARCHIVE_DIR, WORK_DIR, VERIFICATION_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    prepared = prepare_inputs(selected_types)
    commands = [build_command(primary_type, prepared[primary_type], args.resume) for primary_type in PRIMARY_TYPES if primary_type in prepared]

    run_record = {
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(ROOT),
        "manifest": relative(MANIFEST_PATH),
        "selected_types": [primary_type for primary_type in PRIMARY_TYPES if primary_type in selected_types],
        "prepare_only": args.prepare_only,
        "resume": args.resume,
        "commands": commands,
        "configuration": {
            "graft_mode": "joint",
            "n_poses": 3,
            "prescan_angles": "12 x 30 degrees",
            "prescan_contact_weight": 0.0,
            "prescan_pack_weight": 0.2,
            "relax_mode": "openmm_local",
            "stage_structures_archived": True,
        },
    }
    command_record_name = "prepared_commands.json" if args.prepare_only else "run_commands.json"
    version_record_name = "preparation_software_versions.json" if args.prepare_only else "software_versions.json"
    (VERIFICATION_DIR / command_record_name).write_text(json.dumps(run_record, indent=2) + "\n", encoding="utf-8")
    (VERIFICATION_DIR / version_record_name).write_text(json.dumps(get_versions(), indent=2) + "\n", encoding="utf-8")
    (VERIFICATION_DIR / "frozen_file_manifest.json").write_text(json.dumps(capture_frozen_files(), indent=2) + "\n", encoding="utf-8")

    print(f"Prepared {len(prepared)} subset TSV files under {GENERATED_DIR}")
    if args.prepare_only:
        for command in commands:
            print(" ".join(command))
        return 0

    for primary_type, command in zip((t for t in PRIMARY_TYPES if t in prepared), commands):
        print(f"\nRunning archived structural subset: {primary_type}")
        print(" ".join(command))
        completed = subprocess.run(command, cwd=ROOT)
        if completed.returncode != 0:
            print(f"{primary_type} failed with exit code {completed.returncode}", file=sys.stderr)
            return completed.returncode

    print("\nSubset run completed. Next: python scripts/verify_validation_subset.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
