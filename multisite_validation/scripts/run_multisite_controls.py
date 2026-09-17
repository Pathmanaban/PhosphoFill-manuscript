#!/usr/bin/env python3
"""Prepare and run the four PhosphoFill multi-site sensitivity controls."""

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
SOURCE_DIR = ROOT / "inputs" / "source_benchmarks"
CACHE_DIR = ROOT / "inputs" / "cache"
GENERATED_DIR = ROOT / "outputs" / "generated_inputs"
RUNS_DIR = ROOT / "outputs" / "runs"
VERIFICATION_DIR = ROOT / "verification"
FROZEN_DIR = ROOT / "scripts" / "frozen"

PRIMARY_TYPES = ("TPO", "SEP", "PTR")
CONTROL_MODES = ("n_to_c", "c_to_n", "independent", "joint_min")
EXPECTED = {
    "TPO": {"rows": 80, "sites": 169},
    "SEP": {"rows": 32, "sites": 76},
    "PTR": {"rows": 136, "sites": 274},
}
MODE_CONFIG = {
    "n_to_c": {
        "graft_mode": "joint",
        "site_order": "n_to_c",
        "minimization_coupling": "sequential",
        "description": "N-to-C context-aware scan followed by N-to-C sequential OpenMM minimization.",
    },
    "c_to_n": {
        "graft_mode": "joint",
        "site_order": "c_to_n",
        "minimization_coupling": "sequential",
        "description": "C-to-N context-aware scan followed by C-to-N sequential OpenMM minimization.",
    },
    "independent": {
        "graft_mode": "single",
        "site_order": "n_to_c",
        "minimization_coupling": "sequential",
        "description": "Each target is scanned and minimized separately on the original fully stripped structure.",
    },
    "joint_min": {
        "graft_mode": "joint",
        "site_order": "n_to_c",
        "minimization_coupling": "joint",
        "description": "N-to-C context-aware scan followed by one OpenMM minimization containing all modified sites.",
    },
}


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


def prepare_inputs(selected_types: set[str]) -> Dict[str, Path]:
    prepared: Dict[str, Path] = {}
    selection_summary: Dict[str, dict] = {}
    for primary_type in PRIMARY_TYPES:
        if primary_type not in selected_types:
            continue
        source_path = SOURCE_DIR / f"{primary_type}_benchmark.tsv"
        fields, source_rows = read_tsv(source_path)
        selected: List[Dict[str, str]] = []
        for source_index, source_row in enumerate(source_rows, start=1):
            try:
                n_target_sites = int(str(source_row.get("n_target_sites", "0")).strip() or 0)
            except ValueError:
                n_target_sites = 0
            if n_target_sites <= 1:
                continue
            row = dict(source_row)
            pdb_id = str(row.get("PDBID_1", "pdb")).strip().lower()
            chain = str(row.get("CHAINID_1", "chain")).strip()
            positions = str(row.get("newpos_1", "")).strip()
            position_tag = positions.replace(",", "_")
            row["multi_site_case_id"] = f"{primary_type}-{source_index:04d}-{pdb_id}-{chain}-{position_tag}"
            row["primary_type"] = primary_type
            row["source_benchmark_row_index"] = str(source_index)
            selected.append(row)

        observed_sites = sum(int(row["n_target_sites"]) for row in selected)
        expected = EXPECTED[primary_type]
        if len(selected) != expected["rows"] or observed_sites != expected["sites"]:
            raise RuntimeError(
                f"{primary_type} selection mismatch: rows={len(selected)}/{expected['rows']}, "
                f"sites={observed_sites}/{expected['sites']}"
            )
        out_fields = list(fields)
        for extra in ("multi_site_case_id", "primary_type", "source_benchmark_row_index"):
            if extra not in out_fields:
                out_fields.append(extra)
        out_path = GENERATED_DIR / f"{primary_type}_multi_site.tsv"
        write_tsv(out_path, out_fields, selected)
        prepared[primary_type] = out_path
        selection_summary[primary_type] = {
            "source": str(source_path.relative_to(ROOT)),
            "output": str(out_path.relative_to(ROOT)),
            "rows": len(selected),
            "requested_sites": observed_sites,
        }

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "selection_rule": "n_target_sites > 1",
        "totals": {
            "rows": sum(value["rows"] for value in selection_summary.values()),
            "requested_sites": sum(value["requested_sites"] for value in selection_summary.values()),
        },
        "by_primary_type": selection_summary,
        "control_modes": MODE_CONFIG,
    }
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    (GENERATED_DIR / "multi_site_input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return prepared


def relpath(path: Path) -> str:
    return os.path.relpath(path, ROOT)


def build_command(mode: str, primary_type: str, input_path: Path, resume: bool, save_stage_structures: bool) -> List[str]:
    config = MODE_CONFIG[mode]
    run_root = RUNS_DIR / mode
    result_dir = run_root / "results"
    archive_dir = run_root / "archives" / primary_type
    result_dir.mkdir(parents=True, exist_ok=True)
    archive_dir.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        relpath(FROZEN_DIR / "phosphofill_benchmark_controls.py"),
        relpath(input_path),
        "--graft-script",
        relpath(FROZEN_DIR / "graft_phospho_openmm_controls.py"),
        "--report-script",
        relpath(FROZEN_DIR / "phosphofill_report_v4_modes_geom.py"),
        "--output-prefix",
        relpath(result_dir / f"{primary_type}_multi_site"),
        "--cache-dir",
        relpath(CACHE_DIR),
        "--prescan-contact-weight",
        "0.0",
        "--prescan-pack-weight",
        "0.2",
        "--control-mode",
        mode,
        "--graft-mode",
        config["graft_mode"],
        "--site-order",
        config["site_order"],
        "--minimization-coupling",
        config["minimization_coupling"],
        "--archive-root",
        relpath(archive_dir),
        "--checkpoint-every",
        "1",
    ]
    if not save_stage_structures:
        command.append("--no-save-stage-pdbs")
    if resume:
        command.append("--resume")
    return command


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def capture_files() -> dict:
    roots = [SOURCE_DIR, CACHE_DIR, FROZEN_DIR, ROOT / "environment" / "requirements.txt"]
    files: List[Path] = []
    for item in roots:
        if item.is_file():
            files.append(item)
        elif item.is_dir():
            files.extend(path for path in item.rglob("*") if path.is_file())
    records = [
        {"path": relpath(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(set(files))
    ]
    return {"algorithm": "SHA-256", "file_count": len(records), "files": records}


def software_versions() -> dict:
    packages = ("numpy", "Bio", "pandas", "scipy", "matplotlib", "openmm")
    versions: Dict[str, str] = {}
    for package in packages:
        try:
            module = __import__(package)
            versions[package] = str(getattr(module, "__version__", "unknown"))
        except Exception as exc:
            versions[package] = f"unavailable: {type(exc).__name__}"
    return {
        "captured_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_executable": sys.executable,
        "python_version": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "packages": versions,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run multi-site order and coupling controls.")
    parser.add_argument("--modes", nargs="+", choices=CONTROL_MODES, default=list(CONTROL_MODES))
    parser.add_argument("--types", nargs="+", choices=PRIMARY_TYPES, default=list(PRIMARY_TYPES))
    parser.add_argument("--resume", action="store_true", help="Resume successful checkpoints and retry status=ERROR cases.")
    parser.add_argument("--prepare-only", action="store_true", help="Generate inputs and commands without running OpenMM.")
    parser.add_argument(
        "--save-stage-structures",
        action="store_true",
        help="Archive all stage CIFs. Omit for the recommended metrics-first run to avoid several GB of duplicated structures.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    selected_types = set(args.types)
    selected_modes = set(args.modes)
    for directory in (GENERATED_DIR, RUNS_DIR, VERIFICATION_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    prepared = prepare_inputs(selected_types)
    run_items = [
        (mode, primary_type, build_command(mode, primary_type, prepared[primary_type], args.resume, args.save_stage_structures))
        for mode in CONTROL_MODES
        if mode in selected_modes
        for primary_type in PRIMARY_TYPES
        if primary_type in prepared
    ]
    record = {
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "prepare_only": args.prepare_only,
        "resume": args.resume,
        "save_stage_structures": args.save_stage_structures,
        "selected_modes": [mode for mode in CONTROL_MODES if mode in selected_modes],
        "selected_types": [primary_type for primary_type in PRIMARY_TYPES if primary_type in selected_types],
        "mode_definitions": MODE_CONFIG,
        "commands": [command for _, _, command in run_items],
    }
    command_name = "prepared_commands.json" if args.prepare_only else "run_commands.json"
    version_name = "preparation_software_versions.json" if args.prepare_only else "software_versions.json"
    (VERIFICATION_DIR / command_name).write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    (VERIFICATION_DIR / version_name).write_text(json.dumps(software_versions(), indent=2) + "\n", encoding="utf-8")
    (VERIFICATION_DIR / "frozen_file_manifest.json").write_text(json.dumps(capture_files(), indent=2) + "\n", encoding="utf-8")

    print(f"Prepared {len(prepared)} multi-site input TSVs ({sum(EXPECTED[t]['rows'] for t in selected_types)} rows).")
    if args.prepare_only:
        for mode, primary_type, command in run_items:
            print(f"[{mode}/{primary_type}] {' '.join(command)}")
        return 0

    for mode, primary_type, command in run_items:
        print(f"\nRunning {mode} / {primary_type}")
        print(" ".join(command))
        completed = subprocess.run(command, cwd=ROOT)
        if completed.returncode != 0:
            print(f"{mode}/{primary_type} failed with exit code {completed.returncode}", file=sys.stderr)
            return completed.returncode

    print("\nAll requested controls completed.")
    print("Next: python scripts/analyze_multisite_controls.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
