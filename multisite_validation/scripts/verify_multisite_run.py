#!/usr/bin/env python3
"""Audit completeness and internal consistency of multi-site control outputs."""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "outputs" / "generated_inputs"
RUNS_DIR = ROOT / "outputs" / "runs"
VERIFICATION_DIR = ROOT / "verification"
PRIMARY_TYPES = ("TPO", "SEP", "PTR")
MODES = ("n_to_c", "c_to_n", "independent", "joint_min")
EXPECTED_MODE = {
    "n_to_c": ("joint", "n_to_c", "sequential"),
    "c_to_n": ("joint", "c_to_n", "sequential"),
    "independent": ("single", "n_to_c", "sequential"),
    "joint_min": ("joint", "n_to_c", "joint"),
}


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value: object) -> Optional[float]:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def positions(text: object) -> List[int]:
    return [int(chunk.strip()) for chunk in str(text or "").split(",") if chunk.strip()]


def expected_sites() -> Dict[Tuple[str, int], dict]:
    expected: Dict[Tuple[str, int], dict] = {}
    for primary_type in PRIMARY_TYPES:
        for row in read_tsv(INPUT_DIR / f"{primary_type}_multi_site.tsv"):
            case_id = row["multi_site_case_id"]
            for position in positions(row["newpos_1"]):
                expected[(case_id, position)] = {
                    "primary_type": primary_type,
                    "pdb_id": row.get("PDBID_1", ""),
                    "chain": row.get("CHAINID_1", ""),
                }
    return expected


def main() -> int:
    checks: List[dict] = []

    def check(check_id: str, passed: bool, message: str, **details) -> None:
        checks.append({"check_id": check_id, "passed": bool(passed), "message": message, "details": details})

    expected = expected_sites()
    check("expected_site_count", len(expected) == 519, "Generated inputs contain 519 requested multi-site targets.", observed=len(expected))
    mode_counts: Dict[str, dict] = {}

    for mode in MODES:
        observed: Dict[Tuple[str, int], dict] = {}
        error_cases: Dict[str, dict] = {}
        status_counts: Dict[str, int] = {}
        for primary_type in PRIMARY_TYPES:
            path = RUNS_DIR / mode / "results" / f"{primary_type}_multi_site.tsv"
            check(f"result_file:{mode}:{primary_type}", path.exists(), "Control result TSV exists.", path=str(path))
            if not path.exists():
                continue
            for row in read_tsv(path):
                status = str(row.get("status", "")).strip().upper() or "MISSING"
                status_counts[status] = status_counts.get(status, 0) + 1
                case_id = str(row.get("multi_site_case_id", "")).strip()
                pos = as_float(row.get("target_position") or row.get("position"))
                if case_id and pos is not None:
                    observed[(case_id, int(pos))] = row
                elif case_id:
                    error_cases[case_id] = row

                observed_semantics = (
                    str(row.get("run_graft_mode", "")).strip(),
                    str(row.get("run_site_order", "")).strip(),
                    str(row.get("run_minimization_coupling", "")).strip(),
                )
                check(
                    f"semantics:{mode}:{primary_type}:{row.get('benchmark_case_key','row')}",
                    observed_semantics == EXPECTED_MODE[mode],
                    "Recorded execution semantics match the named control.",
                    expected=EXPECTED_MODE[mode],
                    observed=observed_semantics,
                )

        represented = set(observed)
        for key in expected:
            case_id, position = key
            row = observed.get(key)
            if row is None:
                check(
                    f"site_row:{mode}:{case_id}:{position}",
                    case_id in error_cases,
                    "Every expected site has a per-site result or an explicit case-level ERROR record.",
                )
                continue
            if str(row.get("status", "")).strip().upper() != "OK":
                continue
            ranks = [as_float(row.get(f"rank{rank}_postmin_rmsd")) for rank in (1, 2, 3)]
            reported_best = as_float(row.get("top3_best_postmin_rmsd"))
            reported_rank = as_float(row.get("top3_best_postmin_rank"))
            complete = all(value is not None for value in ranks)
            correct_min = complete and reported_best is not None and abs(min(ranks) - reported_best) <= 0.00011
            correct_rank = complete and reported_rank == ranks.index(min(ranks)) + 1
            check(
                f"top3:{mode}:{case_id}:{position}",
                complete and correct_min and correct_rank,
                "Top-1/2/3 RMSDs exist and best-of-3 is their exact minimum.",
                ranks=ranks,
                reported_best=reported_best,
                reported_rank=reported_rank,
            )

        mode_counts[mode] = {
            "expected_sites": len(expected),
            "per_site_rows": len(observed),
            "case_level_errors": len(error_cases),
            "status_counts": status_counts,
        }

    failures = [item for item in checks if not item["passed"]]
    report = {
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "summary": {
            "checks": len(checks),
            "passed": len(checks) - len(failures),
            "failed": len(failures),
            "expected_sites_per_mode": len(expected),
        },
        "mode_counts": mode_counts,
        "checks": checks,
    }
    VERIFICATION_DIR.mkdir(parents=True, exist_ok=True)
    path = VERIFICATION_DIR / "multisite_run_audit.json"
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report["summary"], indent=2))
    print(f"Detailed audit: {path}")
    if failures:
        print("Failed checks:")
        for failure in failures[:50]:
            print(f"- {failure['check_id']}: {failure['message']}")
        if len(failures) > 50:
            print(f"... and {len(failures) - 50} more")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
