#!/usr/bin/env python3
"""Verify the frozen benchmark cohort and independently minimised top-3 fields."""

from __future__ import annotations

import ast
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
PRIMARY = {"TPO": "T", "SEP": "S", "PTR": "Y"}
EXPECTED = {
    "TPO": {
        "raw_rows": 524,
        "statuses": {"OK": 521, "FAILED": 1, "ERROR": 1, "NO_COORDS": 1},
        "primary_ok": 494,
        "single": 332,
        "multi": 162,
    },
    "SEP": {
        "raw_rows": 407,
        "statuses": {"OK": 403, "NO_COORDS": 3, "FAILED": 1},
        "primary_ok": 300,
        "single": 241,
        "multi": 59,
    },
    "PTR": {
        "raw_rows": 503,
        "statuses": {"OK": 500, "NO_COORDS": 2, "FAILED": 1},
        "primary_ok": 490,
        "single": 217,
        "multi": 273,
    },
}
EXPECTED_MEDIANS = {
    "TPO": [0.4225, 0.4162, 0.42395, 0.42395, 0.4158],
    "SEP": [1.2882, 0.4527, 0.5519, 0.5533, 0.4710],
    "PTR": [2.7306, 0.5226, 1.0217, 1.0195, 0.6111],
    "Combined": [1.0389, 0.4720, 0.6057, 0.60565, 0.48745],
}
MEDIAN_COLUMNS = [
    "stage0_seed_median",
    "stage1_oracle_median",
    "stage2_selected_median",
    "top1_postmin_median",
    "top3_best_postmin_median",
]
EXPECTED_PRIORS = {
    "TPO_PRIOR_OG1_P": 1.613,
    "TPO_PRIOR_CB_OG1_P_DEG": 118.8,
    "TPO_PRIOR_CG2_CB_OG1_P_DEG": -58.9,
    "SEP_PRIOR_OG_P": 1.614,
    "SEP_PRIOR_CB_OG_P_DEG": 115.4,
    "SEP_PRIOR_SHARP_MODE_DEG": 66.3,
    "SEP_PRIOR_BROAD_RANGE": [-60.0, 50.0],
    "PTR_PRIOR_OH_P": 1.609,
    "PTR_PRIOR_CZ_OH_P_DEG": 125.1,
    "PTR_PRIOR_DEFAULT_SEED_MODE_DEG": 0.0,
}


def as_float(value: str) -> float:
    return float(str(value).strip())


def site_key(row: dict[str, str]) -> tuple[str, str, str, int, str]:
    return (
        row.get("acc_id", ""),
        row.get("pdb_id", "").lower(),
        row.get("chain", ""),
        int(float(row["position"])),
        row.get("restype", "").upper(),
    )


def source_constants(path: Path) -> dict[str, Any]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    constants: dict[str, Any] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                constants[node.targets[0].id] = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                continue
    return constants


def main() -> None:
    checks: list[dict[str, Any]] = []

    def check(name: str, observed: Any, expected: Any) -> None:
        if isinstance(observed, float) or isinstance(expected, float):
            passed = math.isclose(float(observed), float(expected), rel_tol=0.0, abs_tol=1e-8)
        else:
            passed = observed == expected
        checks.append({"name": name, "passed": passed, "observed": observed, "expected": expected})

    total_primary = total_single = total_multi = 0
    rank_consistency_failures = 0
    best_rank_failures = 0
    stage3_alias_failures = 0
    phosphate_alias_failures = 0

    for restype, primary_letter in PRIMARY.items():
        path = ROOT / "outputs" / "raw" / f"{restype}_postmin.tsv"
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        status_counts = dict(Counter(row.get("status", "") for row in rows))
        primary = [
            row
            for row in rows
            if row.get("status", "").upper() == "OK" and row.get("restype", "").upper() == primary_letter
        ]
        keys = [site_key(row) for row in primary]
        single = sum(int(float(row["context_n_sites"])) == 1 for row in primary)
        multi = sum(int(float(row["context_n_sites"])) > 1 for row in primary)

        for row in primary:
            values = [as_float(row[f"rank{rank}_postmin_rmsd"]) for rank in (1, 2, 3)]
            reported_best = as_float(row["top3_best_postmin_rmsd"])
            minimum = min(values)
            if not math.isclose(reported_best, minimum, abs_tol=1e-4):
                rank_consistency_failures += 1
            reported_rank = int(float(row["top3_best_postmin_rank"]))
            tied_best = {index + 1 for index, value in enumerate(values) if math.isclose(value, minimum, abs_tol=1e-4)}
            if reported_rank not in tied_best:
                best_rank_failures += 1
            if not math.isclose(
                as_float(row["stage3_post_minimization_rmsd"]), values[0], abs_tol=1e-4
            ):
                stage3_alias_failures += 1
            if not math.isclose(as_float(row["phosphate_rmsd"]), values[0], abs_tol=1e-4):
                phosphate_alias_failures += 1

        expected = EXPECTED[restype]
        check(f"{restype}: raw row count", len(rows), expected["raw_rows"])
        check(f"{restype}: status counts", status_counts, expected["statuses"])
        check(f"{restype}: primary OK count", len(primary), expected["primary_ok"])
        check(f"{restype}: single-site count", single, expected["single"])
        check(f"{restype}: multi-site count", multi, expected["multi"])
        check(f"{restype}: duplicate primary keys", len(keys) - len(set(keys)), 0)
        total_primary += len(primary)
        total_single += single
        total_multi += multi

    check("Combined: primary OK count", total_primary, 1284)
    check("Combined: single-site count", total_single, 790)
    check("Combined: multi-site count", total_multi, 494)
    check("All ranks: top3 value matches minimum", rank_consistency_failures, 0)
    check("All ranks: reported best rank is valid", best_rank_failures, 0)
    check("All ranks: Stage 3 alias equals rank 1", stage3_alias_failures, 0)
    check("All ranks: phosphate RMSD equals rank 1", phosphate_alias_failures, 0)

    summary_path = ROOT / "outputs" / "analysis" / "benchmark_postmin_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    check("Analysis summary: single-site count", summary["single_site_n"], 790)
    check("Analysis summary: multi-site count", summary["multi_site_n"], 494)
    summary_rows = {row["type"]: row for row in summary["single_site_summary"]}
    for label, expected_values in EXPECTED_MEDIANS.items():
        for column, expected_value in zip(MEDIAN_COLUMNS, expected_values):
            check(f"{label}: {column}", float(summary_rows[label][column]), expected_value)
    check(
        "Analysis summary: all three poses independently measured",
        summary["statistics"]["all_three_poses_independently_measured"],
        True,
    )

    constants = source_constants(ROOT / "scripts" / "frozen" / "graft_phospho_openmm_production.py")
    for name, expected in EXPECTED_PRIORS.items():
        observed = constants.get(name)
        if isinstance(observed, tuple):
            observed = list(observed)
        check(f"Source prior: {name}", observed, expected)

    benchmark_source = (
        ROOT / "scripts" / "frozen" / "phosphofill_benchmark_production.py"
    ).read_text(encoding="utf-8")
    reference_index = benchmark_source.find("ref_phosphates[pos] = get_phosphate_coords(ref_res)")
    strip_index = benchmark_source.find("strip_phosphate_from_structure(stripped_structure, chain_id, positions)")
    check(
        "Source: experimental phosphate coordinates captured before stripping",
        reference_index >= 0 and strip_index >= 0 and reference_index < strip_index,
        True,
    )
    check(
        "Source: each rank RMSD is evaluated against ref_coords",
        "phosphate_rmsd(ref_coords, coords)" in benchmark_source,
        True,
    )

    failed = [item for item in checks if not item["passed"]]
    report = {
        "status": "PASS" if not failed else "FAIL",
        "checks_run": len(checks),
        "checks_failed": len(failed),
        "cohort": {
            "primary_ok": total_primary,
            "single_site": total_single,
            "multi_site": total_multi,
        },
        "checks": checks,
    }
    output = ROOT / "verification" / "benchmark_verification.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"{report['status']}: {len(checks)} checks; {len(failed)} failed")
    print(f"Report: {output}")
    if failed:
        for item in failed:
            print(f"FAILED {item['name']}: observed={item['observed']!r}, expected={item['expected']!r}")
        sys.exit(1)


if __name__ == "__main__":
    main()
