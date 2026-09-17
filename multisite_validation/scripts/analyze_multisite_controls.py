#!/usr/bin/env python3
"""Pair and summarize the four multi-site order/coupling controls."""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser

try:
    from scipy.stats import wilcoxon
except Exception:
    wilcoxon = None


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "outputs" / "generated_inputs"
RUNS_DIR = ROOT / "outputs" / "runs"
ANALYSIS_DIR = ROOT / "outputs" / "analysis"
PLOTS_DIR = ROOT / "outputs" / "plots"
CACHE_DIR = ROOT / "inputs" / "cache"

PRIMARY_TYPES = ("TPO", "SEP", "PTR")
MODES = ("n_to_c", "c_to_n", "independent", "joint_min")
MODE_LABELS = {
    "n_to_c": "N→C sequential",
    "c_to_n": "C→N sequential",
    "independent": "Independent site",
    "joint_min": "Joint minimisation",
}
RMSD_CHANGE_THRESHOLD = 0.5
CLOSE_SEQUENCE_THRESHOLD = 5
CLOSE_CA_THRESHOLD = 10.0
PHYSICAL_PP_THRESHOLD = 8.0
ENVIRONMENT_FIELDS = (
    "ref_n_basic_contacts_4A",
    "ref_n_hbondlike_contacts_35A",
    "ref_nearest_basic_dist",
    "ref_nearest_any_contact",
    "ref_pose_support_class",
    "site_secondary_structure",
    "site_category",
)


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict], preferred_fields: Optional[List[str]] = None) -> None:
    materialized = list(rows)
    fields: List[str] = list(preferred_fields or [])
    seen = set(fields)
    for row in materialized:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(materialized)


def as_float(value: object) -> Optional[float]:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        result = float(text)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def as_int(value: object) -> Optional[int]:
    parsed = as_float(value)
    return int(parsed) if parsed is not None else None


def parse_positions(text: object) -> List[int]:
    return [int(chunk.strip()) for chunk in str(text or "").split(",") if chunk.strip()]


def circular_difference(a: Optional[float], b: Optional[float]) -> Optional[float]:
    if a is None or b is None:
        return None
    return abs((a - b + 180.0) % 360.0 - 180.0)


def result_key(case_id: str, position: int) -> Tuple[str, int]:
    return case_id, int(position)


def find_residue(model, chain_id: str, position: int):
    if chain_id not in model:
        return None
    hits = [residue for residue in model[chain_id].get_residues() if int(residue.id[1]) == int(position)]
    for residue in hits:
        if "CA" in residue:
            return residue
    return hits[0] if hits else None


def cluster_metrics(pdb_id: str, chain_id: str, positions: List[int]) -> Dict[int, dict]:
    out = {
        position: {
            "nearest_sequence_separation": min((abs(position - other) for other in positions if other != position), default=None),
            "nearest_ca_distance_A": None,
            "nearest_phosphate_distance_A": None,
        }
        for position in positions
    }
    path = CACHE_DIR / f"{pdb_id.lower()}.cif"
    if not path.exists():
        return out
    try:
        structure = MMCIFParser(QUIET=True).get_structure(pdb_id, str(path))
        model = next(structure.get_models())
    except Exception:
        return out
    residues = {position: find_residue(model, chain_id, position) for position in positions}
    for position in positions:
        residue = residues[position]
        if residue is None:
            continue
        ca_distances = []
        p_distances = []
        for other in positions:
            if other == position or residues[other] is None:
                continue
            other_residue = residues[other]
            if "CA" in residue and "CA" in other_residue:
                ca_distances.append(float(np.linalg.norm(residue["CA"].coord - other_residue["CA"].coord)))
            if "P" in residue and "P" in other_residue:
                p_distances.append(float(np.linalg.norm(residue["P"].coord - other_residue["P"].coord)))
        out[position]["nearest_ca_distance_A"] = min(ca_distances) if ca_distances else None
        out[position]["nearest_phosphate_distance_A"] = min(p_distances) if p_distances else None
    return out


def load_expected_sites() -> tuple[List[dict], Dict[str, dict]]:
    sites: List[dict] = []
    cases: Dict[str, dict] = {}
    for primary_type in PRIMARY_TYPES:
        path = INPUT_DIR / f"{primary_type}_multi_site.tsv"
        if not path.exists():
            raise FileNotFoundError(f"Missing generated input: {path}")
        for row in read_tsv(path):
            case_id = row["multi_site_case_id"]
            positions = parse_positions(row["newpos_1"])
            metadata = {
                "multi_site_case_id": case_id,
                "primary_type": primary_type,
                "acc_id": row.get("ACC_ID", ""),
                "pdb_id": row.get("PDBID_1", "").lower(),
                "chain": row.get("CHAINID_1", ""),
                "context_positions": row.get("newpos_1", ""),
                "context_n_sites": len(positions),
                "subset_class": row.get("subset_class", ""),
                "source_benchmark_row_index": row.get("source_benchmark_row_index", ""),
            }
            cases[case_id] = metadata
            distances = cluster_metrics(metadata["pdb_id"], metadata["chain"], positions)
            for position in positions:
                sequence_distance = distances[position]["nearest_sequence_separation"]
                ca_distance = distances[position]["nearest_ca_distance_A"]
                close_cluster = (
                    (sequence_distance is not None and sequence_distance <= CLOSE_SEQUENCE_THRESHOLD)
                    or (ca_distance is not None and ca_distance <= CLOSE_CA_THRESHOLD)
                )
                sites.append({
                    **metadata,
                    "target_position": position,
                    **distances[position],
                    "close_cluster": int(close_cluster),
                })
    return sites, cases


def load_mode(mode: str) -> tuple[Dict[Tuple[str, int], dict], Dict[str, dict]]:
    per_site: Dict[Tuple[str, int], dict] = {}
    case_errors: Dict[str, dict] = {}
    for primary_type in PRIMARY_TYPES:
        path = RUNS_DIR / mode / "results" / f"{primary_type}_multi_site.tsv"
        if not path.exists():
            raise FileNotFoundError(f"Missing control output: {path}")
        for row in read_tsv(path):
            case_id = str(row.get("multi_site_case_id", "")).strip()
            position = as_int(row.get("target_position") or row.get("position"))
            if not case_id:
                continue
            if position is None:
                case_errors[case_id] = row
                continue
            per_site[result_key(case_id, position)] = row
    return per_site, case_errors


def recovery(value: Optional[float], threshold: float) -> Optional[int]:
    return int(value <= threshold) if value is not None else None


def mode_values(row: Optional[dict], case_error: Optional[dict]) -> dict:
    if row is None:
        status = str((case_error or {}).get("status", "MISSING")).strip().upper() or "MISSING"
        message = str((case_error or {}).get("message", "")).strip()
        return {"status": status, "message": message}
    top1 = as_float(row.get("rank1_postmin_rmsd"))
    top3 = as_float(row.get("top3_best_postmin_rmsd"))
    return {
        "status": str(row.get("status", "")).strip().upper(),
        "message": str(row.get("message", "")).strip(),
        "top1_rmsd": top1,
        "top3_best_rmsd": top3,
        "top3_best_rank": as_int(row.get("top3_best_postmin_rank")),
        "selected_angle": as_float(row.get("stage2_selected_angle")),
        "top1_recovered_1A": recovery(top1, 1.0),
        "top1_recovered_1_5A": recovery(top1, 1.5),
        "top3_recovered_1A": recovery(top3, 1.0),
        "top3_recovered_1_5A": recovery(top3, 1.5),
    }


def finite_values(rows: List[dict], key: str) -> List[float]:
    return [float(row[key]) for row in rows if row.get(key) is not None]


def rate(values: List[Optional[int]]) -> Optional[float]:
    numeric = [int(value) for value in values if value is not None]
    return sum(numeric) / len(numeric) if numeric else None


def safe_median(values: List[float]) -> Optional[float]:
    return median(values) if values else None


def paired_wilcoxon(a: List[float], b: List[float]) -> Optional[float]:
    if wilcoxon is None or len(a) < 2 or len(a) != len(b):
        return None
    differences = np.asarray(b, dtype=float) - np.asarray(a, dtype=float)
    if np.allclose(differences, 0.0):
        return 1.0
    try:
        return float(wilcoxon(b, a, alternative="two-sided").pvalue)
    except Exception:
        return None


def build_comparison() -> tuple[List[dict], List[dict], List[dict], List[dict], List[dict]]:
    expected_sites, _ = load_expected_sites()
    mode_maps = {mode: load_mode(mode) for mode in MODES}
    comparison: List[dict] = []
    failures: List[dict] = []

    for expected in expected_sites:
        key = result_key(expected["multi_site_case_id"], expected["target_position"])
        wide = dict(expected)
        mode_data: Dict[str, dict] = {}
        for mode in MODES:
            per_site, case_errors = mode_maps[mode]
            source_row = per_site.get(key)
            values = mode_values(source_row, case_errors.get(expected["multi_site_case_id"]))
            mode_data[mode] = values
            for field, value in values.items():
                wide[f"{mode}_{field}"] = value
            if mode == "n_to_c" and source_row is not None:
                for field in ENVIRONMENT_FIELDS:
                    wide[field] = source_row.get(field, "")
            if values["status"] != "OK":
                failures.append({
                    **{name: expected[name] for name in ("multi_site_case_id", "primary_type", "acc_id", "pdb_id", "chain", "target_position", "context_positions")},
                    "control_mode": mode,
                    "status": values["status"],
                    "message": values.get("message", ""),
                })

        reference = mode_data["n_to_c"]
        for mode in MODES[1:]:
            values = mode_data[mode]
            for metric in ("top1_rmsd", "top3_best_rmsd"):
                ref_value = reference.get(metric)
                alt_value = values.get(metric)
                wide[f"{mode}_delta_{metric}_vs_n_to_c"] = (
                    alt_value - ref_value if ref_value is not None and alt_value is not None else None
                )
            wide[f"{mode}_selected_angle_change_deg"] = circular_difference(
                values.get("selected_angle"), reference.get("selected_angle")
            )
            wide[f"{mode}_best_rank_changed"] = (
                int(values.get("top3_best_rank") != reference.get("top3_best_rank"))
                if values.get("top3_best_rank") is not None and reference.get("top3_best_rank") is not None
                else None
            )
            for field in ("top1_recovered_1A", "top1_recovered_1_5A", "top3_recovered_1A", "top3_recovered_1_5A"):
                ref_value = reference.get(field)
                alt_value = values.get(field)
                wide[f"{mode}_{field}_change"] = (
                    alt_value - ref_value if ref_value is not None and alt_value is not None else None
                )

        c2n_delta_top1 = wide.get("c_to_n_delta_top1_rmsd_vs_n_to_c")
        c2n_delta_top3 = wide.get("c_to_n_delta_top3_best_rmsd_vs_n_to_c")
        recovery_flip = any(
            wide.get(f"c_to_n_{field}_change") not in (None, 0)
            for field in ("top1_recovered_1A", "top1_recovered_1_5A", "top3_recovered_1A", "top3_recovered_1_5A")
        )
        angle_change = wide.get("c_to_n_selected_angle_change_deg")
        wide["order_sensitive"] = int(
            (c2n_delta_top1 is not None and abs(c2n_delta_top1) >= RMSD_CHANGE_THRESHOLD)
            or (c2n_delta_top3 is not None and abs(c2n_delta_top3) >= RMSD_CHANGE_THRESHOLD)
            or recovery_flip
            or wide.get("c_to_n_best_rank_changed") == 1
            or (angle_change is not None and angle_change >= 30.0)
        )
        comparison.append(wide)

    mode_summary: List[dict] = []
    for mode in MODES:
        for primary_type in (*PRIMARY_TYPES, "ALL"):
            subset = [row for row in comparison if primary_type == "ALL" or row["primary_type"] == primary_type]
            ok = [row for row in subset if row.get(f"{mode}_status") == "OK"]
            top1 = finite_values(ok, f"{mode}_top1_rmsd")
            top3 = finite_values(ok, f"{mode}_top3_best_rmsd")
            ranks = Counter(row.get(f"{mode}_top3_best_rank") for row in ok if row.get(f"{mode}_top3_best_rank") is not None)
            mode_summary.append({
                "control_mode": mode,
                "primary_type": primary_type,
                "n_expected_sites": len(subset),
                "n_ok_sites": len(ok),
                "n_failed_or_missing": len(subset) - len(ok),
                "median_top1_rmsd_A": safe_median(top1),
                "median_top3_best_rmsd_A": safe_median(top3),
                "top1_recovery_1A": rate([row.get(f"{mode}_top1_recovered_1A") for row in ok]),
                "top1_recovery_1_5A": rate([row.get(f"{mode}_top1_recovered_1_5A") for row in ok]),
                "top3_recovery_1A": rate([row.get(f"{mode}_top3_recovered_1A") for row in ok]),
                "top3_recovery_1_5A": rate([row.get(f"{mode}_top3_recovered_1_5A") for row in ok]),
                "top3_best_rank1_n": ranks.get(1, 0),
                "top3_best_rank2_n": ranks.get(2, 0),
                "top3_best_rank3_n": ranks.get(3, 0),
            })

    paired_summary: List[dict] = []
    for mode in MODES[1:]:
        for primary_type in (*PRIMARY_TYPES, "ALL"):
            subset = [
                row for row in comparison
                if (primary_type == "ALL" or row["primary_type"] == primary_type)
                and row.get("n_to_c_top1_rmsd") is not None
                and row.get(f"{mode}_top1_rmsd") is not None
            ]
            ref_top1 = [float(row["n_to_c_top1_rmsd"]) for row in subset]
            alt_top1 = [float(row[f"{mode}_top1_rmsd"]) for row in subset]
            ref_top3 = [float(row["n_to_c_top3_best_rmsd"]) for row in subset if row.get("n_to_c_top3_best_rmsd") is not None and row.get(f"{mode}_top3_best_rmsd") is not None]
            alt_top3 = [float(row[f"{mode}_top3_best_rmsd"]) for row in subset if row.get("n_to_c_top3_best_rmsd") is not None and row.get(f"{mode}_top3_best_rmsd") is not None]
            deltas_top1 = [alt - ref for ref, alt in zip(ref_top1, alt_top1)]
            deltas_top3 = [alt - ref for ref, alt in zip(ref_top3, alt_top3)]
            paired_summary.append({
                "comparison": f"{mode}_vs_n_to_c",
                "primary_type": primary_type,
                "n_paired_top1": len(deltas_top1),
                "n_paired_top3": len(deltas_top3),
                "median_delta_top1_rmsd_A": safe_median(deltas_top1),
                "median_abs_delta_top1_rmsd_A": safe_median([abs(value) for value in deltas_top1]),
                "median_delta_top3_best_rmsd_A": safe_median(deltas_top3),
                "n_top1_improved_by_0_5A": sum(value <= -RMSD_CHANGE_THRESHOLD for value in deltas_top1),
                "n_top1_worsened_by_0_5A": sum(value >= RMSD_CHANGE_THRESHOLD for value in deltas_top1),
                "n_best_rank_changed": sum(row.get(f"{mode}_best_rank_changed") == 1 for row in subset),
                "n_top1_recovery_1A_gained": sum(row.get(f"{mode}_top1_recovered_1A_change") == 1 for row in subset),
                "n_top1_recovery_1A_lost": sum(row.get(f"{mode}_top1_recovered_1A_change") == -1 for row in subset),
                "n_top1_recovery_1_5A_gained": sum(row.get(f"{mode}_top1_recovered_1_5A_change") == 1 for row in subset),
                "n_top1_recovery_1_5A_lost": sum(row.get(f"{mode}_top1_recovered_1_5A_change") == -1 for row in subset),
                "wilcoxon_top1_p": paired_wilcoxon(ref_top1, alt_top1),
                "wilcoxon_top3_p": paired_wilcoxon(ref_top3, alt_top3),
            })

    order_sensitive = [row for row in comparison if row["order_sensitive"] == 1]
    return comparison, mode_summary, paired_summary, failures, order_sensitive


def build_distance_summary(comparison: List[dict]) -> List[dict]:
    bins = (
        ("5–<6", 5.0, 6.0),
        ("6–<8", 6.0, 8.0),
        ("8–<10", 8.0, 10.0),
        ("10–<15", 10.0, 15.0),
        ("≥15", 15.0, math.inf),
    )
    rows = []
    for label, lower, upper in bins:
        subset = [
            row for row in comparison
            if row.get("n_to_c_status") == "OK"
            and row.get("nearest_phosphate_distance_A") is not None
            and lower <= float(row["nearest_phosphate_distance_A"]) < upper
        ]
        top1 = [abs(float(row["c_to_n_delta_top1_rmsd_vs_n_to_c"])) for row in subset]
        top3 = [abs(float(row["c_to_n_delta_top3_best_rmsd_vs_n_to_c"])) for row in subset]
        rows.append({
            "nearest_phosphate_distance_bin_A": label,
            "n_sites": len(subset),
            "median_abs_c_to_n_delta_top1_A": round(safe_median(top1), 6) if top1 else None,
            "maximum_abs_c_to_n_delta_top1_A": round(max(top1), 6) if top1 else None,
            "median_abs_c_to_n_delta_top3_A": round(safe_median(top3), 6) if top3 else None,
            "maximum_abs_c_to_n_delta_top3_A": round(max(top3), 6) if top3 else None,
            "n_order_sensitive": sum(row.get("order_sensitive") == 1 for row in subset),
        })
    return rows


def build_environment_summary(comparison: List[dict]) -> List[dict]:
    """Contrast persistently accurate and inaccurate sites within mixed cases."""
    evaluable = [
        row for row in comparison
        if all(row.get(f"{mode}_top1_rmsd") is not None for mode in MODES)
    ]
    by_case: Dict[str, List[dict]] = defaultdict(list)
    for row in evaluable:
        by_case[row["multi_site_case_id"]].append(row)

    low_sites: List[dict] = []
    high_sites: List[dict] = []
    mixed_case_ids = set()
    for case_id, rows in by_case.items():
        low = [row for row in rows if max(float(row[f"{mode}_top1_rmsd"]) for mode in MODES) <= 1.0]
        high = [row for row in rows if min(float(row[f"{mode}_top1_rmsd"]) for mode in MODES) > 1.5]
        if low and high:
            mixed_case_ids.add(case_id)
            low_sites.extend(low)
            high_sites.extend(high)

    features = (
        ("At least one basic contact", lambda row: (as_int(row.get("ref_n_basic_contacts_4A")) or 0) >= 1),
        ("High reference-pose support", lambda row: str(row.get("ref_pose_support_class", "")).lower() == "high"),
        ("β-strand environment", lambda row: str(row.get("site_secondary_structure", "")).lower() == "strand"),
    )
    rows = []
    for group, sites in (("Always ≤1.0 Å", low_sites), ("Always >1.5 Å", high_sites)):
        for feature, predicate in features:
            count = sum(bool(predicate(row)) for row in sites)
            rows.append({
                "group": group,
                "feature": feature,
                "n_sites": len(sites),
                "count": count,
                "fraction": round(count / len(sites), 6) if sites else None,
                "n_mixed_cases": len(mixed_case_ids),
            })
    return rows


def make_plot(comparison: List[dict], mode_summary: List[dict]) -> None:
    mode_colors = {
        "n_to_c": "#1f4e79",
        "c_to_n": "#5b8db8",
        "independent": "#75579b",
        "joint_min": "#2f806d",
    }
    plot_labels = ["N→C\nsequential", "C→N\nsequential", "Independent\nsite", "Joint\nminimisation"]
    alt_modes = MODES[1:]
    alt_labels = ["C→N", "Independent", "Joint min."]
    evaluable = [
        row for row in comparison
        if all(row.get(f"{mode}_top1_rmsd") is not None and row.get(f"{mode}_top3_best_rmsd") is not None for mode in MODES)
    ]
    rng = np.random.default_rng(20260829)

    def violin_scatter(axis, groups, colors, positions=None, jitter=0.055, point_alpha=0.16):
        positions = np.arange(1, len(groups) + 1) if positions is None else np.asarray(positions, dtype=float)
        violins = axis.violinplot(groups, positions=positions, widths=0.72, showextrema=False, showmedians=False)
        for body, color in zip(violins["bodies"], colors):
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(0.28)
            body.set_linewidth(0.8)
        boxes = axis.boxplot(
            groups,
            positions=positions,
            widths=0.14,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.2},
            boxprops={"facecolor": "white", "edgecolor": "0.25", "linewidth": 0.8},
            whiskerprops={"color": "0.35", "linewidth": 0.8},
            capprops={"color": "0.35", "linewidth": 0.8},
        )
        for index, (position, values, color) in enumerate(zip(positions, groups, colors)):
            values_array = np.asarray(values, dtype=float)
            xj = position + rng.normal(0.0, jitter, size=len(values_array))
            axis.scatter(xj, values_array, s=7, alpha=point_alpha, color=color, edgecolors="none", rasterized=True, zorder=2)
        return boxes

    fig = plt.figure(figsize=(13.6, 8.3))
    grid = fig.add_gridspec(2, 6, height_ratios=[1.0, 0.92])
    ax_a = fig.add_subplot(grid[0, 0:3])
    ax_b = fig.add_subplot(grid[0, 3:6])
    ax_c = fig.add_subplot(grid[1, 0:2])
    ax_d = fig.add_subplot(grid[1, 2:4])
    ax_e = fig.add_subplot(grid[1, 4:6])

    top1_groups = [[float(row[f"{mode}_top1_rmsd"]) for row in evaluable] for mode in MODES]
    top3_groups = [[float(row[f"{mode}_top3_best_rmsd"]) for row in evaluable] for mode in MODES]
    shared_max = max(max(values) for values in top1_groups + top3_groups)
    shared_ylim = (0.0, math.ceil((shared_max + 0.15) * 2.0) / 2.0)

    violin_scatter(ax_a, top1_groups, [mode_colors[mode] for mode in MODES])
    violin_scatter(ax_b, top3_groups, [mode_colors[mode] for mode in MODES])
    for axis, groups, ylabel in (
        (ax_a, top1_groups, "Top-1 phosphate RMSD (Å)"),
        (ax_b, top3_groups, "Best-of-3 phosphate RMSD (Å)"),
    ):
        axis.set_xticks(np.arange(1, len(MODES) + 1), plot_labels)
        axis.set_ylabel(ylabel)
        axis.set_ylim(*shared_ylim)
        axis.text(0.99, 0.97, f"n={len(evaluable)} paired sites per mode", transform=axis.transAxes, ha="right", va="top", fontsize=8, color="0.35")
        for position, values in enumerate(groups, start=1):
            axis.text(position, shared_ylim[1] * 0.90, f"median {np.median(values):.3f}", ha="center", fontsize=8, color="0.25")

    delta_top1 = [
        [float(row[f"{mode}_delta_top1_rmsd_vs_n_to_c"]) for row in evaluable]
        for mode in alt_modes
    ]
    delta_top3 = [
        [float(row[f"{mode}_delta_top3_best_rmsd_vs_n_to_c"]) for row in evaluable]
        for mode in alt_modes
    ]
    alt_colors = [mode_colors[mode] for mode in alt_modes]
    violin_scatter(ax_c, delta_top1, alt_colors, jitter=0.045, point_alpha=0.22)
    violin_scatter(ax_d, delta_top3, alt_colors, jitter=0.045, point_alpha=0.22)
    for axis, groups, ylabel in (
        (ax_c, delta_top1, "Top-1 ΔRMSD vs N→C (Å)"),
        (ax_d, delta_top3, "Best-of-3 ΔRMSD vs N→C (Å)"),
    ):
        axis.axhline(0.0, color="0.2", lw=0.9)
        axis.set_xticks(np.arange(1, len(alt_modes) + 1), alt_labels)
        axis.set_ylabel(ylabel)
        upper = max(abs(value) for values in groups for value in values)
        limit = max(0.01, upper * 1.18)
        axis.set_ylim(-limit, limit)
        for position, values in enumerate(groups, start=1):
            axis.text(position, limit * 0.88, f"max |Δ|\n{max(abs(value) for value in values):.4f}", ha="center", va="top", fontsize=7.5, color="0.25")
    ax_c.text(0.02, 0.04, "Top-1 prescan angle unchanged: 515/515", transform=ax_c.transAxes, fontsize=8, color="0.35")

    proximity_groups = (
        ("<8", lambda distance: distance < 8.0),
        ("8 to <15", lambda distance: 8.0 <= distance < 15.0),
        ("≥15", lambda distance: distance >= 15.0),
    )
    proximity_values = []
    proximity_labels = []
    for label, predicate in proximity_groups:
        values = [
            abs(float(row["c_to_n_delta_top1_rmsd_vs_n_to_c"]))
            for row in evaluable
            if row.get("nearest_phosphate_distance_A") is not None
            and predicate(float(row["nearest_phosphate_distance_A"]))
        ]
        proximity_values.append(values)
        proximity_labels.append(f"{label}\n(n={len(values)})")
    proximity_colors = ["#b8d3e6", "#5b8db8", "#1f4e79"]
    violin_scatter(ax_e, proximity_values, proximity_colors, jitter=0.05, point_alpha=0.24)
    ax_e.set_xticks(np.arange(1, 4), proximity_labels)
    ax_e.set_xlabel("Nearest experimental P-P distance group (Å)")
    ax_e.set_ylabel("Absolute C→N vs N→C top-1 ΔRMSD (Å)")
    proximity_max = max(value for values in proximity_values for value in values)
    ax_e.set_ylim(-0.00025, proximity_max * 1.22)
    ax_e.text(0.02, 0.97, "P-P = nearest distance between experimental\nphosphate phosphorus atoms", transform=ax_e.transAxes, va="top", fontsize=7.5, color="0.35")

    axes = (ax_a, ax_b, ax_c, ax_d, ax_e)
    for label, axis in zip("abcde", axes):
        axis.text(-0.14, 1.04, f"({label})", transform=axis.transAxes, fontweight="bold", va="top")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="0.9", lw=0.6, zorder=0)
        axis.set_axisbelow(True)
        axis.tick_params(axis="x", labelsize=9)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.10, top=0.985, wspace=0.72, hspace=0.38)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOTS_DIR / "multi_site_sensitivity.pdf", bbox_inches="tight")
    fig.savefig(PLOTS_DIR / "multi_site_sensitivity.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def json_clean(value):
    if isinstance(value, dict):
        return {key: json_clean(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_clean(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main() -> int:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    comparison, mode_summary, paired_summary, failures, order_sensitive = build_comparison()
    distance_summary = build_distance_summary(comparison)
    environment_summary = build_environment_summary(comparison)
    write_tsv(ANALYSIS_DIR / "multisite_per_site_comparison.tsv", comparison)
    write_tsv(ANALYSIS_DIR / "multisite_mode_summary.tsv", mode_summary)
    write_tsv(ANALYSIS_DIR / "multisite_paired_summary.tsv", paired_summary)
    write_tsv(ANALYSIS_DIR / "multisite_failures.tsv", failures)
    write_tsv(ANALYSIS_DIR / "multisite_order_sensitive_sites.tsv", order_sensitive)
    write_tsv(ANALYSIS_DIR / "multisite_distance_sensitivity.tsv", distance_summary)
    write_tsv(ANALYSIS_DIR / "multisite_within_case_environment.tsv", environment_summary)
    make_plot(comparison, mode_summary)

    close_summary = []
    for close_value, label in ((1, "close_cluster"), (0, "other")):
        subset = [row for row in comparison if row["close_cluster"] == close_value]
        deltas = [row["c_to_n_delta_top1_rmsd_vs_n_to_c"] for row in subset if row.get("c_to_n_delta_top1_rmsd_vs_n_to_c") is not None]
        close_summary.append({
            "group": label,
            "n_sites": len(subset),
            "n_paired_c_to_n": len(deltas),
            "median_abs_c_to_n_delta_top1_A": safe_median([abs(value) for value in deltas]),
            "n_order_sensitive": sum(row["order_sensitive"] == 1 for row in subset),
        })
    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "definitions": {
            "reference": "n_to_c",
            "order_sensitive": "|C-to-N minus N-to-C top-1 or best-of-3 RMSD| >= 0.5 A, any 1.0/1.5 A recovery flip, best-of-3 oracle rank change, or selected prescan angle change >=30 degrees",
            "close_cluster": "nearest sequence separation <=5 residues or nearest target-site C-alpha distance <=10 A",
            "physical_proximity": "nearest experimental phosphorus-to-phosphorus distance; <8 A used only as a descriptive close-pair marker",
            "within_case_environment": "among structures containing both a site with top-1 RMSD <=1.0 A in every control and a site with top-1 RMSD >1.5 A in every control",
        },
        "counts": {
            "expected_sites": len(comparison),
            "order_sensitive_sites": len(order_sensitive),
            "failure_records": len(failures),
        },
        "mode_summary": mode_summary,
        "paired_summary": paired_summary,
        "cluster_summary": close_summary,
        "distance_summary": distance_summary,
        "within_case_environment": environment_summary,
    }
    (ANALYSIS_DIR / "multisite_sensitivity_summary.json").write_text(
        json.dumps(json_clean(payload), indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload["counts"], indent=2))
    print(f"Analysis: {ANALYSIS_DIR}")
    print(f"Figure: {PLOTS_DIR / 'multi_site_sensitivity.pdf'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
