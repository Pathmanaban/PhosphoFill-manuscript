#!/usr/bin/env python3
"""Numerically verify every archived stage in the structural validation subset."""

from __future__ import annotations

import csv
import itertools
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
from Bio.PDB import MMCIFParser


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "inputs" / "subset_cases.json"
RESULTS_DIR = ROOT / "outputs" / "results"
ARCHIVE_DIR = ROOT / "outputs" / "archives"
VERIFICATION_DIR = ROOT / "verification"
VISUAL_DIR = ROOT / "visual_checks"

PHOSPHATE_ATOMS = ("P", "O1P", "O2P", "O3P")
RESIDUE_LETTER = {"TPO": "T", "SEP": "S", "PTR": "Y"}
PRIORS = {
    "TPO": {"anchor": "OG1", "base": "CB", "length": 1.613, "angle": 118.8, "dihedral_atoms": ("CG2", "CB", "OG1", "P"), "dihedral": -58.9},
    "SEP": {"anchor": "OG", "base": "CB", "length": 1.614, "angle": 115.4, "dihedral_atoms": ("CA", "CB", "OG", "P"), "dihedral": 66.3},
    "PTR": {"anchor": "OH", "base": "CZ", "length": 1.609, "angle": 125.1, "dihedral_atoms": None, "dihedral": None},
}


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: object) -> Optional[float]:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def norm_pdb(value: object) -> str:
    return str(value or "").strip().upper()


def parse_structure(path: Path):
    return MMCIFParser(QUIET=True).get_structure(path.stem, str(path))


def find_residue(structure, chain_id: str, position: int, prefer_phosphate: bool = False):
    model = next(structure.get_models())
    if chain_id not in model:
        return None
    hits = [residue for residue in model[chain_id].get_residues() if int(residue.id[1]) == int(position)]
    if prefer_phosphate:
        for residue in hits:
            if "P" in residue:
                return residue
    return hits[0] if hits else None


def residue_coords(residue, atom_names: Iterable[str]) -> Optional[Dict[str, np.ndarray]]:
    if residue is None:
        return None
    coords: Dict[str, np.ndarray] = {}
    for name in atom_names:
        if name not in residue:
            return None
        coords[name] = np.asarray(residue[name].coord, dtype=float)
    return coords


def phosphate_coords(residue) -> Optional[Dict[str, np.ndarray]]:
    return residue_coords(residue, PHOSPHATE_ATOMS)


def phosphate_rmsd(ref: Dict[str, Iterable[float]], pred: Dict[str, np.ndarray]) -> float:
    ref_p = np.asarray(ref["P"], dtype=float)
    pred_p = np.asarray(pred["P"], dtype=float)
    ref_ox = [np.asarray(ref[name], dtype=float) for name in PHOSPHATE_ATOMS[1:]]
    pred_ox = [np.asarray(pred[name], dtype=float) for name in PHOSPHATE_ATOMS[1:]]
    best = math.inf
    for order in itertools.permutations(range(3)):
        ref_array = np.asarray([ref_p, *ref_ox], dtype=float)
        pred_array = np.asarray([pred_p, *(pred_ox[index] for index in order)], dtype=float)
        value = float(np.sqrt(np.mean(np.sum((ref_array - pred_array) ** 2, axis=1))))
        best = min(best, value)
    return best


def angle_deg(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    v1 = a - b
    v2 = c - b
    cosine = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def dihedral_deg(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    n1 = np.cross(b0, b1)
    n2 = np.cross(b1, b2)
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    b1 /= np.linalg.norm(b1)
    m1 = np.cross(n1, b1)
    return math.degrees(math.atan2(float(np.dot(m1, n2)), float(np.dot(n1, n2))))


def circular_difference(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


class Audit:
    def __init__(self) -> None:
        self.checks: List[dict] = []
        self.site_records: List[dict] = []

    def check(self, check_id: str, passed: bool, message: str, **details) -> None:
        self.checks.append({"check_id": check_id, "passed": bool(passed), "message": message, "details": details})

    @property
    def failures(self) -> List[dict]:
        return [item for item in self.checks if not item["passed"]]


def archive_map(audit: Audit) -> Dict[str, Path]:
    found: Dict[str, Path] = {}
    for meta_path in ARCHIVE_DIR.rglob("benchmark_row.json"):
        try:
            row = read_json(meta_path)
        except OSError as exc:
            audit.check(
                f"archive_metadata_read:{meta_path.parent.parent.name}",
                False,
                "Archive metadata could not be opened.",
                path=str(meta_path),
                error=str(exc),
            )
            continue
        case_id = str(row.get("validation_case_id", "")).strip()
        if not case_id:
            continue
        case_dir = meta_path.parents[1]
        if case_id in found and found[case_id] != case_dir:
            audit.check(f"archive_unique:{case_id}", False, "More than one archive directory carries this case ID.", paths=[str(found[case_id]), str(case_dir)])
        found[case_id] = case_dir
    return found


def result_map(audit: Audit) -> Dict[str, List[Dict[str, str]]]:
    result_rows: Dict[str, List[Dict[str, str]]] = {}
    for primary_type in PRIORS:
        path = RESULTS_DIR / f"{primary_type}_validation_subset.tsv"
        audit.check(f"result_file:{primary_type}", path.exists(), "Result TSV exists.", path=str(path))
        if not path.exists():
            continue
        result_rows[primary_type] = read_tsv(path)
    return result_rows


def matching_result_rows(case: dict, rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    match = case["match"]
    positions = set(case["expected_ok_positions"] + case["expected_no_coords_positions"])
    return [
        row
        for row in rows
        if str(row.get("acc_id", "")).strip() == match["ACC_ID"]
        and norm_pdb(row.get("pdb_id")) == norm_pdb(match["PDBID_1"])
        and str(row.get("chain", "")).strip() == match["CHAINID_1"]
        and int(float(row.get("target_position") or row.get("position"))) in positions
    ]


def check_scan_archive(audit: Audit, case: dict, case_dir: Path) -> None:
    scan_path = case_dir / "stage3" / "prescan_frames_labeled.json"
    audit.check(f"scan_file:{case['case_id']}", scan_path.exists(), "Labeled 12-angle prescan archive exists.", path=str(scan_path))
    if not scan_path.exists():
        return
    data = read_json(scan_path)
    chain = case["match"]["CHAINID_1"]
    for position in case["expected_ok_positions"] + case["expected_no_coords_positions"]:
        site = data.get(f"{chain}:{position}", {})
        angles = sorted(round(float(frame.get("angle", -999)), 6) for frame in site.get("frames", []))
        expected = [float(value) for value in range(0, 360, 30)]
        audit.check(
            f"scan_12x30:{case['case_id']}:{position}",
            angles == expected,
            "Prescan contains exactly 12 orientations at 30-degree increments.",
            observed_angles=angles,
        )


def check_stage0_geometry(audit: Audit, case: dict, case_dir: Path) -> None:
    stage0_path = case_dir / "stage0" / "structure_s0.cif"
    if not stage0_path.exists():
        return
    structure = parse_structure(stage0_path)
    prior = PRIORS[case["primary_type"]]
    chain = case["match"]["CHAINID_1"]
    for position in case["expected_ok_positions"] + case["expected_no_coords_positions"]:
        residue = find_residue(structure, chain, position, prefer_phosphate=True)
        needed = [prior["base"], prior["anchor"], "P"]
        if prior["dihedral_atoms"]:
            needed.extend(prior["dihedral_atoms"])
        names = list(dict.fromkeys(needed))
        coords = residue_coords(residue, names)
        if coords is None:
            audit.check(f"seed_geometry:{case['case_id']}:{position}", False, "Stage 0 lacks atoms required to verify the empirical geometry prior.", atoms=names)
            continue
        length = float(np.linalg.norm(coords[prior["anchor"]] - coords["P"]))
        valence = angle_deg(coords[prior["base"]], coords[prior["anchor"]], coords["P"])
        length_ok = abs(length - prior["length"]) <= 0.015
        angle_ok = abs(valence - prior["angle"]) <= 1.0
        torsion = None
        torsion_ok = True
        if prior["dihedral_atoms"]:
            torsion = dihedral_deg(*(coords[name] for name in prior["dihedral_atoms"]))
            torsion_ok = circular_difference(torsion, prior["dihedral"]) <= 1.5
        audit.check(
            f"seed_geometry:{case['case_id']}:{position}",
            length_ok and angle_ok and torsion_ok,
            "Stage 0 reproduces the configured residue-specific internal-coordinate prior.",
            observed_length=round(length, 4),
            expected_length=prior["length"],
            observed_angle=round(valence, 3),
            expected_angle=prior["angle"],
            observed_dihedral=None if torsion is None else round(torsion, 3),
            expected_dihedral=prior["dihedral"],
        )


def check_rank_archives(audit: Audit, case: dict, case_dir: Path, rows_by_position: Dict[int, Dict[str, str]]) -> None:
    top_meta_path = case_dir / "stage3" / "top3_ranked_poses.json"
    audit.check(f"top3_metadata:{case['case_id']}", top_meta_path.exists(), "Top-3 rank metadata archive exists.", path=str(top_meta_path))
    top_meta = read_json(top_meta_path) if top_meta_path.exists() else {}
    definition = str(top_meta.get("definition", ""))
    audit.check(
        f"independent_definition:{case['case_id']}",
        "independently minimized" in definition,
        "Archive explicitly defines the three outputs as independently minimized rank-matched prescan candidates.",
        definition=definition,
    )

    ref_path = case_dir / "ref" / "reference_sites.json"
    refs = read_json(ref_path) if ref_path.exists() else {}
    chain = case["match"]["CHAINID_1"]
    structures: Dict[int, object] = {}

    for rank in (1, 2, 3):
        structure_path = case_dir / "stage3" / f"phosphoFill{rank}.cif"
        report_path = case_dir / "stage3" / f"graft_report_rank{rank}.tsv"
        audit.check(f"rank_structure:{case['case_id']}:{rank}", structure_path.exists(), "Independently minimized rank structure exists.", path=str(structure_path))
        audit.check(f"rank_report:{case['case_id']}:{rank}", report_path.exists(), "Rank-specific OpenMM graft report exists.", path=str(report_path))
        if structure_path.exists():
            structures[rank] = parse_structure(structure_path)
        if report_path.exists():
            report_rows = read_tsv(report_path)
            relevant = [row for row in report_rows if str(row.get("chain_id", "")).strip() == chain]
            relax_ok = bool(relevant) and all(str(row.get("relax_status", "")).strip().upper() == "OK" for row in relevant)
            backend_ok = bool(relevant) and all(str(row.get("relax_backend", "")).strip() == "openmm_local" for row in relevant)
            audit.check(
                f"rank_openmm:{case['case_id']}:{rank}",
                relax_ok and backend_ok,
                "Rank-specific report confirms successful OpenMM minimization.",
                relax_status=[row.get("relax_status") for row in relevant],
                relax_backend=[row.get("relax_backend") for row in relevant],
            )
        rank_meta = top_meta.get("structures", {}).get(f"phosphoFill{rank}", {})
        audit.check(
            f"rank_postmin_flag:{case['case_id']}:{rank}",
            rank_meta.get("post_minimization") is True,
            "Rank metadata marks the structure as post-minimization.",
        )

    for position in case["expected_ok_positions"]:
        row = rows_by_position.get(position)
        if row is None:
            audit.check(
                f"rank_result_row:{case['case_id']}:{position}",
                False,
                "Cannot validate rank RMSDs because the expected result row is absent.",
            )
            continue
        recomputed: List[float] = []
        ref = refs.get(str(position), {}).get("coords")
        for rank in (1, 2, 3):
            pred_residue = find_residue(structures.get(rank), chain, position, prefer_phosphate=True) if rank in structures else None
            pred = phosphate_coords(pred_residue)
            if ref is None or pred is None:
                audit.check(f"rmsd_recompute:{case['case_id']}:{position}:{rank}", False, "Reference or predicted phosphate coordinates are unavailable.")
                continue
            observed = phosphate_rmsd(ref, pred)
            reported = number(row.get(f"rank{rank}_postmin_rmsd"))
            recomputed.append(observed)
            audit.check(
                f"rmsd_recompute:{case['case_id']}:{position}:{rank}",
                reported is not None and abs(observed - reported) <= 0.0011,
                "Permutation-invariant phosphate RMSD recomputes from the archived rank structure and experimental coordinates.",
                recomputed=round(observed, 5),
                reported=reported,
            )
        reported_best = number(row.get("top3_best_postmin_rmsd"))
        reported_rank = number(row.get("top3_best_postmin_rank"))
        if len(recomputed) == 3:
            best = min(recomputed)
            best_rank = recomputed.index(best) + 1
            audit.check(
                f"top3_minimum:{case['case_id']}:{position}",
                reported_best is not None and abs(best - reported_best) <= 0.0011 and reported_rank == best_rank,
                "top3_best_postmin_rmsd is the minimum of the three independently minimized rank RMSDs.",
                recomputed_best=round(best, 5),
                recomputed_rank=best_rank,
                reported_best=reported_best,
                reported_rank=reported_rank,
            )


def write_pymol_loader(case_dirs: Dict[str, Path], manifest: dict) -> None:
    lines = [
        "# Generated by verify_validation_subset.py",
        "reinitialize",
        "set retain_order, 1",
        "set dash_width, 2.5",
        "bg_color white",
    ]
    for case in manifest["cases"]:
        case_id = case["case_id"]
        case_dir = case_dirs.get(case_id)
        if case_dir is None:
            continue
        safe = case_id.replace("-", "_")
        for label, rel_path in (
            ("ref", Path("input/reference_phospho.cif")),
            ("stripped", Path("input/stripped_input.cif")),
            ("s0", Path("stage0/structure_s0.cif")),
            ("s1", Path("stage1/structure_s1.cif")),
            ("s2", Path("stage2/structure_s2.cif")),
            ("top1", Path("stage3/phosphoFill1.cif")),
            ("top2", Path("stage3/phosphoFill2.cif")),
            ("top3", Path("stage3/phosphoFill3.cif")),
        ):
            path = case_dir / rel_path
            if path.exists():
                lines.append(f'load "{path.as_posix()}", {safe}_{label}')
        positions = "+".join(str(value) for value in case["expected_ok_positions"] + case["expected_no_coords_positions"])
        lines.extend([
            f"hide everything, {safe}_*",
            f"show sticks, {safe}_* and chain {case['match']['CHAINID_1']} and resi {positions}",
            f"color gray70, {safe}_ref",
            f"color cyan, {safe}_s0",
            f"color marine, {safe}_s1",
            f"color orange, {safe}_s2",
            f"color tv_blue, {safe}_top1",
            f"color violet, {safe}_top2",
            f"color salmon, {safe}_top3",
            f"disable {safe}_*",
        ])
    lines.extend(["# Enable one case with: enable <case_id>_*", "zoom"])
    VISUAL_DIR.mkdir(parents=True, exist_ok=True)
    (VISUAL_DIR / "load_archived_subset.pml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    manifest = read_json(MANIFEST_PATH)
    audit = Audit()
    results = result_map(audit)
    case_dirs = archive_map(audit)

    expected_case_ids = {case["case_id"] for case in manifest["cases"]}
    audit.check("archive_case_count", set(case_dirs) == expected_case_ids, "Every and only selected validation case has an archive.", expected=sorted(expected_case_ids), observed=sorted(case_dirs))
    audit.check("case_study_count", sum(bool(case.get("case_study")) for case in manifest["cases"]) == manifest["case_study_case_count"], "Manifest contains all six manuscript case-study targets.")

    observed_site_count = 0
    observed_ok_count = 0
    observed_no_coords_count = 0

    for case in manifest["cases"]:
        case_id = case["case_id"]
        primary_type = case["primary_type"]
        rows = matching_result_rows(case, results.get(primary_type, []))
        positions = case["expected_ok_positions"] + case["expected_no_coords_positions"]
        rows_by_position = {int(float(row.get("target_position") or row.get("position"))): row for row in rows}
        audit.check(f"result_rows:{case_id}", set(rows_by_position) == set(positions), "Result TSV contains exactly the expected target positions for this case.", expected=positions, observed=sorted(rows_by_position))
        observed_site_count += len(rows_by_position)

        for position in case["expected_ok_positions"]:
            row = rows_by_position.get(position, {})
            status = str(row.get("status", "")).strip().upper()
            restype = str(row.get("restype", "")).strip().upper()
            audit.check(f"status_ok:{case_id}:{position}", status == "OK", "Site has a successful benchmark result.", status=status)
            audit.check(f"restype:{case_id}:{position}", restype == RESIDUE_LETTER[primary_type], "Reported target residue type matches the benchmark group.", restype=restype, expected=RESIDUE_LETTER[primary_type])
            observed_ok_count += status == "OK"

        for position in case["expected_no_coords_positions"]:
            row = rows_by_position.get(position, {})
            status = str(row.get("status", "")).strip().upper()
            audit.check(f"status_no_coords:{case_id}:{position}", status == "NO_COORDS", "Missing experimental phosphate coordinates are reported explicitly without discarding valid neighboring sites.", status=status)
            observed_no_coords_count += status == "NO_COORDS"

        case_dir = case_dirs.get(case_id)
        audit.check(f"archive_present:{case_id}", case_dir is not None, "Per-case archive directory exists.")
        if case_dir is None:
            continue

        required_stage_paths = [
            case_dir / "input" / "reference_phospho.cif",
            case_dir / "input" / "stripped_input.cif",
            case_dir / "stage0" / "structure_s0.cif",
            case_dir / "stage1" / "structure_s1.cif",
            case_dir / "stage2" / "structure_s2.cif",
            case_dir / "stage3" / "structure_s3.cif",
        ]
        audit.check(f"stage_files:{case_id}", all(path.exists() for path in required_stage_paths), "Reference, stripped, Stage 0, Stage 1, Stage 2, and Stage 3 structures are archived.", missing=[str(path) for path in required_stage_paths if not path.exists()])

        stripped_path = case_dir / "input" / "stripped_input.cif"
        if stripped_path.exists():
            stripped = parse_structure(stripped_path)
            chain = case["match"]["CHAINID_1"]
            phosphate_absent = True
            for position in positions:
                residue = find_residue(stripped, chain, position, prefer_phosphate=True)
                phosphate_absent = phosphate_absent and (residue is None or "P" not in residue)
            audit.check(f"stripped:{case_id}", phosphate_absent, "Experimental phosphate atoms are absent from every stripped target residue.")

        check_stage0_geometry(audit, case, case_dir)
        check_scan_archive(audit, case, case_dir)
        check_rank_archives(audit, case, case_dir, rows_by_position)

        for position in positions:
            row = rows_by_position.get(position, {})
            audit.site_records.append({
                "case_id": case_id,
                "case_study": bool(case.get("case_study")),
                "primary_type": primary_type,
                "pdb_id": case["match"]["PDBID_1"],
                "chain": case["match"]["CHAINID_1"],
                "position": position,
                "status": row.get("status", "MISSING"),
                "rank1_postmin_rmsd": number(row.get("rank1_postmin_rmsd")),
                "rank2_postmin_rmsd": number(row.get("rank2_postmin_rmsd")),
                "rank3_postmin_rmsd": number(row.get("rank3_postmin_rmsd")),
                "top3_best_postmin_rmsd": number(row.get("top3_best_postmin_rmsd")),
                "archive_case_dir": str(case_dir),
            })

    audit.check("site_count", observed_site_count == manifest["expected_site_count"], "Observed result-site count matches the manifest.", expected=manifest["expected_site_count"], observed=observed_site_count)
    audit.check("ok_count", observed_ok_count == manifest["expected_ok_site_count"], "Successful site count matches the manifest.", expected=manifest["expected_ok_site_count"], observed=observed_ok_count)
    audit.check("no_coords_count", observed_no_coords_count == manifest["expected_no_coords_site_count"], "NO_COORDS site count matches the manifest.", expected=manifest["expected_no_coords_site_count"], observed=observed_no_coords_count)

    report = {
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(ROOT),
        "summary": {
            "checks": len(audit.checks),
            "passed": len(audit.checks) - len(audit.failures),
            "failed": len(audit.failures),
            "cases": len(manifest["cases"]),
            "sites": observed_site_count,
            "ok_sites": observed_ok_count,
            "no_coords_sites": observed_no_coords_count,
            "case_study_cases": sum(bool(case.get("case_study")) for case in manifest["cases"]),
        },
        "configured_priors": PRIORS,
        "checks": audit.checks,
        "sites": audit.site_records,
    }
    VERIFICATION_DIR.mkdir(parents=True, exist_ok=True)
    report_path = VERIFICATION_DIR / "spotcheck_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    write_pymol_loader(case_dirs, manifest)

    summary = report["summary"]
    print(json.dumps(summary, indent=2))
    print(f"Detailed report: {report_path}")
    print(f"PyMOL loader: {VISUAL_DIR / 'load_archived_subset.pml'}")
    if audit.failures:
        print("\nFailed checks:")
        for failure in audit.failures:
            print(f"- {failure['check_id']}: {failure['message']}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
