#!/usr/bin/env python3
"""Rebuild Figure 8 panels A--C from the corrected archived benchmark poses.

The comparator structures and the manually prepared PyMOL panels D--G are
unchanged.  PhosphoFill top-1 and best-of-three values, torsions, and contact
recovery are recalculated from the independently minimised rank archives made
by validation_structure_subset.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-root", required=True)
    parser.add_argument("--subset-root", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def corrected_phosphofill_rows(subset_root: Path) -> pd.DataFrame:
    frames = []
    for restype in ("TPO", "PTR"):
        path = subset_root / "outputs" / "results" / f"{restype}_validation_subset.tsv"
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        frame = frame[frame["status"].str.upper() == "OK"].copy()
        frames.append(frame)
    data = pd.concat(frames, ignore_index=True)

    # The archived paths are stored relative to validation_structure_subset.
    archive_columns = [
        "archive_stage3_structure",
        "archive_stage3_top1_structure",
        "archive_stage3_top2_structure",
        "archive_stage3_top3_structure",
    ]
    for column in archive_columns:
        data[column] = data[column].map(
            lambda value: str((subset_root / value).resolve()) if str(value).strip() else ""
        )

    # Legacy field names consumed by the original case-study analysis.
    data["stage3_post_minimization_rmsd"] = data["rank1_postmin_rmsd"]
    data["top3_best_rmsd"] = data["top3_best_postmin_rmsd"]
    return data


def main() -> None:
    args = parse_args()
    benchmark_root = Path(args.benchmark_root).resolve()
    subset_root = Path(args.subset_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(benchmark_root))
    import case_study_tool_analysis as base

    corrected = corrected_phosphofill_rows(subset_root)
    corrected_path = output_dir / "corrected_phosphofill_case_studies.tsv"
    corrected.to_csv(corrected_path, sep="\t", index=False)

    original_read_tsv = base.read_tsv

    def read_tsv(path):
        normalized = str(path).replace("\\", "/")
        if normalized in {
            "TPO_final_res.tsv",
            "case_study_inputs/PTR_P28482_ERK2_Y187.tsv",
        }:
            return corrected.copy()
        return original_read_tsv(path)

    base.read_tsv = read_tsv
    base.OUT_DIR = output_dir
    base.FIG_DIR = output_dir / "figures"
    for case in base.CASES:
        case["pf_file"] = str(corrected_path)
    base.main()

    import case_study_environment_analysis as environment

    environment.OUT_DIR = output_dir
    environment.FIG_DIR = output_dir / "figures"
    environment.PYMOL_DIR = output_dir / "pymol"
    environment.main()

    import make_cdk2_manuscript_figure_readable as figure

    figure.OUT_DIR = output_dir
    figure.FIG_DIR = output_dir / "figures"
    figure.OUTPUT_PNG = figure.FIG_DIR / "fig8_case_study_corrected.png"
    figure.OUTPUT_TIF = figure.FIG_DIR / "fig8_case_study_corrected.tif"
    figure.OUTPUT_PDF = figure.FIG_DIR / "fig8_case_study_corrected.pdf"
    figure.main()


if __name__ == "__main__":
    main()
