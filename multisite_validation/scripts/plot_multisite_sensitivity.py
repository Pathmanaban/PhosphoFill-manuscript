#!/usr/bin/env python3
"""Rebuild the multi-site sensitivity figure from the frozen source-value TSV."""

from __future__ import annotations

from analyze_multisite_controls import ANALYSIS_DIR, make_plot, read_tsv


def main() -> None:
    comparison = read_tsv(ANALYSIS_DIR / "multisite_per_site_comparison.tsv")
    for row in comparison:
        for key, value in list(row.items()):
            if (
                key == "nearest_phosphate_distance_A"
                or key.endswith("_top1_rmsd")
                or key.endswith("_top3_best_rmsd")
                or "_delta_top1_rmsd_" in key
                or "_delta_top3_best_rmsd_" in key
            ):
                row[key] = None if str(value).strip() == "" else float(value)
    make_plot(comparison, [])


if __name__ == "__main__":
    main()

