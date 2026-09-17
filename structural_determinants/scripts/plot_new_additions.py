#!/usr/bin/env python3
"""Review figure for paired anchor distances and multi-site benchmark context."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


COLORS = {"TPO": "#1565c0", "SEP": "#2e7d32", "PTR": "#6a1b9a"}
TYPES = ("TPO", "SEP", "PTR")


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main() -> None:
    parent = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path,
                        default=parent / "outputs" / "new_additions")
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[4] / "output" / "pdf"
                        / "structural_determinants_additions_review.pdf")
    args = parser.parse_args()
    paired = read_tsv(args.analysis_dir / "paired_anchor_basic_sites.tsv")
    multi = read_tsv(args.analysis_dir / "multisite_joined_environment.tsv")
    summary = json.loads((args.analysis_dir / "paired_and_multisite_summary.json")
                         .read_text(encoding="utf-8"))
    rng = np.random.default_rng(20260914)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42})
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.2, 5.2),
                                     gridspec_kw={"width_ratios": [1, 1.18]})
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.255, top=0.85,
                        wspace=0.27)

    left.axvline(0, color="#333333", lw=1, ls="--", zorder=0)
    for index, residue in enumerate(TYPES):
        rows = [row for row in paired if row["comparison"].startswith(residue)]
        values = np.array([float(row["difference_A"]) for row in rows])
        ypos = 2 - index
        left.scatter(values, ypos + rng.uniform(-0.20, 0.20, size=len(values)),
                     s=15, color=COLORS[residue], alpha=0.38, linewidths=0,
                     zorder=2)
        stats = summary["paired"]["by_residue"][residue]
        median = stats["paired_difference_median_A"]
        lo, hi = stats["paired_difference_95pct_protein_cluster_bootstrap_ci_A"]
        left.plot([lo, hi], [ypos, ypos], color="#111111", lw=2.5, zorder=3)
        left.scatter([median], [ypos], color="#111111", marker="D", s=36,
                     zorder=4)
    left.set_yticks([2, 1, 0],
                    [f"{r} (n={summary['paired']['by_residue'][r]['complete_anchor_distance_pairs_n']})"
                     for r in TYPES])
    left.set_xlim(-9.0, 7.1)
    left.set_ylim(-0.55, 2.55)
    left.set_xlabel("Matched-site modified - unmodified distance (Å)")
    left.set_title("a  Same bridging-oxygen marker", loc="left", weight="bold")
    left.grid(axis="x", color="#e5e5e5", lw=0.65)
    left.set_axisbelow(True)

    positions = []
    data = []
    fills = []
    labels = []
    for index, residue in enumerate(TYPES):
        center = 1 + index * 2.6
        for offset, group in ((-0.34, "salt_like"), (0.34, "other")):
            rows = [row for row in multi if row["restype_3"] == residue
                    and row["contact_group"] == group]
            vals = np.array([float(row["rank1_postmin_rmsd"]) for row in rows])
            pos = center + offset
            positions.append(pos)
            data.append(vals)
            fills.append(COLORS[residue] if group == "salt_like" else "#ffffff")
            labels.append(f"{group}\n{len(vals)}")
            right.scatter(pos + rng.uniform(-0.12, 0.12, len(vals)), vals,
                          s=8, color=COLORS[residue],
                          alpha=0.18 if group == "salt_like" else 0.25,
                          linewidths=0, zorder=1)
    boxes = right.boxplot(data, positions=positions, widths=0.42,
                          showfliers=False, patch_artist=True,
                          medianprops={"color": "#111111", "linewidth": 1.7},
                          whiskerprops={"color": "#555555"},
                          capprops={"color": "#555555"},
                          boxprops={"edgecolor": "#444444", "linewidth": 0.9})
    for box, fill in zip(boxes["boxes"], fills):
        box.set_facecolor(fill)
        box.set_alpha(0.66 if fill != "#ffffff" else 1)
    right.axhline(1.0, color="#666666", lw=0.85, ls="--", zorder=0)
    right.set_xlim(0.1, 7.1)
    right.set_ylim(0, 4.0)
    right.set_ylabel("PhosphoFill top-1 phosphate RMSD (Å)")
    right.set_xticks([1, 3.6, 6.2], TYPES)
    right.set_title("b  Multi-site experimental contact classes", loc="left",
                    weight="bold")
    right.grid(axis="y", color="#e5e5e5", lw=0.65)
    right.set_axisbelow(True)
    right.legend(handles=[Patch(facecolor="#777777", edgecolor="#444444",
                                label="Salt-bridge-like"),
                          Patch(facecolor="white", edgecolor="#444444",
                                label="Other audit classes")],
                 loc="upper right", frameon=False, fontsize=8)
    for pos, label in zip(positions, labels):
        right.text(pos, -0.24, f"n={label.split(chr(10))[-1]}", ha="center", va="top",
                   fontsize=7.5, color="#555555", clip_on=False)

    fig.text(0.075, 0.105,
             "a: one median per protein-position and state; black diamond and bar = median paired difference and 95% protein-cluster bootstrap CI.",
             fontsize=7.7, color="#3b3b3b")
    fig.text(0.075, 0.070,
             "b: 490/494 multi-site benchmark records matched by accession, PDB, chain, position, and residue type; contacts are defined in the crystal reference.",
             fontsize=7.7, color="#3b3b3b")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
